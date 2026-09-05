#ifdef DWE_WITH_GDAL

#include "io_gdal.hpp"

#include <gdal_priv.h>
#include <cpl_conv.h>
#include <ogr_spatialref.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <limits>

namespace dwe {

namespace {

struct GdalInit
{
    GdalInit() { GDALAllRegister(); }
};
static GdalInit g_gdal_init;

// Split a CSV line, honouring double-quoted fields. The NLCD class names in
// manning_n_2d.csv contain commas ("Developed, Open Space"), so a naive split on
// ',' silently shifts the n column into the name column and every roughness
// value comes back wrong.
std::vector<std::string> split_csv(const std::string& line)
{
    std::vector<std::string> out;
    std::string cur;
    bool q = false;
    for (char c : line) {
        if (c == '"')            q = !q;
        else if (c == ',' && !q) { out.push_back(cur); cur.clear(); }
        else                     cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

} // namespace

// -----------------------------------------------------------------------------
bool raster_meta(const std::string& path, RasterMeta& m)
{
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(path.c_str(), GA_ReadOnly));
    if (!ds) return false;
    double gt[6];
    if (ds->GetGeoTransform(gt) != CE_None) { GDALClose(ds); return false; }
    m.nx = ds->GetRasterXSize();
    m.ny = ds->GetRasterYSize();
    m.dx = gt[1];
    m.x0 = gt[0];
    m.y0 = gt[3];
    int has = 0;
    m.nodata = ds->GetRasterBand(1)->GetNoDataValue(&has);
    m.has_nodata = has != 0;
    const char* wkt = ds->GetProjectionRef();
    m.crs = wkt ? wkt : "";
    GDALClose(ds);
    return true;
}

// -----------------------------------------------------------------------------
bool load_terrain(const std::string& path,
                  double target_dx,
                  Agg agg,
                  Grid& g,
                  std::vector<double>& z,
                  std::vector<char>& valid,
                  std::string* err)
{
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(path.c_str(), GA_ReadOnly));
    if (!ds) { if (err) *err = "cannot open " + path; return false; }

    double gt[6];
    if (ds->GetGeoTransform(gt) != CE_None) {
        if (err) *err = "no geotransform in " + path;
        GDALClose(ds); return false;
    }
    const int    snx = ds->GetRasterXSize();
    const int    sny = ds->GetRasterYSize();
    const double sdx = gt[1];
    if (sdx <= 0.0) { if (err) *err = "non-positive pixel size"; GDALClose(ds); return false; }

    // Integer block factor. Refusing a non-integer ratio keeps the aggregation
    // exact -- a fractional factor would need interpolation, which reintroduces
    // the smoothing bias this function documents.
    int f = int(std::lround(target_dx / sdx));
    if (f < 1) f = 1;
    const double eff_dx = f * sdx;

    GDALRasterBand* band = ds->GetRasterBand(1);
    int has_nd = 0;
    const double nd = band->GetNoDataValue(&has_nd);

    const int nx = snx / f;
    const int ny = sny / f;
    if (nx < 2 || ny < 2) { if (err) *err = "target_dx too coarse"; GDALClose(ds); return false; }

    g = Grid(nx, ny, eff_dx, gt[0], gt[3]);
    z.assign(g.ncell(), 0.0);
    valid.assign(g.ncell(), 1);

    std::vector<float> buf(std::size_t(snx) * f);

    for (int j = 0; j < ny; ++j) {
        // read f source rows at a time
        if (band->RasterIO(GF_Read, 0, j * f, snx, f,
                           buf.data(), snx, f, GDT_Float32, 0, 0) != CE_None) {
            if (err) *err = "RasterIO failed";
            GDALClose(ds); return false;
        }
        for (int i = 0; i < nx; ++i) {
            double acc = 0.0;
            double vmin =  std::numeric_limits<double>::infinity();
            double vmax = -std::numeric_limits<double>::infinity();
            int cnt = 0;
            for (int jj = 0; jj < f; ++jj) {
                const float* row = buf.data() + std::size_t(jj) * snx;
                for (int ii = 0; ii < f; ++ii) {
                    const double v = row[i * f + ii];
                    if (has_nd && v == nd)             continue;
                    if (!std::isfinite(v))             continue;
                    acc += v; vmin = std::min(vmin, v); vmax = std::max(vmax, v); ++cnt;
                }
            }
            const std::size_t k = g.k(i, j);
            if (cnt == 0) { z[k] = 0.0; valid[k] = 0; continue; }
            switch (agg) {
            case Agg::Mean: z[k] = acc / cnt; break;
            case Agg::Min:  z[k] = vmin;      break;
            case Agg::Max:  z[k] = vmax;      break;
            }
        }
    }
    GDALClose(ds);
    return true;
}

// -----------------------------------------------------------------------------
bool load_manning_direct(const std::string& n_path,
                         const Grid& g,
                         std::vector<double>& n,
                         double n_default,
                         std::string* err)
{
    n.assign(g.ncell(), n_default);

    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(n_path.c_str(), GA_ReadOnly));
    if (!ds) { if (err) *err = "cannot open " + n_path; return false; }
    double gt[6];
    if (ds->GetGeoTransform(gt) != CE_None) {
        if (err) *err = "no geotransform in " + n_path;
        GDALClose(ds); return false;
    }
    const int lnx = ds->GetRasterXSize(), lny = ds->GetRasterYSize();
    std::vector<float> nn(std::size_t(lnx) * lny);
    if (ds->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, lnx, lny,
                                       nn.data(), lnx, lny, GDT_Float32, 0, 0) != CE_None) {
        if (err) *err = "RasterIO failed on " + n_path;
        GDALClose(ds); return false;
    }
    int has_nd = 0;
    const double nd = ds->GetRasterBand(1)->GetNoDataValue(&has_nd);
    GDALClose(ds);

    const double inv_dx = 1.0 / gt[1];
    const double inv_dy = 1.0 / gt[5];

    std::size_t outside = 0, nodata = 0, bad = 0;
    for (int j = 0; j < g.ny; ++j) {
        const double y = g.yc(j);
        for (int i = 0; i < g.nx; ++i) {
            const double x = g.xc(i);
            const int ci = int((x - gt[0]) * inv_dx);
            const int cj = int((y - gt[3]) * inv_dy);
            const std::size_t k = g.k(i, j);
            if (ci < 0 || ci >= lnx || cj < 0 || cj >= lny) { ++outside; continue; }
            const double v = double(nn[std::size_t(cj) * lnx + ci]);
            if (has_nd && v == nd) { ++nodata; continue; }
            // A roughness raster that is silently wrong is worse than one that
            // fails loudly: reject anything outside a physically sane band
            // rather than routing a flood over n = 0 or n = 5.
            if (!(v > 0.005 && v < 0.5)) { ++bad; continue; }
            n[k] = v;
        }
    }
    if (outside || nodata || bad)
        std::printf("  roughness: %zu cells outside the raster, %zu nodata, %zu out of range"
                    " -- left at n = %.3f\n", outside, nodata, bad, n_default);
    return true;
}

bool load_manning(const std::string& landcover_path,
                  const std::string& csv_path,
                  const Grid& g,
                  std::vector<double>& n,
                  double n_default,
                  std::string* err)
{
    n.assign(g.ncell(), n_default);

    // --- the lookup table ------------------------------------------------
    std::map<int, double> table;
    {
        std::ifstream f(csv_path);
        if (!f) { if (err) *err = "cannot open " + csv_path; return false; }
        std::string line;
        bool header = true;
        while (std::getline(f, line)) {
            if (line.empty()) continue;
            if (header) { header = false; continue; }
            const auto c = split_csv(line);
            if (c.size() < 3) continue;
            try { table[std::stoi(c[0])] = std::stod(c[2]); } catch (...) {}
        }
    }
    if (table.empty()) { if (err) *err = "no rows parsed from " + csv_path; return false; }

    // --- sample the raster by coordinate ---------------------------------
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(landcover_path.c_str(), GA_ReadOnly));
    if (!ds) { if (err) *err = "cannot open " + landcover_path; return false; }
    double gt[6];
    if (ds->GetGeoTransform(gt) != CE_None) {
        if (err) *err = "no geotransform in " + landcover_path;
        GDALClose(ds); return false;
    }
    const int lnx = ds->GetRasterXSize(), lny = ds->GetRasterYSize();
    std::vector<int32_t> lc(std::size_t(lnx) * lny);
    if (ds->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, lnx, lny,
                                       lc.data(), lnx, lny, GDT_Int32, 0, 0) != CE_None) {
        if (err) *err = "RasterIO failed on land cover";
        GDALClose(ds); return false;
    }
    int has_nd = 0;
    const double nd = ds->GetRasterBand(1)->GetNoDataValue(&has_nd);
    GDALClose(ds);

    const double inv_dx = 1.0 / gt[1];
    const double inv_dy = 1.0 / gt[5];      // negative

    std::size_t unmapped = 0, outside = 0;
    for (int j = 0; j < g.ny; ++j) {
        const double y = g.yc(j);
        for (int i = 0; i < g.nx; ++i) {
            const double x = g.xc(i);
            const int ci = int((x - gt[0]) * inv_dx);
            const int cj = int((y - gt[3]) * inv_dy);
            const std::size_t k = g.k(i, j);
            if (ci < 0 || ci >= lnx || cj < 0 || cj >= lny) { ++outside; continue; }
            const int code = lc[std::size_t(cj) * lnx + ci];
            if (has_nd && double(code) == nd) { ++unmapped; continue; }
            const auto it = table.find(code);
            if (it == table.end()) { ++unmapped; continue; }
            n[k] = it->second;
        }
    }
    if (unmapped || outside)
        std::printf("  land cover: %zu cells unmapped, %zu outside the raster"
                    " -- left at n = %.3f\n", unmapped, outside, n_default);
    return true;
}

// -----------------------------------------------------------------------------
bool write_gtiff(const std::string& path,
                 const Grid& g,
                 const std::vector<double>& a,
                 const std::string& crs,
                 double nodata)
{
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!drv) return false;
    char* opts[] = { const_cast<char*>("COMPRESS=DEFLATE"),
                     const_cast<char*>("TILED=YES"), nullptr };
    GDALDataset* ds = drv->Create(path.c_str(), g.nx, g.ny, 1, GDT_Float32, opts);
    if (!ds) return false;

    double gt[6] = { g.x0, g.dx, 0.0, g.y0, 0.0, -g.dx };
    ds->SetGeoTransform(gt);

    OGRSpatialReference srs;
    if (crs.rfind("EPSG:", 0) == 0) srs.importFromEPSG(std::atoi(crs.c_str() + 5));
    else                            srs.importFromWkt(crs.c_str());
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    if (wkt) { ds->SetProjection(wkt); CPLFree(wkt); }

    std::vector<float> row(g.nx);
    GDALRasterBand* b = ds->GetRasterBand(1);
    b->SetNoDataValue(nodata);
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) row[i] = float(a[g.k(i, j)]);
        b->RasterIO(GF_Write, 0, j, g.nx, 1, row.data(), g.nx, 1, GDT_Float32, 0, 0);
    }
    GDALClose(ds);
    return true;
}

} // namespace dwe

#endif // DWE_WITH_GDAL
