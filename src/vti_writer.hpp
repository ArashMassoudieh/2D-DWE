// =============================================================================
// vti_writer.hpp -- VTK ImageData (.vti) output, written directly as XML.
//
// WHY NOT LINK VTK FOR THIS
//   A structured grid maps exactly onto vtkImageData, and the .vti container is
//   plain XML with ASCII arrays. Writing it by hand keeps the solver free of any
//   VTK link dependency, so the core builds and runs anywhere, and the files
//   open natively in ParaView and in any VTK reader.
//
//   Nothing in this project links VTK. Results are opened in ParaView.
//
// WHY .vti AND NOT .vtp
//   .vtp is PolyData -- an explicit list of points and connectivity. This model
//   is a uniform structured grid, so ImageData stores no geometry at all beyond
//   an origin and a spacing, which is both far smaller on disk and the form
//   every ParaView filter (Glyph, Contour, Threshold, Slice) handles fastest.
//   Writing the same field as PolyData would mean emitting 9.2 M explicit point
//   coordinates per output to describe a grid already defined by three numbers.
//   ParaView reads the vectors below and draws them as arrows either way.
//
// ARRAYS WRITTEN
//   depth       h                 [m]   water depth
//   bed         z                 [m]   terrain
//   head        z + h             [m]   piezometric head = water surface
//   total_head  z + h + |u|^2/2g  [m]   energy head
//   velocity    (ux, uy, 0)     [m/s]   VECTOR -- glyph this directly
//   speed       |u|              [m/s]   for colouring and thresholds
//
// ORIENTATION
//   Row 0 of the arrays is the NORTH row (grid.hpp). VTK ImageData indexes +y
//   upward, so rows are emitted bottom-up and the origin is the SOUTH-WEST
//   corner. Get this wrong and every result is mirrored.
// =============================================================================
#pragma once

#include <cstdio>
#include <string>
#include <vector>

#include "grid.hpp"
#include "fields.hpp"

#include <cmath>

namespace dwe {

inline bool write_vti(const std::string& path,
                      const Fields& F,
                      double h_min = 1.0e-3)
{
    const Grid& g = F.g;
    std::FILE* fp = std::fopen(path.c_str(), "w");
    if (!fp) return false;

    const double oy = g.y0 - g.ny * g.dx;   // south edge

    std::fprintf(fp,
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
        "  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"%.6f %.6f 0\" "
        "Spacing=\"%.6f %.6f 1\">\n"
        "    <Piece Extent=\"0 %d 0 %d 0 0\">\n"
        "      <CellData Scalars=\"depth\" Vectors=\"velocity\">\n",
        g.nx, g.ny, g.x0, oy, g.dx, g.dx, g.nx, g.ny);

    auto emit = [&](const char* name, const std::vector<double>& a) {
        std::fprintf(fp, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n",
                     name);
        for (int j = g.ny - 1; j >= 0; --j) {          // bottom-up for VTK
            for (int i = 0; i < g.nx; ++i)
                std::fprintf(fp, "%.4f ", float(a[g.k(i, j)]));
            std::fputc('\n', fp);
        }
        std::fprintf(fp, "        </DataArray>\n");
    };

    // ---- cell-centred velocity ------------------------------------------
    // q lives on faces; average the two faces bounding each cell and divide by
    // depth. Below h_min the division is meaningless -- a vanishing depth makes
    // u = q/h blow up -- so dry and near-dry cells are reported as zero rather
    // than as enormous spurious arrows that would wreck any glyph scaling.
    //
    // SIGN: qy is positive SOUTHWARD (grid.hpp -- fy(i,j) is the north edge and
    // its left neighbour is the north cell). VTK's +y is north because the rows
    // below are emitted bottom-up, so the y component is NEGATED on the way out.
    std::vector<double> vx(g.ncell(), 0.0), vy(g.ncell(), 0.0),
                        sp(g.ncell(), 0.0), th(g.ncell(), 0.0), hd(g.ncell(), 0.0);
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            const std::size_t k = g.k(i, j);
            hd[k] = F.H(k);
            if (F.h[k] <= h_min) { th[k] = hd[k]; continue; }
            const double inv = 1.0 / F.h[k];
            vx[k] =  0.5 * (F.qx[g.fx(i, j)] + F.qx[g.fx(i + 1, j)]) * inv;
            vy[k] = -0.5 * (F.qy[g.fy(i, j)] + F.qy[g.fy(i, j + 1)]) * inv;
            sp[k] = std::sqrt(vx[k] * vx[k] + vy[k] * vy[k]);
            th[k] = hd[k] + sp[k] * sp[k] / (2.0 * 9.80665);
        }
    }

    emit("depth",      F.h);
    emit("bed",        F.z);
    emit("head",       hd);
    emit("total_head", th);
    emit("speed",      sp);

    // ---- the vector, 3 components interleaved ----------------------------
    std::fprintf(fp, "        <DataArray type=\"Float32\" Name=\"velocity\" "
                     "NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int j = g.ny - 1; j >= 0; --j) {
        for (int i = 0; i < g.nx; ++i) {
            const std::size_t k = g.k(i, j);
            std::fprintf(fp, "%.4f %.4f 0 ", float(vx[k]), float(vy[k]));
        }
        std::fputc('\n', fp);
    }
    std::fprintf(fp, "        </DataArray>\n");

    std::fprintf(fp,
        "      </CellData>\n"
        "    </Piece>\n"
        "  </ImageData>\n"
        "</VTKFile>\n");

    std::fclose(fp);
    return true;
}

} // namespace dwe
