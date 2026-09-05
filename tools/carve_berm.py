#!/usr/bin/env python3
"""
Carve the proposed C6 boundary berm into the 1 m terrain, for re-running the
2-D model with the berm in place.

    python3 tools/carve_berm.py [top_width_ft] [side_slope_z]      defaults 6, 3

Geometry (from the berm design, see Report_for_C6/Figures/NM15-097_berm_profile.png):
  * centerline  outputs_2D/berm_centerline.shp  -- 2 segments, 3,117 ft
  * Station 0-315 ft: FLAT crest at 6153.4 ft -- the mean crest of the existing
    south-west berm.  Carving here brings that berm up to a uniform designed
    section and closes the low spot at its south-east terminus.
  * Station 315-3,115 ft: crest on a single constant grade
        crest(sta) = 6153.4 - 0.00200 * (sta - 315)      [ft, NAVD88]
    continuous with the flat reach, no step.
  Terrain is only raised, so where the existing berm already exceeds the design
  crest nothing changes.
  * Trapezoidal: TOP_W ft crest, side slopes z H : 1 V both faces.

The terrain is only ever RAISED: new = max(old, berm surface).  Where the berm
surface falls below existing ground the side slope has met grade and nothing
changes, so no explicit toe clipping is needed.

Writes
  outputs_2D/cliped_dem_berm.tif    terrain with the berm, model-ready
  outputs_2D/berm_raise.tif         raise above existing ground, for QC
"""
import sys
import numpy as np
from osgeo import gdal, ogr
from scipy.spatial import cKDTree

gdal.UseExceptions(); ogr.UseExceptions()

FT      = 3.280839895
TOP_W   = float(sys.argv[1]) if len(sys.argv) > 1 else 6.0     # ft
SIDE_Z  = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0     # H:V
SITE    = "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097"
DEM_IN  = SITE + "/outputs_2D/cliped_dem.tif"
DEM_OUT = SITE + "/outputs_2D/cliped_dem_berm.tif"
RAISE   = SITE + "/outputs_2D/berm_raise.tif"
CL      = SITE + "/outputs_2D/berm_centerline.shp"

BUILD_FROM_FT = 0.0          # carve the whole alignment
GRADE         = -0.00200     # ft per ft
CREST_AT_315  = 6153.4       # ft NAVD88 at station 315
STA_REF       = 315.0        # flat at CREST_AT_315 upstream of this station

ds  = gdal.Open(DEM_IN)
gt  = ds.GetGeoTransform(); prj = ds.GetProjection()
NX, NY = ds.RasterXSize, ds.RasterYSize
Z   = ds.GetRasterBand(1).ReadAsArray().astype("f8")          # metres
nd  = ds.GetRasterBand(1).GetNoDataValue()

# ---- centerline, densified, station in feet -------------------------------
src = ogr.Open(CL); lyr = src.GetLayer(0)
feat = lyr.GetNextFeature(); geom = feat.GetGeometryRef().Clone()
P = np.array([geom.GetPoint(i)[:2] for i in range(geom.GetPointCount())])
seg = np.hypot(*np.diff(P, axis=0).T); s = np.r_[0, np.cumsum(seg)]
tt = np.arange(0, s[-1], 0.5)                                  # 0.5 m spacing
cx = np.interp(tt, s, P[:, 0]); cy = np.interp(tt, s, P[:, 1])
sta_ft = tt * FT
keep = sta_ft >= BUILD_FROM_FT
cx, cy, sta_ft = cx[keep], cy[keep], sta_ft[keep]
crest_ft = np.where(sta_ft <= STA_REF,
                    CREST_AT_315,                                   # flat reach
                    CREST_AT_315 + GRADE * (sta_ft - STA_REF))      # constant grade
print("berm: station %.0f to %.0f ft" % (sta_ft[0], sta_ft[-1]))
print("      flat at %.2f ft to station %.0f, then %.3f%% down to %.2f ft"
      % (CREST_AT_315, STA_REF, 100 * GRADE, crest_ft[-1]))
print("section: %.0f ft top width, %.0f:1 side slopes" % (TOP_W, SIDE_Z))

# ---- cells within reach of the alignment ----------------------------------
maxh_ft   = float(np.max(crest_ft - Z[
    ((cy - gt[3]) / gt[5]).astype(int), ((cx - gt[0]) / gt[1]).astype(int)] * FT))
half_w_ft = TOP_W / 2 + SIDE_Z * max(maxh_ft, 0.0) + 5.0
half_w_m  = half_w_ft / FT
print("search half-width %.1f ft (max height %.1f ft)" % (half_w_ft, maxh_ft))

i0 = max(int((cx.min() - half_w_m - gt[0]) / gt[1]) - 2, 0)
i1 = min(int((cx.max() + half_w_m - gt[0]) / gt[1]) + 2, NX)
j0 = max(int((cy.max() + half_w_m - gt[3]) / gt[5]) - 2, 0)
j1 = min(int((cy.min() - half_w_m - gt[3]) / gt[5]) + 2, NY)
jj, ii = np.mgrid[j0:j1, i0:i1]
px = gt[0] + (ii + 0.5) * gt[1]
py = gt[3] + (jj + 0.5) * gt[5]

tree = cKDTree(np.c_[cx, cy])
dist_m, idx = tree.query(np.c_[px.ravel(), py.ravel()], k=1)
dist_ft = dist_m * FT
crest_here = crest_ft[idx]

berm_ft = np.where(dist_ft <= TOP_W / 2,
                   crest_here,
                   crest_here - (dist_ft - TOP_W / 2) / SIDE_Z)
berm_m = (berm_ft / FT).reshape(px.shape)

sub = Z[j0:j1, i0:i1]
new = np.maximum(sub, berm_m)
if nd is not None:
    new = np.where(sub == nd, sub, new)
raise_m = new - sub

Znew = Z.copy(); Znew[j0:j1, i0:i1] = new
Rall = np.zeros_like(Z); Rall[j0:j1, i0:i1] = raise_m

cell_area_ft2 = (gt[1] * FT) ** 2
vol_cy = (raise_m * FT * cell_area_ft2).sum() / 27.0
print("\nfootprint  %.2f ac  (%d cells raised)"
      % ((raise_m > 0.003).sum() * cell_area_ft2 / 43560.0, (raise_m > 0.003).sum()))
print("max raise  %.2f ft" % (raise_m.max() * FT))
print("fill       %s cy" % format(round(vol_cy), ","))

drv = gdal.GetDriverByName("GTiff")
for arr, path, unit in ((Znew, DEM_OUT, "m"), (Rall * FT, RAISE, "ft")):
    o = drv.Create(path, NX, NY, 1, gdal.GDT_Float32,
                   options=["COMPRESS=DEFLATE", "TILED=YES"])
    o.SetGeoTransform(gt); o.SetProjection(prj)
    b = o.GetRasterBand(1); b.WriteArray(arr.astype("f4"))
    if unit == "m" and nd is not None:
        b.SetNoDataValue(nd)
    o.FlushCache(); o = None
    print("wrote", path)
