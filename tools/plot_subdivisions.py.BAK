#!/usr/bin/env python3
"""
Maximum-depth map clipped to the NHA subdivisions, with structures classified by
the depth of water on the ground at their location.

    python3 tools/plot_subdivisions.py out_c6_1m 1

Also prints the per-subdivision statistics that appear in the report table.
Depths are GROUND depths: structures are absent from the terrain model and no
finished-floor elevations exist, so a structure in a depth class is not a
flooded building.
"""
import sys
RUN = sys.argv[1] if len(sys.argv) > 1 else "out_c6_1m"
DXM = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
TITLE = sys.argv[4] if len(sys.argv) > 4 else (
    "C10" if "nm141" in RUN else "C6 Coyote Canyon")
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, LightSource
from matplotlib.patches import Patch
from osgeo import gdal, ogr
gdal.UseExceptions(); ogr.UseExceptions()
FT=3.280839895
OUT="/home/arash/Projects/2D-DWE/"+RUN
r=gdal.Open(OUT+"/depth_max.tif"); b=r.GetRasterBand(1)
a=b.ReadAsArray().astype(float); nd=b.GetNoDataValue()
a=np.where((a==nd)|~np.isfinite(a),0.0,a)*FT
gt=r.GetGeoTransform(); ny,nx=a.shape
ext=[gt[0],gt[0]+nx*gt[1],gt[3]-ny*gt[1],gt[3]]
DEMPATH = sys.argv[3] if len(sys.argv) > 3 else \
    "/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs_2D/cliped_dem.tif"
t=gdal.Warp("",DEMPATH,
            format="MEM",width=nx,height=ny,resampleAlg="average")
zb=t.GetRasterBand(1).ReadAsArray().astype(float)

_s=ogr.Open("/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs/05_Subdivisions_clip.gpkg")
sly=_s.GetLayer()
polys=[]
for f in sly:
    g=f.GetGeometryRef()
    if g is None: continue
    e=g.GetEnvelope()
    if e[1]<ext[0] or e[0]>ext[1] or e[3]<ext[2] or e[2]>ext[3]: continue
    polys.append((f.GetField("proj_num"),g.Clone()))
_p=ogr.Open("/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs/00_NHA_Structures_clip.gpkg")
ply=_p.GetLayer()
sx,sy,sd=[],[],[]
for f in ply:
    g=f.GetGeometryRef()
    if g is None: continue
    x,y=g.GetX(),g.GetY()
    i=int((x-ext[0])/gt[1]); j=int((ext[3]-y)/gt[1])
    if 0<=i<nx and 0<=j<ny: sx.append(x); sy.append(y); sd.append(a[j,i])
sx,sy,sd=np.array(sx),np.array(sy),np.array(sd)

xs=[p[1].GetEnvelope() for p in polys]
bx0=min(e[0] for e in xs)-260; bx1=max(e[1] for e in xs)+260
by0=min(e[2] for e in xs)-260; by1=max(e[3] for e in xs)+260

fig,ax=plt.subplots(figsize=(14,10.5),dpi=150)
ls=LightSource(azdeg=315,altdeg=45)
ax.imshow(ls.hillshade(zb,vert_exag=3,dx=DXM,dy=DXM),cmap="gray",extent=ext,alpha=.9,zorder=1)
lv=[0.5,1,2,3,4,6]
cm=ListedColormap(["#bfe3f7","#7fc4e8","#3d94cf","#1c5fa8","#0d2f6e"])
im=ax.imshow(np.ma.masked_less(a,0.5),cmap=cm,norm=BoundaryNorm(lv,cm.N),extent=ext,alpha=.85,zorder=2)

_LBL=[]
for pn,g in polys:
    for k in range(g.GetGeometryCount()):
        gg=g.GetGeometryRef(k)
        ring=gg.GetGeometryRef(0) if gg.GetGeometryCount() else gg
        p=np.array(ring.GetPoints())
        if p.size:
            ax.plot(p[:,0],p[:,1],c="#111",lw=2.4,zorder=6)
            ax.plot(p[:,0],p[:,1],c="#ffd400",lw=1.2,zorder=6)
    # Labels were plain text sitting above each envelope with nothing tying them
    # to their polygon; where parcels are adjacent -- which is the normal case --
    # they land in a row over the wrong ones. Anchor each to its polygon's
    # representative point and stagger the callouts to the side with a leader.
    e=g.GetEnvelope()
    try:
        rp=g.PointOnSurface(); ax_, ay_ = rp.GetX(), rp.GetY()
    except Exception:
        ax_, ay_ = (e[0]+e[1])/2, (e[2]+e[3])/2
    # Put each label at ITS OWN parcel's height in the left margin, so leaders run
    # roughly horizontally and do not cross one another. Stacking them all in one
    # corner made three leaders sweep across the whole map and intersect.
    lx = bx0 + 0.11*(bx1-bx0)
    ly = ay_
    span = by1 - by0
    for prev in _LBL:
        if abs(ly - prev) < 0.062*span:
            ly = prev - 0.062*span
    _LBL.append(ly)
    ax.annotate(str(pn), xy=(ax_, ay_), xytext=(lx, ly),
                ha="center", va="center", fontsize=10.5, weight="bold", zorder=11,
                arrowprops=dict(arrowstyle="-|>", lw=1.3, color="k",
                                shrinkA=0, shrinkB=2),
                bbox=dict(fc="#ffd400", ec="k", lw=.8, pad=2.2))

dry=sd<=0.5
ax.scatter(sx[dry],sy[dry],s=26,c="#ffffff",ec="#333",lw=.7,zorder=7,label="Structure, below 0.5 ft (%d)" % int(dry.sum()))
ax.scatter(sx[~dry],sy[~dry],s=90,c="#e8112d",ec="k",lw=1.1,marker="s",zorder=8,
           label="Structure, above 0.5 ft (%d)" % int((~dry).sum()))
cb=fig.colorbar(im,ax=ax,shrink=.85,pad=.012,ticks=lv)
cb.set_label("Maximum depth over the 9-hour event (ft)",fontsize=10)
ax.set_xlim(bx0,bx1); ax.set_ylim(by0,by1); ax.set_aspect("equal")
ax.set_title(TITLE + " — flooding WITHIN the subdivisions\n"
             "2D-DWE, %.0f m grid, 100-yr 6-hr transient hydrograph — PRELIMINARY" % DXM,fontsize=12.5)
ax.set_xlabel("UTM Zone 12N Easting (m)"); ax.set_ylabel("Northing (m)")
h,l=ax.get_legend_handles_labels()
h.append(Patch(fc="none",ec="#ffd400",lw=2,label="NHA subdivision boundary"))
ax.legend(handles=h,loc="lower left",fontsize=9.5,framealpha=.93)
bl=bx0+90
ax.plot([bl,bl+304.8],[by0+70,by0+70],c="k",lw=4,zorder=9,solid_capstyle="butt")
ax.text(bl+152,by0+95,"1,000 ft",ha="center",fontsize=9,zorder=9)
plt.tight_layout(); plt.savefig(OUT+"/c6_subdivision_flooding.png",bbox_inches="tight")
print("wrote map")


# ---- statistics for the report table -------------------------------------
from osgeo import gdal as _g
acre = gt[1]*gt[1]/4046.8564224
print("%-12s %8s %10s %8s %9s" % ("subdivision","acres","ac >1 ft","% area","max ft"))
for pn, g in polys:
    mem = ogr.GetDriverByName("Memory").CreateDataSource("m")
    l2 = mem.CreateLayer("p", sly.GetSpatialRef(), ogr.wkbPolygon)
    ft_ = ogr.Feature(l2.GetLayerDefn()); ft_.SetGeometry(g); l2.CreateFeature(ft_)
    msk = _g.GetDriverByName("MEM").Create("", nx, ny, 1, _g.GDT_Byte)
    msk.SetGeoTransform(gt); msk.SetProjection(r.GetProjection())
    _g.RasterizeLayer(msk, [1], l2, burn_values=[1])
    m = msk.GetRasterBand(1).ReadAsArray().astype(bool)
    if not m.sum(): continue
    d = a[m]
    print("%-12s %8.0f %10.1f %7.1f%% %9.2f"
          % (pn, m.sum()*acre, (d>1).sum()*acre, 100*(d>1).mean(), d.max()))
print("structures: %d in domain, %d >0.5 ft, %d >1 ft, deepest %.2f ft"
      % (len(sd), int((sd>0.5).sum()), int((sd>1).sum()), sd.max()))
