import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, LightSource
from osgeo import gdal, ogr
gdal.UseExceptions(); ogr.UseExceptions()

OUT="/home/arash/Projects/2D-DWE/out_c6_10m"
SHP="/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs"
DEM="/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097/outputs_2D/cliped_dem.tif"
FT=3.280839895

d=gdal.Open(f"{OUT}/depth_max.tif"); bnd=d.GetRasterBand(1)
a=bnd.ReadAsArray().astype(float); nd=bnd.GetNoDataValue()
a=np.where(a==nd,0.0,a)
gt=d.GetGeoTransform(); ny,nx=a.shape
ext=[gt[0],gt[0]+nx*gt[1],gt[3]+ny*gt[5],gt[3]]
t=gdal.Warp("",DEM,format="MEM",width=nx,height=ny,resampleAlg="average")
zb=t.GetRasterBand(1).ReadAsArray().astype(float)
wse=zb+a                       # peak water surface

def samp(arr,x,y):
    i=np.clip(((x-gt[0])/gt[1]).astype(int),0,nx-1)
    j=np.clip(((y-gt[3])/gt[5]).astype(int),0,ny-1)
    return arr[j,i]

# --- the barrier polyline, densified ------------------------------------
_ds1=ogr.Open(f"{SHP}/natural_barrier.shp"); ly=_ds1.GetLayer()
_f1=ly.GetNextFeature(); pts=np.array(_f1.GetGeometryRef().GetPoints())[:,:2]
seg=np.r_[0,np.cumsum(np.hypot(*np.diff(pts,axis=0).T))]
s=np.arange(0,seg[-1],5.0)
bx=np.interp(s,seg,pts[:,0]); by=np.interp(s,seg,pts[:,1])
zc=samp(zb,bx,by); wc=samp(wse,bx,by); hc=samp(a,bx,by)

_ds2=ogr.Open(f"{SHP}/divide_limits.shp"); dl=_ds2.GetLayer()
_fs=[f for f in dl]
dpts=np.array([f.GetGeometryRef().GetPoint()[:2] for f in _fs])

fig=plt.figure(figsize=(16,11.5),dpi=140)
gs=fig.add_gridspec(2,1,height_ratios=[2.05,1],hspace=.16)

ax=fig.add_subplot(gs[0])
ls=LightSource(azdeg=315,altdeg=45)
ax.imshow(ls.hillshade(zb,vert_exag=3,dx=10,dy=10),cmap="gray",extent=ext,alpha=.85,zorder=1)
lv=[0.5,1,2,3,5,8,12,25]
cm=ListedColormap(["#c6e5f5","#8ec9ec","#4ba3d8","#2272b5","#1a4f95","#0d2f6e","#5a1a6e"])
im=ax.imshow(np.ma.masked_less(a*FT,0.5),cmap=cm,norm=BoundaryNorm(lv,cm.N),
             extent=ext,alpha=.82,zorder=2)
ax.plot(pts[:,0],pts[:,1],c="#ff9500",lw=2.6,zorder=5,label="Natural barrier")
ax.scatter(dpts[:,0],dpts[:,1],s=130,c="#e8112d",marker="o",ec="k",lw=1.2,zorder=6,
           label="Controlling-divide limits")
ov=hc*FT>0.5
ax.scatter(bx[ov],by[ov],s=13,c="#e8112d",zorder=7,label="Barrier overtopped")
cb=fig.colorbar(im,ax=ax,shrink=.9,pad=.012,ticks=lv); cb.set_label("Maximum depth (ft)",fontsize=10)
ax.set_title("C6 Coyote Canyon — 2D-DWE local-inertial solver, 100-yr peak 32,628 cfs\n"
             "10 m grid, 1-hour run (NOT at steady state) — PRELIMINARY",fontsize=12)
ax.set_xlabel("UTM Zone 12N Easting (m)");ax.set_ylabel("Northing (m)")
ax.legend(loc="lower left",fontsize=9,framealpha=.92)
sx,sy=ext[0]+300,ext[2]+180
ax.plot([sx,sx+1609.34],[sy,sy],c="k",lw=4,zorder=8,solid_capstyle="butt")
ax.text(sx+804,sy+60,"1 mile",ha="center",fontsize=9,zorder=8)
ax.annotate("",xy=(ext[1]-250,ext[3]-200),xytext=(ext[1]-250,ext[3]-520),
            arrowprops=dict(arrowstyle="-|>",lw=2,color="k"),zorder=8)
ax.text(ext[1]-250,ext[3]-150,"N",ha="center",fontsize=11,weight="bold",zorder=8)
ax.set_aspect("equal")

ax2=fig.add_subplot(gs[1])
sf=s*FT
ax2.fill_between(sf,zc*FT,zc.min()*FT-8,color="#c8b89a",zorder=2,label="Barrier crest (10 m model grid, MEAN-aggregated)")
ax2.plot(sf,zc*FT,c="#6b5b3e",lw=1.6,zorder=3)
ax2.plot(sf,wc*FT,c="#1a4f95",lw=2.0,zorder=4,label="Peak water surface")
ax2.fill_between(sf,zc*FT,wc*FT,where=hc*FT>0.5,color="#e8112d",alpha=.45,zorder=3,
                 label="Overtopping")
ax2.set_xlim(0,sf[-1]); ax2.set_ylim(zc.min()*FT-8,max(wc.max(),zc.max())*FT+6)
ax2.set_xlabel("Distance along natural barrier, south → north (ft)")
ax2.set_ylabel("Elevation (ft NAVD88)")
ax2.grid(alpha=.3,zorder=1); ax2.legend(fontsize=9,loc="upper left")
pc=100*ov.sum()/ov.size
ax2.set_title(f"Barrier profile — overtopped along {pc:.0f}% of its {sf[-1]:,.0f} ft length; "
              f"max head over crest {np.nanmax(hc)*FT:.1f} ft",fontsize=11)
plt.savefig(f"{OUT}/c6_2d_preliminary.png",bbox_inches="tight");print("wrote map")
print(f"barrier: {sf[-1]:,.0f} ft long, overtopped {pc:.0f}%, max depth over crest {np.nanmax(hc)*FT:.2f} ft")
print(f"crest elev {zc.min()*FT:.1f}-{zc.max()*FT:.1f} ft, peak WSE {wc.min()*FT:.1f}-{wc.max()*FT:.1f} ft")
