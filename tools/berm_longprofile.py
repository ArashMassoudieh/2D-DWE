#!/usr/bin/env python3
"""Longitudinal profile of the C6 berm: existing ground, carved crest, flood WSE.
Writes outputs_2D/berm_profile.csv and Report_for_C6/Figures/NM15-097_berm_longprofile.png"""
import numpy as np, csv, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from osgeo import gdal, ogr
gdal.UseExceptions(); ogr.UseExceptions()
FT=3.280839895
SITE='/home/arash/Dropbox/Chloeta/NHA/WS3_GIS/NM15-097'; DWE='/home/arash/Projects/2D-DWE/'
d0=gdal.Open(SITE+'/outputs_2D/cliped_dem.tif'); gt=d0.GetGeoTransform()
Z0=d0.GetRasterBand(1).ReadAsArray().astype('f8')
d1=gdal.Open(SITE+'/outputs_2D/cliped_dem_berm.tif'); Z1=d1.GetRasterBand(1).ReadAsArray().astype('f8')
NY,NX=Z0.shape
def load(p):
    ds=gdal.Open(p); b=ds.GetRasterBand(1); a=b.ReadAsArray().astype('f8'); nd=b.GetNoDataValue()
    return np.where(a==nd,np.nan,a) if nd is not None else a
W6=load(DWE+'out_c6_1m_rain/wse_max.tif'); W24=load(DWE+'out_c6_1m_24hr/wse_max_snap.tif')
D6=load(DWE+'out_c6_1m_rain/depth_max.tif'); D24=load(DWE+'out_c6_1m_24hr/depth_max_snap.tif')
src=ogr.Open(SITE+'/outputs_2D/berm_centerline.shp'); lyr=src.GetLayer(0)
f=lyr.GetNextFeature(); g=f.GetGeometryRef()
P=np.array([g.GetPoint(i)[:2] for i in range(g.GetPointCount())])
seg=np.hypot(*np.diff(P,axis=0).T); s=np.r_[0,np.cumsum(seg)]
tt=np.arange(0,s[-1],1.0); cx=np.interp(tt,s,P[:,0]); cy=np.interp(tt,s,P[:,1]); sft=tt*FT
def samp(A,xx,yy):
    i=((xx-gt[0])/gt[1]).astype(int); j=((yy-gt[3])/gt[5]).astype(int)
    ok=(i>=0)&(i<NX)&(j>=0)&(j<NY)
    o=np.full(len(xx),np.nan); o[ok]=A[j[ok],i[ok]]; return o
dx=np.gradient(cx); dy=np.gradient(cy); L=np.hypot(dx,dy); nxv,nyv=-dy/L,dx/L
off=np.arange(-2,2.1,0.5)
ground=np.nanmax(np.vstack([samp(Z0,cx+o*nxv,cy+o*nyv) for o in off]),axis=0)*FT
crest =np.nanmax(np.vstack([samp(Z1,cx+o*nxv,cy+o*nyv) for o in off]),axis=0)*FT
def fw(dep,wse):
    best=np.full(len(cx),np.nan)
    for o in range(2,41,2):
        d=samp(dep,cx+o*nxv,cy+o*nyv); w=samp(wse,cx+o*nxv,cy+o*nyv)
        best=np.fmax(best,np.where((d*FT)>=0.5,w,np.nan))
    return best*FT
wse=np.fmax(fw(D6,W6),fw(D24,W24))
with open(SITE+'/outputs_2D/berm_profile.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['station_ft','easting','northing','existing_ground_ft',
                                  'berm_crest_ft','height_ft','flood_wse_ft','freeboard_ft'])
    for k in range(len(sft)):
        wv='' if not np.isfinite(wse[k]) else round(wse[k],2)
        fb='' if not np.isfinite(wse[k]) else round(crest[k]-wse[k],2)
        w.writerow([round(sft[k],1),round(cx[k],2),round(cy[k],2),round(ground[k],2),
                    round(crest[k],2),round(crest[k]-ground[k],2),wv,fb])
fig,(ax,ax2)=plt.subplots(2,1,figsize=(14,8.4),height_ratios=[2.1,1.0],sharex=True)
ax.fill_between(sft,ground.min()-3,ground,color='#c8b89a',zorder=1)
ax.fill_between(sft,ground,crest,where=crest>ground+0.02,color='#e07000',alpha=0.5,zorder=2,label='Berm fill (new)')
ax.plot(sft,ground,color='#5a4a32',lw=1.9,zorder=4,label='Existing ground / existing berm')
ax.plot(sft,crest,color='#8a3800',lw=2.3,zorder=5,label='Design crest (carved)')
ax.plot(sft,wse,color='#0050b3',lw=1.7,zorder=6,label='Flood water surface (envelope, both storms)')
for xv in (315,587): ax.axvline(xv,color='k',lw=1.0,ls=':',zorder=7)
ax.annotate('flat 6153.4 ft',(157,6156.6),ha='center',fontsize=9,style='italic')
ax.annotate('constant grade $-$0.200%',(1850,6156.6),ha='center',fontsize=9.5,style='italic')
ax.set_ylim(ground.min()-3,6158.0); ax.set_ylabel('Elevation (ft, NAVD88)',fontsize=10)
ax.grid(alpha=0.25,lw=0.5); ax.legend(loc='upper right',fontsize=9,framealpha=0.94)
ax.set_title('C6 berm — longitudinal profile, full alignment carved (3,117 ft, 6 ft crest, 3:1)',
             fontsize=12.5,fontweight='bold')
fb=crest-wse
ax2.axhline(2.0,color='#c1121f',lw=1.3,ls='--',label='2 ft design freeboard')
ax2.plot(sft,fb,color='#0050b3',lw=1.9,label='Crest above the flood water surface')
ax2.fill_between(sft,2.0,fb,where=fb<2.0,color='#c1121f',alpha=0.3)
for xv in (315,587): ax2.axvline(xv,color='k',lw=1.0,ls=':')
ax2.set_xlim(0,sft[-1]); ax2.set_ylim(0,max(6,np.nanmax(fb)*1.1))
ax2.set_xlabel('Station along the centerline (ft)  —  south-west end at 0',fontsize=10)
ax2.set_ylabel('Freeboard (ft)',fontsize=10); ax2.grid(alpha=0.25,lw=0.5); ax2.legend(fontsize=9,ncol=2)
plt.tight_layout()
plt.savefig(SITE+'/Report_for_C6/Figures/NM15-097_berm_longprofile.png',dpi=200,bbox_inches='tight')
ok=np.isfinite(fb)
print('crest height above ground: max %.2f ft'%(crest-ground).max())
print('freeboard where loaded: min %.2f  mean %.2f ft; below 2 ft at %d of %d stations'%(
    np.nanmin(fb[ok]),np.nanmean(fb[ok]),int((fb[ok]<1.995).sum()),int(ok.sum())))
print('wrote berm_profile.csv and NM15-097_berm_longprofile.png')
