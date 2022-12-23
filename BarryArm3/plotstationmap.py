import pygmt
import pandas as pd
import numpy as np
name = "Seward"

gtevlat=64.745 #FourMiles 2022
gtevlon=-149.126 #FourMiles 2022


fig = pygmt.Figure()
pygmt.config(MAP_FRAME_TYPE="plain",FORMAT_GEO_MAP="ddd",
#MAP_GRID_CROSS_SIZE_PRIMARY=0.3,
#MAP_GRID_CROSS_SIZE_SECONDARY=0.3,
MAP_TICK_LENGTH_PRIMARY=0.15,
MAP_TICK_LENGTH_SECONDARY=0.10,
FONT_TITLE=8,
FONT_HEADING=8,
FONT_LABEL=8,
FONT_TAG=8,
FONT_ANNOT_PRIMARY=8,
FONT_ANNOT_SECONDARY=6,
COLOR_NAN="white")


stdf = pd.read_csv("stafile.csv",index_col=None,keep_default_na=False)
evdf = pd.read_csv("locfile.csv",index_col=None,keep_default_na=False)
ledf = pd.read_csv("locerror_0.9.csv",index_col=None,keep_default_na=False)

#FOR A ZOOM IN MAP
#latmin = (np.min(gtevlat))-0.5
#latmax = (np.max(gtevlat))+0.5
#lonmin = (np.min(gtevlon))-1.0
#lonmax = (np.max(gtevlon))+1.0

#FOR A ZOOM OUT MAP
latmin = (np.min(gtevlat))-2.0
latmax = (np.max(gtevlat))+2.0
lonmin = (np.min(gtevlon))-4.0
lonmax = (np.max(gtevlon))+4.0

region=[lonmin, lonmax, latmin, latmax]
val=str((lonmin+lonmax)/2)
projection='S'+ val + '/90/2i'

lonlatgrid = np.load("lonlatgridfile.npy")
stacked_strength = np.load("strengthfile.npy")
y= lonlatgrid[0]
x= lonlatgrid[1]
z = stacked_strength
xx, yy, zz = x.flatten(), y.flatten(), z.flatten()
region2=[np.amin(x), np.amax(x), np.amin(y), np.amax(y)]

ds = pygmt.xyz2grd(x=xx, y=yy, z=zz,region=region2, spacing=(0.02,0.01))
grid = pygmt.datasets.load_earth_relief(resolution="15s", region=region,)

fig.basemap(region=region, projection=projection, frame=["WSne","xaf", "yaf"])
fig.grdimage(grid="/Users/ezgikarasozen/Documents/GitHub/barry/grid_search/results/gmt_files/aec_zoomed.grd",projection=projection,cmap='/Users/ezgikarasozen/Documents/GitHub/barry/grid_search/results/gmt_files/gray_blue_AEC_simple.cpt',shading="+a45+nt0.5")
pygmt.makecpt(cmap='hot',series='0/1/0.1',continuous=True,reverse=True,)
#fig.coast(water="white", borders="1/0.5p", shorelines="1/0.5p", frame=['+t"0.5 deg."'])
fig.coast(water="white", borders="1/0.5p", shorelines="1/0.5p")
fig.grdimage(grid=ds,projection=projection,cmap=True,transparency="30",)
fig.colorbar(position="JMR+o0.5c/0c+w4c",box=True,frame=["x+lMax. stack"],)
fig.plot(x=stdf.longitude, y=stdf.latitude, style="i0.2c", color="blue", pen="0.5p,black")
fig.plot(x=evdf.longitude, y=evdf.latitude, style="a0.25c", color="darkred", pen="0.5p,black")
fig.plot(x=gtevlon, y=gtevlat, style="a0.25c", color="darkgreen", pen="0.5p,black")
#fig.plot(x=gtevlon1, y=gtevlat1, style="a0.25c", color="darkgreen", pen="0.5p,black")
#fig.plot(x=gtevlon2, y=gtevlat2, style="a0.25c", color="darkblue", pen="0.5p,black")
#rectangle = [[ledf.rblx[0],ledf.rbly[0],ledf.rurx[0],ledf.rury[0]]]
#fig.plot(data=rectangle, style="r+s", pen="0.5p,black")
#fig.legend(position="jBR+jBR+o0.4c", box="+gwhite+p1p", S=0.3)


fig.show()
fig.savefig(name + "_grid_wstations.pdf")
