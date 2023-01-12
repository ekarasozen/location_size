from obspy.clients.fdsn import Client
#client = Client(base_url="https://earthquake.alaska.edu", timeout=600)
client = Client("IRIS")
from obspy.clients.iris import Client
client_distaz = Client()
from obspy import UTCDateTime
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy import read
import numpy as np
from numpy import loadtxt
import pandas as pd
from obspy import Stream
from obspy import Trace
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import datetime
import matplotlib.dates as mdates
from pyproj import Geod
from matplotlib.transforms import blended_transform_factory

def plotGrid(lonlatgridfile,stacked_strengthfile):
    fig, (ax1) = plt.subplots(1,1)
    lonlatgrid = np.load(lonlatgridfile)
    stacked_strength = np.load(stacked_strengthfile)
    yy= lonlatgrid[0]
    xx= lonlatgrid[1]
    #plt.plot(xx, yy, marker='.', color='k', linestyle='none')
    #plt.show()
    maxamp = abs(stacked_strength).max()/2         # these two lines just adjust the color scale
    print(abs(stacked_strength).max())
    minamp = abs(stacked_strength).min()
 #  #t, f = np.meshgrid(st[0].times("matplotlib"), freqs)
#    im1 = ax1.pcolormesh(xx,yy, np.transpose(stacked_strength),shading='nearest')
    im1 = ax1.pcolormesh(xx,yy, stacked_strength,shading='nearest',cmap=cm.hot_r,vmin=0,vmax=1)
    ax111 = fig.add_axes([0.92, 0.11, 0.01, 0.77])
    fig.colorbar(im1, cax=ax111)    
    #ax1.set_ylabel('freq. (Hz)')
    #ax1.set_xlim(timevector[0], timevector[-1])
    #ax1.set_ylim(freqs[-1], freqs[0])
    #ax1.xaxis.set_major_formatter(date_format)
    #ax1.xaxis_date()
    #ax1.axes.xaxis.set_ticklabels([])
    #ax2.set_xlabel("Time [UTC]" )
    #ax1.text(0.01, 0.88, text, transform=ax2.transAxes, fontsize=10, fontweight='bold', color='white')
    fig.savefig('figures/stacked_strength', bbox_inches='tight')
    
    
def plotRecSecNorm(st,stafile,gtevlat,gtevlon):
    st = read (st)                 
  #  print(st.__str__(extended=True))
    df = pd.read_csv(stafile,index_col=None,keep_default_na=False)
    nos = len(df) #number of stations
    g = Geod(ellps='WGS84')
    station = df['station']
    stalat = df['latitude']
    stalon = df['longitude']
  #  print(station)
    for s in range(nos):
        azimuth1, azimuth2, distance_2d = g.inv(stalon[s], stalat[s], gtevlon, gtevlat) 
        distm = distance_2d
        st[s].stats.distance = distm 
       # print(st[s].stats.distance)
    st.sort(keys=['distance'])
    print(st.__str__(extended=True))
    fig = plt.figure()
    st.plot(type="section",plot_dx=20e3, 
#    st.plot(type="section",plot_dx=20e3,norm_method='stream',  #if you want un-normalized record section
      time_down=True, linewidth=.75, grid_linewidth=.25,
      show=False, fig=fig)
    ax = fig.axes[0]
    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for tr in st:
        ax.text(tr.stats.distance / 1e3, 1.0, tr.stats.station, rotation=270,
            va="bottom", ha="center", transform=transform, zorder=10)
    fig.savefig('figures/record_section_norm', bbox_inches='tight')


def plotWfUnnorm(st,stafile,gtevlat,gtevlon): #function to plot all waveforms at one subplot; alternative way to see un-normalized record sections
    st = read (st)                 
    print(st.__str__(extended=True))
    df = pd.read_csv(stafile,index_col=None,keep_default_na=False)
    nos = len(df) #number of stations
    g = Geod(ellps='WGS84')
    station = df['station']
    stalat = df['latitude']
    stalon = df['longitude']
   # print(station)
    for s in range(nos):
        azimuth1, azimuth2, distance_2d = g.inv(stalon[s], stalat[s], gtevlon, gtevlat) 
        distm = distance_2d
       # print(distm,station[s])
        st[s].stats.distance = distm 
    st.sort(keys=['distance'])
    print(st.__str__(extended=True))
    fig = plt.figure()
    tv = st[0].times("matplotlib")
    fig = plt.figure(figsize=(12, 24))
    for s in range(nos):
        tr = st[s]
        txt = (tr.stats.station + "." + tr.stats.channel)
        ax = fig.add_subplot(1, 1,1)
        ax.plot(tr.times("matplotlib"), tr.data-(s*0.0000002), "k-", linewidth=0.8)
        ax.set_xlim(tv[0], tv[-1]) 
        ax.text(0.88, 0.905-(s*(1/(nos+5))), txt, transform=ax.transAxes, fontsize=8, fontweight='bold', verticalalignment='top')
  #      ax.axvline(starttime, lw=0.8, c='darkblue', ls='--', label='event onset')
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        ax.xaxis_date()
    fig.autofmt_xdate()
    fig.savefig('figures/record_section_unnorm', bbox_inches='tight')

