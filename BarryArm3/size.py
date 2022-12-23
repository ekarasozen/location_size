from obspy.clients.fdsn import Client
from obspy import UTCDateTime
#client = Client(base_url="https://earthquake.alaska.edu", timeout=600)
client = Client("IRIS")
from obspy.clients.iris import Client
client_distaz = Client()
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy import Stream
from obspy import Trace
from obspy import read, read_inventory
from obspy.signal.filter import envelope
import numpy as np
from numpy import loadtxt
import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Geod
import time
import matplotlib.cm as cm

g = Geod(ellps='WGS84')
evdf = pd.read_csv("locfile.csv",index_col=None,keep_default_na=False)
stdf = pd.read_csv("stafile.csv",index_col=None,keep_default_na=False)
station = stdf['station']
stalat = stdf['latitude']
stalon = stdf['longitude']
#print(station)

st = read("before_ts_wfs.mseed")  
#print(st.__str__(extended=True))
#st.normalize()
npts_all = [len(tr) for tr in st]
npts = min(npts_all)
data = np.array([tr.data[:npts] for tr in st])
dataenv = np.array([envelope(tr.data[:npts]) for tr in st])
noss = data.shape[0] # same with nos = len(st), using this for consistency 
for s in range(noss):
    #print(station[s],st[s])
    absmaxamp = np.max(np.abs(data[s])) #max of the abs value of the stack
    envabsmaxamp = np.max(np.abs(dataenv[s])) #max of the abs value of the stack
    maxamp = np.max((data[s])) #max of the abs value of the stack
    minamp = np.min((data[s])) #max of the abs value of the stack
    #print(maxamp,minamp,absmaxamp,envabsmaxamp)
    azimuth, azimuth, distance_2d = g.inv(evdf.longitude[0],evdf.latitude[0],stalon[s],stalat[s])
    dist = float(distance_2d)/1000
    print("BA3",absmaxamp,maxamp,minamp,envabsmaxamp,dist,station[s])