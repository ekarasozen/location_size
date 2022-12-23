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

a = 1.4719977544177063 
b = 36.31567785721356

ampdf =  pd.read_csv('amplitudefile.csv', index_col=False,)
ampdf['distance (deg)'] = ampdf['dist']/111

print(ampdf.info())
print(ampdf)
abs_amp=ampdf.loc[(ampdf['distance (deg)'] >= 1), 'absmaxamp']


#log10 (y) = b + a*log10(x)
mean = np.mean(abs_amp)
print (mean)

volume = np.exp((np.log(mean)*a)+b)

#volume = np.exp((np.log(mean) - b)/a)
print(volume)