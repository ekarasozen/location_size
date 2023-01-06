from obspy.clients.fdsn import Client
from obspy import UTCDateTime
client = Client("IRIS")
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



####################################################################################################################################################################################
#makestalist
#This is an informal script to aid in the creation of sta_file. Unlike the other programs, this is a script with hardwired values that users are likely to manipulate directly in order to assemble the specific set of stations they want. 
####################################################################################################################################################################################

def makestalist(network,station,channel,location,datetime):
    start_time = time.time()
    df =  pd.DataFrame(columns = ['network','station','channel','location','latitude','longitude','elevation'])
    starttime = UTCDateTime(datetime)  
    inv = client.get_stations(network=network,station=station,channel=channel,location=location,starttime=starttime,level='response')
    noi = len(inv) #number of networks
    nol = len(df) #number of locations in df. this is done to write different network rows properly into dataframe. 
    for i in range(noi):
        net = inv[i]
        netcod = inv[i].code
        print(netcod)
        nos = len(net)  #number of stations in each network
        for s in range(nos):
            stacod=net[s].code 
            stalat=net[s].latitude
            stalon=net[s].longitude
            staelv=net[s].elevation
            chn = net[s] 
            noc = len(chn) #number of channels in each station
            for c in range(noc):
               cha = chn[c].code
               loc = chn[c].location_code
            data = netcod,stacod,cha,loc,stalat,stalon,staelv
            df.loc[s+nol] = data 
        nol = len(df)
    df.append(df, ignore_index = True)        
    print(df)
    df.to_csv('output/stafile.csv',index=False)
    print("--- %s seconds ---" % (time.time() - start_time))
   

####################################################################################################################################################################################
#createttgrid
#This code defines a grid of source locations and pre-computes travel times from each source to all stations. This is a core computational function that users should not need to alter. 

def creategrid(lonmin,lonmax,lonnum,latmin,latmax,latnum,wavespeed,stafile):
    start_time = time.time()
    df = pd.read_csv(stafile,index_col=None,keep_default_na=False)
    nos = len(df) #number of stations
    g = Geod(ellps='WGS84')
    station = df['station']
    stalat = df['latitude']
    stalon = df['longitude']
    grdpts_x = np.array([i for i in np.arange(lonmin,lonmax,lonnum)],dtype=np.float32)
    grdpts_y = np.array([j for j in np.arange(latmin,latmax,latnum)],dtype=np.float32)
    yy, xx = np.meshgrid(grdpts_y,grdpts_x)
    lonlatgrid = np.array([yy, xx])
    distgrid = np.zeros((len(yy), len(xx[0]), nos)) #longitude rows are needed, latitude columns
    ttgrid = np.zeros((len(yy), len(xx[0]), nos))
    for s in range(nos):
        for j in range(len(yy)):
            for i in range(len(xx[0])):
                azimuth1, azimuth2, distance_2d = g.inv(stalon[s], stalat[s], xx[j,i], yy[j,i]) 
                distm = distance_2d
                distkm = float(distm)/1000
                distgrid[j,i,s] = distkm #distances are saved in meters
                tt = np.divide(float(distkm),float(wavespeed)) #calculate the traveltime
                ttgrid[j,i,s] = tt #save to traveltime table
    np.save("output/distgridfile", distgrid)  
    np.save("output/ttgridfile", ttgrid)  
    np.save("output/lonlatgridfile", lonlatgrid)
    print("--- %s seconds ---" % (time.time() - start_time))
    
####################################################################################################################################################################################
#prepwaveforms
#this code loads waveforms and prepares them for the stacking process. This is a core computational function that users should not need to alter. 
####################################################################################################################################################################################

def prepwaveforms(stafile,datetime,duration):
    start_time = time.time()
    s1 = UTCDateTime(datetime)  
    st = Stream() #initial stream
    df = pd.read_csv(stafile,index_col=None,keep_default_na=False)
    nos = len(df) #number of stations
    network = df['network']
    station = df['station']
    location = df['location']
    channel = df['channel']
    for s in range(nos):
        try:
            st += client.get_waveforms(network=network[s], station=station[s], location=location[s], channel=channel[s], starttime=s1, endtime=s1+duration, attach_response=True)
            #print(network[s],station[s],location[s],channel[s])
        except FDSNNoDataException:
            print('No data available for request. Creating a blank trace for this station: ' + station[s])
            tr = Trace()
            tr.stats.starttime=s1 
            tr.stats.network = network[s] 
            tr.stats.station = station[s]
            tr.stats.location = location[s]
            tr.stats.channel = channel[s]
            tr.stats.sampling_rate = 50 
            tr.stats.npts=duration*tr.stats.sampling_rate
            tr.data=np.zeros(tr.stats.npts)
            st += tr            
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=0.05, type='cosine')
    st.filter('bandpass', freqmin=0.01, freqmax=0.05)
    for tr in st:
        is_all_zero = np.all((tr.data == 0))
        if not is_all_zero:
      #  if tr.data != []: #creating empty traces didn't work because mseed doesn't save them. 
            tr.remove_response(output='DISP')
    print(st.__str__(extended=True))
    st.write('output/before_ts_wfs.mseed', format="MSEED")  
    print("--- %s seconds ---" % (time.time() - start_time))

####################################################################################################################################################################################
#locmethod
#This code includes different location methodologies that the locate function calls.
####################################################################################################################################################################################

def locmethod(st,method=1):
    print(st)
    if method ==1 or method == 2: 
        st.normalize()
    npts_all = [len(tr) for tr in st]
    npts = min(npts_all)
    if method ==2: 
        data = np.array([envelope(tr.data[:npts]) for tr in st])
    else:
        data = np.array([tr.data[:npts] for tr in st])
    nots = data.shape[1] #number of time samples
    noss = data.shape[0] # same with nos = len(st), using this for consistency 
    if method == 1: #amplitude stacking  
        name = "Amplitude stacking"
        print("Location method: ",name)
        stack = np.mean(data, axis=0)
        trs = Trace() #save stack as a trace
        trs.data=stack
        trs.stats.starttime=st[0].stats.starttime 
        trs.stats.station = "STCK"
        trs.stats.sampling_rate = 50 
        st += trs
        maxpower = np.max(np.abs(stack)) #max of the abs value of the stack
        print(maxpower)
    if method == 2: #amplitude envelope stacking, add envelope - this method does not work yet!!
        name = "Envelope stacking"
        print("Location method: ",name)
        stack = np.mean(data, axis=0)
        trs = Trace() #save stack as a trace
        trs.data=stack
        trs.stats.starttime=st[0].stats.starttime 
        trs.stats.station = "STCK"
        trs.stats.sampling_rate = 50 
        st += trs
        maxpower = np.max(np.abs(stack)) #max of the abs value of the stack
        print(maxpower)
    if method == 3: #semblance equation #1
        name = "Semblance eqn #1"
        print("Location method: ",name)
        snum = np.zeros((1, nots)) #semblance eqn numerator
        sden = np.zeros((1, nots)) # semblance eqn denumerator
        for ts in range(nots):
           sem1 = np.square(np.sum(data[:,ts]))
           snum[0,ts] = sem1
        sumn = np.sum(snum) #sum of numerator
        for ts in range(nots):
           sem2= np.sum(np.square(data[:,ts]))
           sden[0,ts] = sem2
        sumd = noss * np.sum(sden) #sum of denumerator
        maxpower = np.divide(sumn,sumd) #i.e. semblance
        print(maxpower)
    if method == 4: #semblance equation #2
        name = "Semblance eqn #2"
        print("Location method: ",name)
        sigm = np.zeros((noss, 1)) # sigma
        sden = np.zeros((1, nots)) # semblance eqn denumerator
        for s in range(noss):
           sig = np.sqrt((1/nots)*(np.sum(np.square(data[s,:]))))
           sigm[s,0] = sig
        #print(data[:,0]/sigm[:,0])
        for ts in range(nots):
           sem1= np.square(np.sum(data[:,ts]/sigm[:,0]))
           sden[0,ts] = sem1
        maxpower = 1/(nots * np.square(noss)) * np.sum(sden) #i.e. semblance
        print(maxpower)
    if method == 5: #semblance equation #3
        name = "Semblance eqn Ripepe"
        print("Location method: ",name)
        osum = np.zeros((noss-1)) #outer summation
        for i in range(noss-1): #should lead (noss!/2!*(noss-2)!) combinations, e.g. for 12 traces = 66 combinations
            gamma = np.zeros((noss)) #fix this range, it works but could be better
            for j in range(i+1, noss):
               cov1 = np.cov(data[i,:],data[j,:])
               cov = cov1[1,0]
               std = np.std(data[i,:])*np.std(data[j,:])
               gamma[j] = np.divide(cov,std)
            isum = np.sum(gamma) #inner summation
            osum[i] = isum # outer summation
        maxpower = np.divide(np.sum(osum),np.sum(range(1,noss))) #range (1,noss) means 1,2,3,4,5,6,7,8,9,10,11 i.e. noss-1
        print(maxpower)
    return maxpower


####################################################################################################################################################################################
#shiftstack
#This code iterates through the i x j grid of source locations. For each grid point the code shifts and stacks waveforms using the pre-computed travel times. For each grid point, some type of stack amplitude measure is retained. Eventually this function should also include an estimate of the location error (perhaps based on the second derivative at the stack maximum?) This is a core computational function that users should not be tweaking. 
####################################################################################################################################################################################


def locate(ttgridfile,before_ts_wfs,lonlatgridfile,method=1):
    start_time = time.time()
    df =  pd.DataFrame(columns = ['latitude','longitude','spower','rtime'])
    ttgrid = np.load(ttgridfile)
    lonlen = ttgrid.shape[0]
    latlen = ttgrid.shape[1]
    lonlatgrid = np.load(lonlatgridfile)
    strength = np.zeros((lonlen, latlen)) 
    maxpowerall = 0
    for i in range(lonlen):
        for j in range(latlen):
            st = Stream() #initial stream
            st = read(before_ts_wfs)
            nos = len(st) #number of stations in the stream
            for s in range(nos):
                maxtt = np.max(ttgrid[i,j,:])
                st[s].trim((st[s].stats.starttime+np.around(ttgrid[i,j,s],2)),(np.around(maxtt-ttgrid[i,j,s],2))) #does maxtt needed???
                st[s].stats.starttime = st[s].stats.starttime-ttgrid[i,j,s] # this is needed for record section. not for the stacking
            maxpower = locmethod(st,method)                 
            strength[i,j] = maxpower #save to traveltime table
            if maxpower > maxpowerall:
               maxpowerall = maxpower
            else: 
               continue
    print(maxpowerall)
    maxinx=np.where(strength == maxpowerall)
    cordinx = list(zip(maxinx[0], maxinx[1]))
    for cord in cordinx:
        evlat = lonlatgrid[0][cord]
        evlon = lonlatgrid[1][cord]
    print(cordinx)
    rtime = (time.time() - start_time)
    location=evlat,evlon,maxpowerall,rtime
    df.loc[0] = location 
    print(df)
    print("--- %s seconds ---" % (time.time() - start_time))
    df.to_csv('output/locfile.csv',index=False)
    np.save("output/strengthfile", strength)


####################################################################################################################################################################################
#error
#This code calculates  absolute and relative location errors
####################################################################################################################################################################################

def error(gtevlat,gtevlon,conf,locfile,lonlatgridfile,stacked_strengthfile):

    g = Geod(ellps='WGS84')
    #conf=0.90
    #gtevlat=58.6948   #Taku River
    #gtevlon=-133.3919 #Taku River
    	
    evdf = pd.read_csv(locfile,index_col=None,keep_default_na=False)
    errdf =  pd.DataFrame(columns = ['spower','threshold','loc_error','confidence','err1','err2','rblx','rbly','rurx','rury'])
    #rblx: rectangle bottom left X
    #rblx: rectangle bottom left Y
    #rblx: rectangle upper right X
    #rblx: rectangle upper right X

    #plot the entire map
    
    lonlatgrid = np.load(lonlatgridfile)
    stacked_strength = np.load(stacked_strengthfile)
    yy= lonlatgrid[0]
    xx= lonlatgrid[1]
    
    
    fig, (ax) = plt.subplots(1,1)
    im = ax.pcolormesh(xx,yy, stacked_strength,shading='nearest',cmap=cm.hot_r, vmin=0, vmax=1)
    ax110 = fig.add_axes([0.92, 0.11, 0.01, 0.77])
    fig.colorbar(im, cax=ax110)    
    fig.savefig('output/stacked_strength_map_' + str(conf) + '.png', bbox_inches='tight')
    power = np.amax(stacked_strength)
    threshold = power*conf
    stacked_strength[stacked_strength < threshold] = np.NaN 
    
    notnan = np.argwhere(~np.isnan(stacked_strength)) #get the indexes where there are numbers
    
    sub_lonlatgrid = [] #put those indexes into an array
    cordinx = list(zip(notnan[:,0], notnan[:,1]))
    for cord in cordinx:
        lat = lonlatgrid[0][cord]
        lon = lonlatgrid[1][cord]
        sub_lonlatgrid.append([lon,lat,cord[0],cord[1]])
    sub_lonlatgrid = np.array(sub_lonlatgrid)
    
    maxx_inx=np.argwhere(sub_lonlatgrid == np.max(sub_lonlatgrid[:,0])) #maxx
    x1 = sub_lonlatgrid[maxx_inx[0][0]][0:2]
    maxx = int(sub_lonlatgrid[maxx_inx[0][0]][2:3])  #for the figure below
    
    minx_inx=np.argwhere(sub_lonlatgrid == np.min(sub_lonlatgrid[:,0])) #maxx
    x2 =  sub_lonlatgrid[minx_inx[0][0]][0:2]
    minx =  int(sub_lonlatgrid[minx_inx[0][0]][2:3])  #for the figure below
    
    miny_inx=np.argwhere(sub_lonlatgrid == np.min(sub_lonlatgrid[:,1])) #maxy
    y1 = sub_lonlatgrid[miny_inx[0][0]][0:2]
    miny = int(sub_lonlatgrid[miny_inx[0][0]][3:4]) #for the figure below
    
    maxy_inx=np.argwhere(sub_lonlatgrid == np.max(sub_lonlatgrid[:,1])) #maxy
    y2 = sub_lonlatgrid[maxy_inx[0][0]][0:2]
    maxy = int(sub_lonlatgrid[maxy_inx[0][0]][3:4])  #for the figure below
    
    rect_x = x1[0],x2[0]
    rect_y = y1[1],y2[1]
    rect_blx = np.min(rect_x) #rectangle bottom left x
    rect_bly = np.min(rect_y) #rectangle bottom left y
    rect_urx = np.max(rect_x) #rectangle upper right x
    rect_ury = np.max(rect_y) #rectangle upper right y
    
    azimuth1, azimuth1, distance_2d1 = g.inv(x1[0], x1[1], x2[0], x2[1])
    err1 = float(distance_2d1)/1000
    print(err1)
    
    azimuth2, azimuth2, distance_2d2 = g.inv(y1[0], y1[1], y2[0], y2[1])
    err2 = float(distance_2d2)/1000
    print(err2)
    
    #absolute location error, nothing to do with previous calclation
    azimuth3, azimuth3, distance_2d3 = g.inv(evdf.longitude[0],evdf.latitude[0],gtevlon,gtevlat)
    gterr = float(distance_2d3)/1000
    loc_error = power,threshold,gterr,conf,err1,err2,rect_blx,rect_bly,rect_urx,rect_ury
    errdf.loc[0] = loc_error 
    print(errdf)            
    errdf.to_csv('output/locerror_' + str(conf) + '.csv',index=False)
    
    #plot the isolated peak
   
    #sub_strength = stacked_strength[minx:maxx,miny:maxy]
    #sub_yy= yy[minx:maxx,miny:maxy]
    #sub_xx= xx[minx:maxx,miny:maxy]
    #
    #fig1, (ax1) = plt.subplots(1,1)
    #im1 = ax1.pcolormesh(xx,yy, stacked_strength,shading='nearest',cmap=cm.hot_r,vmin=0, vmax=1)
    #ax111 = fig1.add_axes([0.92, 0.11, 0.01, 0.77])
    #fig1.colorbar(im1, cax=ax111)    
    #fig1.savefig('stacked_strength_nan_map_' +  str(conf) + '.png', bbox_inches='tight')
    #
    #fig2, (ax2) = plt.subplots(1,1)
    #im2 = ax2.pcolormesh(sub_xx,sub_yy, sub_strength,shading='nearest',cmap=cm.hot_r,vmin=0, vmax=1)
    #ax112 = fig2.add_axes([0.92, 0.11, 0.01, 0.77])
    #fig2.colorbar(im2, cax=ax112)    
    ##np.savetxt("sub_lonlatgrid", sub_lonlatgrid)
    #
    #fig2.savefig('sub_strength_map_' + str(conf) + '.png', bbox_inches='tight')
    
####################################################################################################################################################################################
#size
#This code calculates the volume of the event
####################################################################################################################################################################################

def size(stafile,locfile,before_ts_wfs):
    g = Geod(ellps='WGS84')
    a = 1.4719977544177063 
    b = 36.31567785721356
    sizedf =  pd.DataFrame(columns = ['volume'])
    ampdf =  pd.DataFrame(columns = ['absmaxamp','maxamp','minamp','envabsmaxamp','dist','station[s]'])
    evdf = pd.read_csv(locfile,index_col=None,keep_default_na=False)
    stdf = pd.read_csv(stafile,index_col=None,keep_default_na=False)
    station = stdf['station']
    stalat = stdf['latitude']
    stalon = stdf['longitude']
    st = read(before_ts_wfs)
    npts_all = [len(tr) for tr in st]
    npts = min(npts_all)
    data = np.array([tr.data[:npts] for tr in st])
    dataenv = np.array([envelope(tr.data[:npts]) for tr in st])
    noss = data.shape[0] # same with nos = len(st), using this for consistency 
    for s in range(noss):
        absmaxamp = np.max(np.abs(data[s])) #max of the abs value of the stack
        envabsmaxamp = np.max(np.abs(dataenv[s])) #max of the abs value of the stack
        maxamp = np.max((data[s])) #max of the abs value of the stack
        minamp = np.min((data[s])) #max of the abs value of the stack
        azimuth, azimuth, distance_2d = g.inv(evdf.longitude[0],evdf.latitude[0],stalon[s],stalat[s])
        dist = float(distance_2d)/1000
        amplitude = absmaxamp,maxamp,minamp,envabsmaxamp,dist,station[s]
        ampdf.loc[s] = amplitude
    ampdf['distance (deg)'] = ampdf['dist']/111    
    abs_amp=ampdf.loc[(ampdf['distance (deg)'] >= 1), 'absmaxamp']        
    #log10 (y) = b + a*log10(x)
    mean = np.mean(abs_amp)
    volume = np.exp((np.log(mean)*a)+b)    
    sizedf.loc[0] = volume 
    print(sizedf)
    sizedf.to_csv('output/sizefile.csv',index=False)

        