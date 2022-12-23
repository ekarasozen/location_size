import os
import grid_search_functions as grid
import plot_functions as plot

#grid.makestalist("AK,AT,AV,CN","BAE,BAT,BMR,CAPN,CUT,DHY,DIV,EYAK,FID,FIRE,GHO,GLB,GLI,GOAT,HARP,HIN,KLU,KNK,L22K,M23K,P23K,PS11,PS12,PWL,Q23K,RAG,RC01,SAW,SCM,SKN,SLK,SWD,VMT,WAT6,WAT7,PMR,N25K,STLK,WACK,WAT1","BHZ,HHZ","*","2021-08-09 07:45:14") 
#grid.creategrid(-149.9,-145.9,0.02,60.2,62.2,0.01,3.4,"stafile.csv")
#grid.prepwaveforms("stafile.csv","2021-08-09 07:45:14",180)  
#plot.plotRecSecNorm("before_ts_wfs.mseed","stafile.csv",61.242,-147.94)
#plot.plotWfUnnorm("before_ts_wfs.mseed","stafile.csv",61.242,-147.94)
#grid.locate("ttgridfile.npy","before_ts_wfs.mseed","lonlatgridfile.npy",1)
#grid.error(61.242,-147.94,0.90,"locfile.csv","lonlatgridfile.npy","strengthfile.npy")
####plot.plotGrid("lonlatgridfile.npy","strengthfile.npy")  #no need to use