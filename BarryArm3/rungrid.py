import os
import grid_search_functions as grid
import plot_functions as plot

if not os.path.exists('output'):
   os.makedirs('output')

if not os.path.exists('figures'):
   os.makedirs('figures')
   
#grid.makestalist("AK,AT,AV,CN","BAE,BAT,BMR,CAPN,CUT,DHY,DIV,EYAK,FID,FIRE,GHO,GLB,GLI,GOAT,HARP,HIN,KLU,KNK,L22K,M23K,P23K,PS11,PS12,PWL,Q23K,RAG,RC01,SAW,SCM,SKN,SLK,SWD,VMT,WAT6,WAT7,PMR,N25K,STLK,WACK,WAT1","BHZ,HHZ","*","2021-08-09 07:45:14") 
grid.creategrid(-149.9,-145.9,0.02,60.2,62.2,0.01,3.4,"output/stafile.csv")
#grid.prepwaveforms("output/stafile.csv","2021-08-09 07:45:14",180)  
#plot.plotRecSecNorm("output/before_ts_wfs.mseed","output/stafile.csv",61.242,-147.94)
#plot.plotWfUnnorm("output/before_ts_wfs.mseed","output/stafile.csv",61.242,-147.94)
#grid.locate("output/ttgridfile.npy","output/before_ts_wfs.mseed","output/lonlatgridfile.npy",1)
#grid.error(61.242,-147.94,0.90,"output/locfile.csv","output/lonlatgridfile.npy","output/strengthfile.npy")
#grid.size("output/stafile.csv","output/locfile.csv","output/before_ts_wfs.mseed")
####plot.plotGrid("output/lonlatgridfile.npy","output/strengthfile.npy")  #no need to use