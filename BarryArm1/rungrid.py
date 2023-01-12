import os
import sys
sys.path.append("..")
import grid_search_functions as grid
import plot_functions as plot

if not os.path.exists('output'):
   os.makedirs('output')
if not os.path.exists('figures'):
   os.makedirs('figures')

grid.makestalist("AK,AV,AT,CN","BAE,BAW,BMR,BRLK,BRSE,CAPN,CUT,DHY,DIV,EYAK,FID,GHO,GLI,GOAT,HARP,HIN,KLU,KNK,M23K,P23K,PS11,PS12,PWL,Q23K,RAG,RC01,SAW,SCM,SKN,SLK,SSN,SWD,VMT,WAT6,WAT7,PMR,N25K,SPCN,SPCP,STLK,WAT1","BHZ,HHZ","*","2020-10-05 05:02:14") 
grid.creategrid(-150.2,-146.2,0.02,60.2,62.2,0.01,3.4,"output/stafile.csv")
grid.prepwaveforms("output/stafile.csv","2020-10-05 05:02:14",180)
plot.plotRecSecNorm("output/before_ts_wfs.mseed","output/stafile.csv",61.153,-148.163)
plot.plotWfUnnorm("output/before_ts_wfs.mseed","output/stafile.csv",61.153,-148.163)
grid.locate("output/ttgridfile.npy","output/before_ts_wfs.mseed","output/lonlatgridfile.npy",1)
grid.error(61.153,-148.163,0.90,"output/locfile.csv","output/lonlatgridfile.npy","output/strengthfile.npy")
grid.size("output/stafile.csv","output/locfile.csv","output/before_ts_wfs.mseed")
####plot.plotGrid("output/lonlatgridfile.npy","output/strengthfile.npy")  #no need to use