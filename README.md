README file for landslide location and size estimation codes. 

These codes are described in Karasozen & West (submitted in 2023), Toward the rapid assessment of potential tsunamigenic landslides. 

This folder contains 3 landslide examples and each example includes: 

1. a control script, rungdrid.py
2. a functions script, grid_seach_functions.py, that stores core scripts that prepare waveforms, create the grid and solve for location and size.  
3. a plotting script, plot_functions.py, that includes necessary plotting codes. 

rungrid.py can be executed by uncommenting each line and running rungrid.py at every step. Theses scripts will create two folders:

1. an output folder that contains results in obspy, numpy and pandas format, 
2. a figures folder.  

The code is not optimized for speed. 

Latest version: 6 January 2023. 
