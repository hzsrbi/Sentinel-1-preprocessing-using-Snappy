### This repository contains Python script for Sentinel-1 image pre-processing using snappy. 

1. Before running the code, you need to install the [Sentinel Toolbox Application (SNAP)](https://step.esa.int/main/download/snap-download/) and [configure Python to use the SNAP-Python (snappy) interface](https://senbox.atlassian.net/wiki/spaces/SNAP/pages/50855941/Configure+Python+to+use+the+SNAP-Python+snappy+interface)

2. The code reads in unzipped Sentinel-1 GRD products (EW and IW modes).
    - Sentinel-1 images can be downloaded from:
      - [Sentinels Scientific Data Hub](https://scihub.copernicus.eu/dhus/#/home)  or
      - [Alaska Satellite Facility](https://vertex.daac.asf.alaska.edu/#)
3. Sentinel-1 pre-processing steps:
    
    (1) Apply orbit file  
    (2) Thermal noise removal  
    (3) Radiometric calibration  
    (4) Speckle filtering  
    (5) Terrain correction  
    (6) Subset

   
   This is the general pre-processing steps for Sentinel-1. Since each step is written in a separate function, it can be cutomized based on user needs (Tips: If you would like to modify the parameters, you can refer to the graph builder file (.xml) which can be created in the Graph Builder of SNAP software.).
   IW images are downsampled from 10m to 40m (the same resolution as EW images) in the terrain correction step.




