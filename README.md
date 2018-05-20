# moose:  multi-objective optimization of subensembles

This repository contains scripts for recreating the analysis and figures in [Langenbrunner and Neelin (2017, GRL)](https://doi.org/10.1002/2017GL075226).

## Jupyter Notebook and appropriate packages
To run the scripts as-is, you'll need a working installation of Jupyter Notebook and a couple required Python packages for reading NetCDF data and plotting maps:
* netCDF4
* basemap

## ```SCRIPTS``` and ```DATA``` directories

To reproduce Figs. 1–4 in Langenbrunner and Neelin (2017), access the appropriate notebooks via the ```SCRIPTS``` directory.  This is downloaded if you clone the ```moose``` repository.

You'll also need to download the ```DATA``` folder, which is stored as a zipped directory on Google Drive [here](https://drive.google.com/file/d/1MNcl8nYhaHLvkQ1kfOdQCnG8tSSVcks4/view?usp=sharing).  Unzip this directory into the moose repository to give scripts access to all the necessary data.  __Note:  ```DATA.zip``` is ~550 MB, and unzipped it expands to ~1.5 GB.__

The ```DATA``` folder contains:
* December-January-February average precipitation, 200 hPa wind, and skin temperature for 36 CMIP5 models in the ```DATA/cmip5_data``` directory, for:
  * 1980-2010 (historical period climatologies)
  * 2070-2100 (end-of-century change calculations)
* Precipitation, 200 hPa wind, and SST observations during DJF 1980-2010 in the ```DATA/obs_data``` directory
* Pre-calculated Pareto front information to reproduce Figs. 1–4 in the ```DATA/subensemble_data```
* Some miscellaneous ```.npy``` and NetCDF4 files for applying land masks and other final steps
