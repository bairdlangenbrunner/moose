import numpy
import scipy.stats
import sys
from netCDF4 import Dataset
import netCDF4
import numpy.linalg
import datetime
import scipy.signal
import itertools
import numpy.random
import time

##########################################################################################
##########################################################################################
##########################################################################################

season_names = ['djf','mam','jja','son','ondjfm']

season='djf'
#season='mam'
#season='jja'
#season='son'
#season='annual'

# ua_corr_lat_lo, ua_corr_lat_hi, ua_corr_lon_lo, ua_corr_lon_hi = 30,40,170.,205.
# ts_corr_lat_lo, ts_corr_lat_hi, ts_corr_lon_lo, ts_corr_lon_hi = -25,-15,200.,230.

ncfile = Dataset('/ninod/baird/cmip5/obs/ERSSTv4/sst.mnmean.v4_invertlat_72x144regrid.nc', 'r', format='NETCDF4')
sst_data_orig = ncfile.variables['sst'][:]
sst_lat = ncfile.variables['lat'][:]
sst_lon = ncfile.variables['lon'][:]

global_nlat, global_nlon = sst_data_orig.shape[1:3]
global_lat_vals = sst_lat[:]
global_lon_vals = sst_lon[:]

##########################################################################################
##########################################################################################
##########################################################################################

# (YYYY,MM,DD)

#hist_start = netCDF4.datetime(1969,12,1)
#hist_end = netCDF4.datetime(2000,2,28)

hist_start = netCDF4.datetime(1979,12,1)
hist_end = netCDF4.datetime(2010,2,28)

sst_histmonths_list = []
sst_histclim_list = []
hist_stdevs_list = []

# OPEN HISTORICAL PERIOD PR DATA

file_name = '/ninod/baird/cmip5/obs/ERSSTv4/sst.mnmean.v4_invertlat_72x144regrid.nc'

ncfile = Dataset(file_name, 'r', format='NETCDF4')
sst_hist_data = ncfile.variables['sst'][:,:,:]
time_variable = ncfile.variables['time']
date_start = netCDF4.date2num(hist_start, time_variable.units, time_variable.calendar)
date_end = netCDF4.date2num(hist_end, time_variable.units, time_variable.calendar)
time_variable_converted = netCDF4.num2date(time_variable[:], time_variable.units, calendar='standard')
model_time = time_variable[:]

ts_lat_lo, ts_lat_hi, ts_lon_lo, ts_lon_hi = -25,-15,200.,230.; region = 'correlation_region'

# pull out only lat lon that are relevant
sst_lat_inds = numpy.where((global_lat_vals>=ts_lat_lo) & (global_lat_vals<=ts_lat_hi))[0]
sst_lon_inds = numpy.where((global_lon_vals>=ts_lon_lo) & (global_lon_vals<=ts_lon_hi))[0]
sst_hist_data = sst_hist_data[:, sst_lat_inds[0]:(sst_lat_inds[-1]+1), sst_lon_inds[0]:(sst_lon_inds[-1]+1)]

if season=='djf':		
	time_indices = numpy.array([(t.month in [12,1,2])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	sst_hist_data_seas = sst_hist_data[time_indices,:,:]

regional_nlat, regional_nlon = sst_hist_data_seas.shape[1:3]
sst_histclim_data = numpy.mean(sst_hist_data_seas, axis=0)
sst_histclim_stdev = numpy.std(sst_hist_data_seas, axis=0, ddof=1)
sst_hist_data_seas_means = sst_hist_data_seas.reshape((-1,3,regional_nlat,regional_nlon)).mean(axis=1)

print(sst_hist_data_seas_means.shape)

sst_hist_data_seas_spatial_means = numpy.mean(sst_hist_data_seas_means, axis=2)
sst_hist_data_seas_spatial_means = numpy.mean(sst_hist_data_seas_spatial_means, axis=1)

print(sst_hist_data_seas_spatial_means.shape)

N = sst_hist_data_seas_spatial_means.size
sst_stdev = numpy.std(sst_hist_data_seas_spatial_means)
std_err = sst_stdev/numpy.sqrt(N)
print(std_err)