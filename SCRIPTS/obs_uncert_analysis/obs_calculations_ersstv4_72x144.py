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

hist_start = netCDF4.datetime(1969,12,1)
hist_end = netCDF4.datetime(2000,2,28)

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
model_time = time_variable[:]

time_variable_converted = netCDF4.num2date(time_variable[:], time_variable.units, calendar='standard')
#timespan = model_time[numpy.where((model_time[:]>=date_start)&(model_time<date_end))]
#sst_hist_data = sst_hist_data[numpy.where((model_time[:]>=date_start)&(model_time<date_end))[0], :,:]
ncfile.close()

if season=='djf':		
	time_indices = numpy.array([(t.month in [12,1,2])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	sst_hist_data_seas = sst_hist_data[time_indices,:,:]

elif season=='mam':
	time_indices = numpy.array([(t.month in [3,4,5])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	sst_hist_data_seas = sst_hist_data[time_indices,:,:]

elif season=='jja':
	time_indices = numpy.array([(t.month in [6,7,8])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	sst_hist_data_seas = sst_hist_data[time_indices,:,:]
	
elif season=='son':
	time_indices = numpy.array([(t.month in [9,10,11])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	sst_hist_data_seas = sst_hist_data[time_indices,:,:]

elif season=='annual':
	sst_hist_data_seas = sst_hist_data[:,:,:]

sst_histclim_data = numpy.mean(sst_hist_data_seas, axis=0)

# save hist clim
#filename = '/ninod/baird/cmip5/obs_calculations/pr/' + 'obs_GPCP_96x144_sst_1979-09_climatology_' + season + '.nc'
filename = '/ninod/baird/cmip5/obs_calculations/ERSSTv4/' + 'obs_ERSSTv4_72x144_SST_1970-2000_climatology_' + season + '.nc'
ncfile = Dataset(filename, 'w', format='NETCDF4')
lat = ncfile.createDimension('lat', global_nlat)
lon = ncfile.createDimension('lon', global_nlon)	
lats = ncfile.createVariable('lat', 'f4', ('lat',))
lons = ncfile.createVariable('lon', 'f4', ('lon',))
sst_data = ncfile.createVariable('PRECT', 'f4', ('lat','lon'))
lats[:] = global_lat_vals
lons[:] = global_lon_vals
sst_data[:] = sst_histclim_data
lats.units = 'degrees_north'
lons.units = 'degrees_east'
sst_data.units = 'deg C'
ncfile.history = 'Created ' + time.ctime(time.time())
ncfile.close()
print(filename, "saved")