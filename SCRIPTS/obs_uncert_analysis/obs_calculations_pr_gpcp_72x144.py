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

#season='djf'
#season='mam'
#season='jja'
#season='son'
season='annual'

ncfile = Dataset('/ninod/baird/cmip5/obs/pr_gpcp/precip.mon.mean_197901-201507_invertlat.nc', 'r', format='NETCDF4')
pr_data_orig = ncfile.variables['precip'][:]
pr_lat = ncfile.variables['lat'][:]
pr_lon = ncfile.variables['lon'][:]

global_nlat, global_nlon = pr_data_orig.shape[1:3]
global_lat_vals = pr_lat[:]
global_lon_vals = pr_lon[:]

##########################################################################################
##########################################################################################
##########################################################################################

# (YYYY,MM,DD)
#hist_start = netCDF4.datetime(1979,1,1)
#hist_end = netCDF4.datetime(2009,12,31)

#hist_start = netCDF4.datetime(1986,1,1)
#hist_end = netCDF4.datetime(2005,12,31)

hist_start = netCDF4.datetime(1970,1,1)
hist_end = netCDF4.datetime(2000,12,31)

pr_histmonths_list = []
pr_histclim_list = []
hist_stdevs_list = []

# OPEN HISTORICAL PERIOD PR DATA

file_name = '/ninod/baird/cmip5/obs/pr_gpcp/precip.mon.mean_197901-201507_invertlat.nc'

ncfile = Dataset(file_name, 'r', format='NETCDF4')
pr_hist_data = ncfile.variables['precip'][:,:,:]
time_variable = ncfile.variables['time']
date_start = netCDF4.date2num(hist_start, time_variable.units, time_variable.calendar)
date_end = netCDF4.date2num(hist_end, time_variable.units, time_variable.calendar)
model_time = time_variable[:]

time_variable_converted = netCDF4.num2date(time_variable[:], time_variable.units, calendar='standard')
#timespan = model_time[numpy.where((model_time[:]>=date_start)&(model_time<date_end))]
#pr_hist_data = pr_hist_data[numpy.where((model_time[:]>=date_start)&(model_time<date_end))[0], :,:]
ncfile.close()

if season=='djf':		
	time_indices = numpy.array([(t.month in [12,1,2])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	pr_hist_data_seas = pr_hist_data[time_indices,:,:]

elif season=='mam':
	time_indices = numpy.array([(t.month in [3,4,5])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	pr_hist_data_seas = pr_hist_data[time_indices,:,:]

elif season=='jja':
	time_indices = numpy.array([(t.month in [6,7,8])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	pr_hist_data_seas = pr_hist_data[time_indices,:,:]
	
elif season=='son':
	time_indices = numpy.array([(t.month in [9,10,11])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	pr_hist_data_seas = pr_hist_data[time_indices,:,:]

elif season=='annual':
	pr_hist_data_seas = pr_hist_data[:,:,:]

pr_histclim_data = numpy.mean(pr_hist_data_seas, axis=0)

# save hist clim
#filename = '/ninod/baird/cmip5/obs_calculations/pr/' + 'obs_GPCP_96x144_pr_1979-09_climatology_' + season + '.nc'
filename = '/ninod/baird/cmip5/obs_calculations/pr_gpcp/' + 'obs_GPCP_72x144_PRECT_1970-2000_climatology_' + season + '.nc'
ncfile = Dataset(filename, 'w', format='NETCDF4')
lat = ncfile.createDimension('lat', global_nlat)
lon = ncfile.createDimension('lon', global_nlon)	
lats = ncfile.createVariable('lat', 'f4', ('lat',))
lons = ncfile.createVariable('lon', 'f4', ('lon',))
pr_data = ncfile.createVariable('PRECT', 'f4', ('lat','lon'))
lats[:] = global_lat_vals
lons[:] = global_lon_vals
pr_data[:] = pr_histclim_data
lats.units = 'degrees_north'
lons.units = 'degrees_east'
pr_data.units = 'none'
ncfile.history = 'Created ' + time.ctime(time.time())
ncfile.close()
print(filename, "saved")