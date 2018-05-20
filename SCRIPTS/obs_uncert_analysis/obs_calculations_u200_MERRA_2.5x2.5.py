#!/home/baird/anaconda3_new/bin/python
###!/Users/baird/anaconda/bin/python

import numpy
from netCDF4 import Dataset
import netCDF4
import datetime
import time

##########################################################################################
##########################################################################################
##########################################################################################

season_names = ['djf','mam','jja','son','ondjfm']

season='djf'
#season='jja'
#season='annual'

##########################################################################################
##########################################################################################
##########################################################################################

# (YYYY,MM,DD)
#hist_start = netCDF4.datetime(1979,1,1)
#hist_end = netCDF4.datetime(2009,12,31)

hist_start = netCDF4.datetime(1970,1,1)
hist_end = netCDF4.datetime(2000,12,31)

u200_histmonths_list = []
u200_histclim_list = []
hist_stdevs_list = []

#ncfile = Dataset('/Users/baird/google_drive/_data_original/era_interim/u200_MERRA_monthly_197901-201412_invertlat_96x144regrid.nc', 'r', format='NETCDF4')
ncfile = Dataset('/ninod/baird/cmip5/obs/MERRA/u200/MERRA_197901-201412_u200_monthly_2.5x2.5regrid.nc', 'r', format='NETCDF4')
u200_hist_data = ncfile.variables['u200'][:]
u200_lat = ncfile.variables['lat'][:]
u200_lon = ncfile.variables['lon'][:]

global_nlat, global_nlon = u200_hist_data.shape[1:3]
global_lat_vals = u200_lat[:]
global_lon_vals = u200_lon[:]

time_variable = ncfile.variables['time']
date_start = netCDF4.date2num(hist_start, time_variable.units, time_variable.calendar)
date_end = netCDF4.date2num(hist_end, time_variable.units, time_variable.calendar)
model_time = time_variable[:]

#timespan = model_time[numpy.where((model_time[:]>=date_start)&(model_time<date_end))]
#u200_hist_data = u200_hist_data[numpy.where((model_time[:]>=date_start)&(model_time<date_end))[0], :,:]
#ncfile.close()

time_variable_converted = netCDF4.num2date(time_variable[:], time_variable.units, time_variable.calendar)

if season=='djf':		
	time_indices = numpy.array([(t.month in [12,1,2])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	u200_hist_data_seas = u200_hist_data[time_indices,:,:]

elif season=='jja':
	time_indices = numpy.array([(t.month in [6,7,8])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	u200_hist_data_seas = u200_hist_data[time_indices,:,:]

elif season=='annual':
	u200_hist_data_seas = u200_hist_data[:,:,:]

u200_histclim_data = numpy.mean(u200_hist_data_seas, axis=0)

# save hist clim

#filename = '/Users/baird/google_drive/_data_analyzed/obs_calculations/u200_MERRA/' + 'obs_MERRA_96x144_u200_1986-05_climatology_' + season + '.nc'
filename = '/ninod/baird/cmip5/obs_calculations/u200_MERRA/' + 'obs_MERRA_2.5x2.5_u200_1970-00_climatology_' + season + '.nc'
ncfile = Dataset(filename, 'w', format='NETCDF4')
lat = ncfile.createDimension('lat', global_nlat)
lon = ncfile.createDimension('lon', global_nlon)	
lats = ncfile.createVariable('lat', 'f4', ('lat',))
lons = ncfile.createVariable('lon', 'f4', ('lon',))
u200_data = ncfile.createVariable('u200', 'f4', ('lat','lon'))
lats[:] = global_lat_vals
lons[:] = global_lon_vals
u200_data[:] = u200_histclim_data
lats.units = 'degrees_north'
lons.units = 'degrees_east'
u200_data.units = 'm s-1'
ncfile.history = 'Created ' + time.ctime(time.time())
ncfile.close()
print(filename, "saved")