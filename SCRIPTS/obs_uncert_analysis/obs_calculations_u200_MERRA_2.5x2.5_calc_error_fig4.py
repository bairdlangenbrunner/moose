import numpy
import netCDF4
import datetime
import time
import scipy.stats

##########################################################################################
##########################################################################################
##########################################################################################

season_names = ['djf','mam','jja','son','ondjfm']

season='djf'
#season='jja'
#season='annual'

# ua_corr_lat_lo, ua_corr_lat_hi, ua_corr_lon_lo, ua_corr_lon_hi = 30,40,170.,205.
# ts_corr_lat_lo, ts_corr_lat_hi, ts_corr_lon_lo, ts_corr_lon_hi = -25,-15,200.,230.

##########################################################################################
##########################################################################################
##########################################################################################

# (YYYY,MM,DD)

hist_start = netCDF4.datetime(1979,12,1)
hist_end = netCDF4.datetime(2010,2,28)

sst_histmonths_list = []
sst_histclim_list = []
hist_stdevs_list = []

#ncfile = Dataset('/Users/baird/google_drive/_data_original/era_interim/u200_MERRA_monthly_197901-201412_invertlat_96x144regrid.nc', 'r', format='NETCDF4')
ncfile = netCDF4.Dataset('/ninod/baird/cmip5/obs/MERRA/u200/MERRA_197901-201412_u200_monthly_2.5x2.5regrid.nc', 'r', format='NETCDF4')
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

ua_lat_lo, ua_lat_hi, ua_lon_lo, ua_lon_hi = 30,40,170.,205.; region = 'correlation_region'

# pull out only lat lon that are relevant
u200_lat_inds = numpy.where((global_lat_vals>=ua_lat_lo) & (global_lat_vals<=ua_lat_hi))[0]
u200_lon_inds = numpy.where((global_lon_vals>=ua_lon_lo) & (global_lon_vals<=ua_lon_hi))[0]
u200_hist_data = u200_hist_data[:, u200_lat_inds[0]:(u200_lat_inds[-1]+1), u200_lon_inds[0]:(u200_lon_inds[-1]+1)]

if season=='djf':		
	time_indices = numpy.array([(t.month in [12,1,2])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	u200_hist_data_seas = u200_hist_data[time_indices,:,:]

regional_nlat, regional_nlon = u200_hist_data_seas.shape[1:3]
u200_histclim_data = numpy.mean(u200_hist_data_seas, axis=0)
u200_histclim_stdev = numpy.std(u200_hist_data_seas, axis=0, ddof=1)
u200_hist_data_seas_means = u200_hist_data_seas.reshape((-1,3,regional_nlat,regional_nlon)).mean(axis=1)

print(u200_hist_data_seas_means.shape)

u200_hist_data_seas_spatial_means = numpy.mean(u200_hist_data_seas_means, axis=2)
u200_hist_data_seas_spatial_means = numpy.mean(u200_hist_data_seas_spatial_means, axis=1)

print(u200_hist_data_seas_spatial_means.shape)

N = u200_hist_data_seas_spatial_means.size
u200_stdev = numpy.std(u200_hist_data_seas_spatial_means)
std_err = u200_stdev/numpy.sqrt(N)
print(std_err)