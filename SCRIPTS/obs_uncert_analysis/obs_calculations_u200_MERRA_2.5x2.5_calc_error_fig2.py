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

#pr_lat_lo, pr_lat_hi, pr_lon_lo, pr_lon_hi = 30., 45., 232.5, 248; region = 'CA'
#ts_lat_lo, ts_lat_hi, ts_lon_lo, ts_lon_hi = -30., 10., 155., 270.; region = 'tropacific'
ua_lat_lo, ua_lat_hi, ua_lon_lo, ua_lon_hi = 20., 50., 170., 250.; region = 'midlatpacific'

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
print(u200_hist_data_seas.shape)

seas_rmse_values = numpy.zeros((u200_hist_data_seas_means.shape[0]))
for seas_clim_idx in range(u200_hist_data_seas_means.shape[0]):
	seas_rmse_values[seas_clim_idx] = numpy.sqrt(numpy.mean( (u200_hist_data_seas_means[seas_clim_idx,:,:]-u200_histclim_data)**2. ) )

rval, pval = scipy.stats.pearsonr(seas_rmse_values[1:], seas_rmse_values[:-1])

N_full = u200_hist_data_seas_means.shape[0]
N_adjusted = N_full*((1-rval)/(1+rval))

print( 1./numpy.sqrt(N_adjusted) * numpy.std(seas_rmse_values) )
print(numpy.std(seas_rmse_values))

#  1.54937455102 m s-1