import numpy
import scipy.stats
import sys
import netCDF4
import datetime
import time
import scipy.stats

##########################################################################################
##########################################################################################
##########################################################################################

season_names = ['djf','mam','jja','son','ondjfm']

season='djf'
#season='mam'
#season='jja'
#season='son'
#season='annual'

ncfile = netCDF4.Dataset('/ninod/baird/cmip5/obs/pr_gpcp/precip.mon.mean_197901-201507_invertlat.nc', 'r', format='NETCDF4')
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

#hist_start = netCDF4.datetime(1969,12,1)
#hist_end = netCDF4.datetime(2000,2,28)

hist_start = netCDF4.datetime(1979,12,1)
hist_end = netCDF4.datetime(2010,2,28)

sst_histmonths_list = []
sst_histclim_list = []
hist_stdevs_list = []

# OPEN HISTORICAL PERIOD PR DATA

file_name = '/ninod/baird/cmip5/obs/pr_gpcp/precip.mon.mean_197901-201507_invertlat.nc'

ncfile = netCDF4.Dataset(file_name, 'r', format='NETCDF4')
pr_hist_data = ncfile.variables['precip'][:,:,:]
time_variable = ncfile.variables['time']
date_start = netCDF4.date2num(hist_start, time_variable.units, time_variable.calendar)
date_end = netCDF4.date2num(hist_end, time_variable.units, time_variable.calendar)
model_time = time_variable[:]
time_variable_converted = netCDF4.num2date(time_variable[:], time_variable.units, calendar='standard')

pr_lat_lo, pr_lat_hi, pr_lon_lo, pr_lon_hi = 30., 45., 232.5, 248; region = 'CA'
#pr_lat_lo, pr_lat_hi, pr_lon_lo, pr_lon_hi = -30., 10., 155., 270.; region = 'tropacific'
#ua_lat_lo, ua_lat_hi, ua_lon_lo, ua_lon_hi = 20., 50., 170., 250.; region = 'midlatpacific'

# pull out only lat lon that are relevant
pr_lat_inds = numpy.where((global_lat_vals>=pr_lat_lo) & (global_lat_vals<=pr_lat_hi))[0]
pr_lon_inds = numpy.where((global_lon_vals>=pr_lon_lo) & (global_lon_vals<=pr_lon_hi))[0]
pr_hist_data = pr_hist_data[:, pr_lat_inds[0]:(pr_lat_inds[-1]+1), pr_lon_inds[0]:(pr_lon_inds[-1]+1)]

if season=='djf':		
	time_indices = numpy.array([(t.month in [12,1,2])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	pr_hist_data_seas = pr_hist_data[time_indices,:,:]

regional_nlat, regional_nlon = pr_hist_data_seas.shape[1:3]
pr_histclim_data = numpy.mean(pr_hist_data_seas, axis=0)
pr_histclim_stdev = numpy.std(pr_hist_data_seas, axis=0, ddof=1)
pr_hist_data_seas_means = pr_hist_data_seas.reshape((-1,3,regional_nlat,regional_nlon)).mean(axis=1)

print(pr_hist_data_seas_means.shape)
print(pr_hist_data_seas.shape)

seas_rmse_values = numpy.zeros((pr_hist_data_seas_means.shape[0]))
for seas_clim_idx in range(pr_hist_data_seas_means.shape[0]):
	seas_rmse_values[seas_clim_idx] = numpy.sqrt(numpy.mean( (pr_hist_data_seas_means[seas_clim_idx,:,:]-pr_histclim_data)**2. ) )

rval, pval = scipy.stats.pearsonr(seas_rmse_values[1:], seas_rmse_values[:-1])

N_full = pr_hist_data_seas_means.shape[0]
N_adjusted = N_full*((1-rval)/(1+rval))

print( 1./numpy.sqrt(N_adjusted) * numpy.std(seas_rmse_values) )
print(numpy.std(seas_rmse_values))

#  0.0816906403487 mm day^-1