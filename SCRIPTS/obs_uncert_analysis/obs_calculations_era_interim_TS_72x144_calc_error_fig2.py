import numpy
from netCDF4 import Dataset
import netCDF4
import datetime
import time
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as mp

##########################################################################################
##########################################################################################
##########################################################################################

season='djf'

##########################################################################################
##########################################################################################
##########################################################################################

hist_start = netCDF4.datetime(1969,12,1)
hist_end = netCDF4.datetime(2000,2,28)

skt_histmonths_list = []
skt_histclim_list = []
hist_stdevs_list = []

ncfile = Dataset('/ninod/baird/cmip5/obs/era_interim/skt_erai_monthly_197901-201412_invertlat_72x144regrid.nc', 'r', format='NETCDF4')
skt_full = ncfile.variables['skt']
skt_hist_data = ncfile.variables['skt'][:]# * skt_scale_factor) + skt_add_offset
skt_lat = ncfile.variables['lat'][:]
skt_lon = ncfile.variables['lon'][:]

global_nlat, global_nlon = skt_hist_data.shape[1:3]
global_lat_vals = skt_lat[:]
global_lon_vals = skt_lon[:]

time_variable = ncfile.variables['time']
date_start = netCDF4.date2num(hist_start, time_variable.units, time_variable.calendar)
date_end = netCDF4.date2num(hist_end, time_variable.units, time_variable.calendar)
time_variable_converted = netCDF4.num2date(time_variable[:], time_variable.units, calendar='standard')

#pr_lat_lo, pr_lat_hi, pr_lon_lo, pr_lon_hi = 30., 45., 232.5, 248; region = 'CA'
ts_lat_lo, ts_lat_hi, ts_lon_lo, ts_lon_hi = -30., 10., 155., 270.; region = 'tropacific'
#ua_lat_lo, ua_lat_hi, ua_lon_lo, ua_lon_hi = 20., 50., 170., 250.; region = 'midlatpacific'

# pull out only lat lon that are relevant
skt_lat_inds = numpy.where((global_lat_vals>=ts_lat_lo) & (global_lat_vals<=ts_lat_hi))[0]
skt_lon_inds = numpy.where((global_lon_vals>=ts_lon_lo) & (global_lon_vals<=ts_lon_hi))[0]
skt_hist_data = skt_hist_data[:, skt_lat_inds[0]:(skt_lat_inds[-1]+1), skt_lon_inds[0]:(skt_lon_inds[-1]+1)]

if season=='djf':		
	time_indices = numpy.array([(t.month in [12,1,2])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	skt_hist_data_seas = skt_hist_data[time_indices,:,:]

regional_nlat, regional_nlon = skt_hist_data_seas.shape[1:3]

skt_histclim_data = numpy.mean(skt_hist_data_seas, axis=0)
skt_histclim_stdev = numpy.std(skt_hist_data_seas, axis=0, ddof=1)

skt_hist_data_seas_means = skt_hist_data_seas.reshape((-1,3,regional_nlat,regional_nlon)).mean(axis=1)

print(skt_hist_data_seas_means.shape)

fig = mp.figure(figsize=(5,10))
ax = fig.add_subplot(211)
c=ax.contourf(skt_hist_data_seas_means[0,:,:])
mp.colorbar(c)
ax = fig.add_subplot(212)
c=ax.contourf(skt_histclim_data)
mp.colorbar(c)

mp.show()


exit()
filename = '/ninod/baird/cmip5/obs_calculations/ppe_obs_comparisons/' + 'obs_era_interim_72x144_skt_1970-2000_climatology_' + season + '.nc'
ncfile = Dataset(filename, 'w', format='NETCDF4')
lat = ncfile.createDimension('lat', global_nlat)
lon = ncfile.createDimension('lon', global_nlon)	
lats = ncfile.createVariable('lat', 'f4', ('lat',))
lons = ncfile.createVariable('lon', 'f4', ('lon',))
skt_data = ncfile.createVariable('skt', 'f4', ('lat','lon'))
lats[:] = global_lat_vals
lons[:] = global_lon_vals
skt_data[:] = skt_histclim_data
lats.units = 'degrees_north'
lons.units = 'degrees_east'
skt_data.units = skt_full.units
ncfile.history = 'Created ' + time.ctime(time.time())
ncfile.close()
print(filename, "saved")