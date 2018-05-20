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
season='jja'
season='annual'

##########################################################################################
##########################################################################################
##########################################################################################

# (YYYY,MM,DD)
#hist_start = netCDF4.datetime(1979,1,1)
#hist_end = netCDF4.datetime(2009,12,31)

hist_start = netCDF4.datetime(1970,1,1)
hist_end = netCDF4.datetime(2000,12,31)

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
model_time = time_variable[:]

#timespan = model_time[numpy.where((model_time[:]>=date_start)&(model_time<date_end))]
skt_hist_data = skt_hist_data[numpy.where((model_time[:]>=date_start)&(model_time<date_end))[0], :,:]

print(skt_hist_data.shape)



if season=='djf':
	skt_hist_data_dec = skt_hist_data[11::12,:,:] # make sure these indices correspond to the start_date month
	skt_hist_data_jan = skt_hist_data[12::12,:,:] # dec WAS 11, 
	skt_hist_data_feb = skt_hist_data[13::12,:,:]
	skt_hist_data_seas = numpy.zeros((skt_hist_data_dec.shape[0]+skt_hist_data_jan.shape[0]+skt_hist_data_feb.shape[0], global_nlat, global_nlon))
	skt_hist_data_seas[0::3,:,:] = skt_hist_data_dec
	skt_hist_data_seas[1::3,:,:] = skt_hist_data_jan
	skt_hist_data_seas[2::3,:,:] = skt_hist_data_feb

elif season=='jja':
	skt_hist_data_jun = skt_hist_data[5::12,:,:]
	skt_hist_data_jul = skt_hist_data[6::12,:,:]
	skt_hist_data_aug = skt_hist_data[7::12,:,:]
	skt_hist_data_seas = numpy.zeros((skt_hist_data_jun.shape[0]+skt_hist_data_jul.shape[0]+skt_hist_data_aug.shape[0], global_nlat, global_nlon))
	skt_hist_data_seas[0::3,:,:] = skt_hist_data_jun
	skt_hist_data_seas[1::3,:,:] = skt_hist_data_jul
	skt_hist_data_seas[2::3,:,:] = skt_hist_data_aug

elif season=='annual':
	skt_hist_data_seas = skt_hist_data[:,:,:]

skt_histclim_data = numpy.mean(skt_hist_data_seas, axis=0)

skt_histclim_stdev = numpy.std(skt_hist_data_seas, axis=0, ddof=1)

# filename = '/ninod/baird/cmip5/obs_calculations/ppe_obs_comparisons/' + 'obs_era_interim_72x144_skt_1960-1990_stdev_' + season + '.nc'
# ncfile = Dataset(filename, 'w', format='NETCDF4')
# lat = ncfile.createDimension('lat', global_nlat)
# lon = ncfile.createDimension('lon', global_nlon)	
# lats = ncfile.createVariable('lat', 'f4', ('lat',))
# lons = ncfile.createVariable('lon', 'f4', ('lon',))
# skt_data = ncfile.createVariable('skt', 'f4', ('lat','lon'))
# lats[:] = global_lat_vals
# lons[:] = global_lon_vals
# skt_data[:] = skt_histclim_stdev
# lats.units = 'degrees_north'
# lons.units = 'degrees_east'
# skt_data.units = skt_full.units
# ncfile.history = 'Created ' + time.ctime(time.time())
# ncfile.close()
# print(filename, "saved")

#cf = mp.contourf(skt_lon, skt_lat, skt_histclim_data)
#cbar = mp.colorbar()
#mp.show()
#exit()

# save hist clim

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