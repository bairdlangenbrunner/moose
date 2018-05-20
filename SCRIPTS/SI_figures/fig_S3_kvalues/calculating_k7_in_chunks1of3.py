import numpy
import netCDF4
#import matplotlib.pyplot as mp
#import matplotlib.colors as mc
#import matplotlib.cm as cm
#import mpl_toolkits.mplot3d
#import matplotlib
#import scipy.ndimage
import datetime

import itertools
import random
import numpy.random
import scipy.stats
import os

model_names = numpy.array(( 'ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1-m', 'bcc-csm1-1', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5', 'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'EC-EARTH', 'FGOALS-g2', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M', 'GISS-E2-H', 'GISS-E2-R', 'HadGEM2-AO', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'MIROC5', 'MIROC-ESM-CHEM', 'MIROC-ESM', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'NorESM1-ME', 'NorESM1-M' ))
nmods = len(model_names)


# In[3]:

season_names = ["djf","mam","jja","son","annual","ondjfm"]

pr_lat_lo, pr_lat_hi, pr_lon_lo, pr_lon_hi = 30., 45., 232.5, 248; region = 'CA'

ts_lat_lo, ts_lat_hi, ts_lon_lo, ts_lon_hi = -30., 10., 155., 270.; region = 'tropacific'

ua_lat_lo, ua_lat_hi, ua_lon_lo, ua_lon_hi = 20., 50., 170., 250.; region = 'midlatpacific'

tel_lat_lo, tel_lat_hi, tel_lon_lo, tel_lon_hi = 25., 60., 220., 255.; region = 'CA'

season='djf'; SEASON='DJF'
#season='jja'; SEASON='JJA'
#season='annual'; SEASON='annual'


# In[4]:

# OPEN TS DATASET
ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cmip5_calculations/ts/djf/NorESM1-M_ts_1980-10_climatology_djf.nc', 'r', format='NETCDF4')

ts_data_orig = ncfile.variables['ts'][:]
ts_lat = ncfile.variables['lat'][:]
ts_lon = ncfile.variables['lon'][:]

# pull out lat/lon indices
ts_lat_inds = numpy.where((ts_lat>=ts_lat_lo) & (ts_lat<=ts_lat_hi))[0]
ts_lon_inds = numpy.where((ts_lon>=ts_lon_lo) & (ts_lon<=ts_lon_hi))[0]
ts_regional_lat_vals = ts_lat[ts_lat_inds[0]:(ts_lat_inds[-1]+1)]
ts_regional_lon_vals = ts_lon[ts_lon_inds[0]:(ts_lon_inds[-1]+1)]    

ts_data = ts_data_orig[ts_lat_inds[0]:(ts_lat_inds[-1]+1), ts_lon_inds[0]:(ts_lon_inds[-1]+1)]
ts_regional_nlat, ts_regional_nlon = ts_data.shape
global_nlat, global_nlon = ts_data_orig.shape[0:2]
global_lat_vals = ts_lat[:]
global_lon_vals = ts_lon[:]


# In[5]:

# OPEN TS OBSERVATIONS
# OPEN TS OBSERVATIONS
ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/obs_calculations/ERSSTv4/obs_ERSSTv4_72x144_SST_1980-2010_climatology_'+season+'.nc', 'r', format='NETCDF4')
obs_field_ts = ncfile.variables['sst'][ts_lat_inds[0]:(ts_lat_inds[-1]+1), ts_lon_inds[0]:(ts_lon_inds[-1]+1)]+273.15

# OPEN PR OBSERVATIONS
ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/obs_calculations/pr_gpcp/obs_GPCP_72x144_PRECT_1980-2010_climatology_'+season+'.nc', 'r', format='NETCDF4')
pr_lat = ncfile.variables['lat'][:]
pr_lon = ncfile.variables['lon'][:]
pr_lat_inds = numpy.where((pr_lat>=pr_lat_lo) & (pr_lat<=pr_lat_hi))[0]
pr_lon_inds = numpy.where((pr_lon>=pr_lon_lo) & (pr_lon<=pr_lon_hi))[0]
obs_field_pr = ncfile.variables['PRECT'][pr_lat_inds[0]:(pr_lat_inds[-1]+1), pr_lon_inds[0]:(pr_lon_inds[-1]+1)]
pr_regional_nlat, pr_regional_nlon = obs_field_pr.shape

pr_regional_lat_vals = pr_lat[pr_lat_inds[0]:(pr_lat_inds[-1]+1)]
pr_regional_lon_vals = pr_lon[pr_lon_inds[0]:(pr_lon_inds[-1]+1)] 

# OPEN UA OBSERVATIONS
ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/obs_calculations/u200_MERRA/obs_MERRA_2.5x2.5_u200_1980-10_climatology_'+season+'.nc', 'r', format='NETCDF4')
ua_lat = ncfile.variables['lat'][:]
ua_lon = ncfile.variables['lon'][:]
ua_lat_inds = numpy.where((ua_lat>=ua_lat_lo) & (ua_lat<=ua_lat_hi))[0]
ua_lon_inds = numpy.where((ua_lon>=ua_lon_lo) & (ua_lon<=ua_lon_hi))[0]
obs_field_ua = ncfile.variables['u200'][ua_lat_inds[0]:(ua_lat_inds[-1]+1), ua_lon_inds[0]:(ua_lon_inds[-1]+1)]
ua_regional_nlat, ua_regional_nlon = obs_field_ua.shape

ua_regional_lat_vals = ua_lat[ua_lat_inds[0]:(ua_lat_inds[-1]+1)]
ua_regional_lon_vals = ua_lon[ua_lon_inds[0]:(ua_lon_inds[-1]+1)] 

# OPEN PR TELECONNECTIONS OBSERVATIONS
ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/obs_calculations/pr_teleconnections/pr_teleconnections_GPCP_ERSSTv4_nino34_1979-2010_djf.nc', 'r', format='NETCDF4')
tel_lat = ncfile.variables['lat'][:]
tel_lon = ncfile.variables['lon'][:]
obs_field_tel = ncfile.variables['pr'][(tel_lat>=tel_lat_lo)&(tel_lat<=tel_lat_hi),:][:,(tel_lon>=tel_lon_lo)&(tel_lon<=tel_lon_hi)]
tel_regional_lat_vals = tel_lat[(tel_lat>=tel_lat_lo)&(tel_lat<=tel_lat_hi)]
tel_regional_lon_vals = tel_lon[(tel_lon>=tel_lon_lo)&(tel_lon<=tel_lon_hi)]
tel_regional_nlat, tel_regional_nlon = obs_field_tel.shape


# In[6]:

# set up data
model_data_hist_pr = numpy.zeros((len(model_names), pr_regional_nlat, pr_regional_nlon))
model_data_eoc_pr = numpy.zeros((len(model_names), pr_regional_nlat, pr_regional_nlon))
model_data_hist_pr_LENS = numpy.zeros((40, pr_regional_nlat, pr_regional_nlon))

for i in range(nmods):
    #print("opening model", model_names[i])
    modelname = model_names[i]
    # OPEN HISTORICAL FIELDS
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cmip5_calculations/pr/'+season+'/'+modelname+'_pr_1980-10_climatology_'+season+'.nc', 'r', format='NETCDF4')
    model_data_hist_pr[i,:,:] = ncfile.variables['pr'][pr_lat_inds[0]:(pr_lat_inds[-1]+1), pr_lon_inds[0]:(pr_lon_inds[-1]+1)]
    ncfile.close()
    # OPEN FUTURE CHANGE FIELDS
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cmip5_calculations/pr/'+season+'/'+modelname+'_pr_2070-99_climatology_'+season+'.nc', 'r', format='NETCDF4')
    model_data_eoc_pr[i,:,:] = ncfile.variables['pr'][pr_lat_inds[0]:(pr_lat_inds[-1]+1), pr_lon_inds[0]:(pr_lon_inds[-1]+1)]
    ncfile.close()

LENS_names = ['{:02d}'.format(i) for i in range(1,36)] + ['{:03d}'.format(i) for i in range(101,106)]
for i in range(len(LENS_names)): # 40
    member_name = LENS_names[i]
    # get convective precipitation
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cesm_LENS_calculations/PRECC/'+season+'/'+member_name + '_PRECC_1980-10_climatology_'+season+'_2.5x2.5regrid.nc', 'r', format='NETCDF4')
    precc_temp = ncfile.variables['PRECC'][pr_lat_inds[0]:(pr_lat_inds[-1]+1), pr_lon_inds[0]:(pr_lon_inds[-1]+1)]
    # get large-scale precipitation
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cesm_LENS_calculations/PRECL/'+season+'/'+member_name + '_PRECL_1980-10_climatology_'+season+'_2.5x2.5regrid.nc', 'r', format='NETCDF4')
    precl_temp = ncfile.variables['PRECL'][pr_lat_inds[0]:(pr_lat_inds[-1]+1), pr_lon_inds[0]:(pr_lon_inds[-1]+1)]
    # add together
    model_data_hist_pr_LENS[i,:,:] = precc_temp + precl_temp


# In[7]:

# IMPORT TS DATA
model_data_hist_ts = numpy.zeros((len(model_names), ts_regional_nlat, ts_regional_nlon))
model_data_eoc_ts = numpy.zeros((len(model_names), ts_regional_nlat, ts_regional_nlon))
model_data_hist_ts_LENS = numpy.zeros((40, ts_regional_nlat, ts_regional_nlon))

for i in range(nmods):
    #print("opening model", model_names[i])
    modelname = model_names[i]
    # OPEN HISTORICAL FIELDS
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cmip5_calculations/ts/'+season+'/'+modelname+'_ts_1980-10_climatology_'+season+'.nc', 'r', format='NETCDF4')
    model_data_hist_ts[i,:,:] = ncfile.variables['ts'][ts_lat_inds[0]:(ts_lat_inds[-1]+1), ts_lon_inds[0]:(ts_lon_inds[-1]+1)]
    ncfile.close()
    # OPEN FUTURE CHANGE FIELDS
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cmip5_calculations/ts/'+season+'/'+modelname+'_ts_2070-99_climatology_'+season+'.nc', 'r', format='NETCDF4')
    model_data_eoc_ts[i,:,:] = ncfile.variables['ts'][ts_lat_inds[0]:(ts_lat_inds[-1]+1), ts_lon_inds[0]:(ts_lon_inds[-1]+1)]
    ncfile.close()

LENS_names = ['{:02d}'.format(i) for i in range(1,36)] + ['{:03d}'.format(i) for i in range(101,106)]
for i in range(len(LENS_names)): # 40
    member_name = LENS_names[i]
    # get convective precipitation
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cesm_LENS_calculations/TS/'+season+'/'+member_name + '_TS_1980-10_climatology_'+season+'_2.5x2.5regrid.nc', 'r', format='NETCDF4')
    model_data_hist_ts_LENS[i,:,:] = ncfile.variables['TS'][ts_lat_inds[0]:(ts_lat_inds[-1]+1), ts_lon_inds[0]:(ts_lon_inds[-1]+1)]    


# In[ ]:




# In[8]:

# IMPORT UA200 DATA
model_data_hist_ua = numpy.zeros((len(model_names), ua_regional_nlat, ua_regional_nlon))
model_data_eoc_ua = numpy.zeros((len(model_names), ua_regional_nlat, ua_regional_nlon))
model_data_hist_ua_LENS = numpy.zeros((40, ua_regional_nlat, ua_regional_nlon))

for i in range(nmods):
    #print("opening model", model_names[i])
    modelname = model_names[i]
    # OPEN HISTORICAL FIELDS
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cmip5_calculations/ua200/'+season+'/'+modelname+'_ua200_1980-10_climatology_'+season+'.nc', 'r', format='NETCDF4')
    model_data_hist_ua[i,:,:] = ncfile.variables['ua'][ua_lat_inds[0]:(ua_lat_inds[-1]+1), ua_lon_inds[0]:(ua_lon_inds[-1]+1)]
    ncfile.close()
    # OPEN FUTURE CHANGE FIELDS
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cmip5_calculations/ua200/'+season+'/'+modelname+'_ua200_2070-99_climatology_'+season+'.nc', 'r', format='NETCDF4')
    model_data_eoc_ua[i,:,:] = ncfile.variables['ua'][ua_lat_inds[0]:(ua_lat_inds[-1]+1), ua_lon_inds[0]:(ua_lon_inds[-1]+1)]
    ncfile.close()

LENS_names = ['{:02d}'.format(i) for i in range(1,36)] + ['{:03d}'.format(i) for i in range(101,106)]
for i in range(len(LENS_names)): # 40
    member_name = LENS_names[i]
    # get convective precipitation
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cesm_LENS_calculations/U200/'+season+'/'+member_name + '_U_1980-10_climatology_'+season+'_2.5x2.5regrid.nc', 'r', format='NETCDF4')
    model_data_hist_ua_LENS[i,:,:] = ncfile.variables['U'][ua_lat_inds[0]:(ua_lat_inds[-1]+1), ua_lon_inds[0]:(ua_lon_inds[-1]+1)]    


# In[9]:

# IMPORT TELECONNECTIONS DATA
model_data_hist_tel = numpy.zeros((len(model_names), tel_regional_nlat, tel_regional_nlon))

for i in range(nmods):
    #print("opening model", model_names[i])
    modelname = model_names[i]
    # OPEN HISTORICAL FIELDS
    ncfile = netCDF4.Dataset('/Volumes/data1/baird/subensembles_data/cmip5_calculations/pr_teleconnections/'+modelname+'_pr_teleconnections_nino34_1961-90_'+season+'.nc', 'r', format='NETCDF4')
    model_data_hist_tel[i,:,:] = ncfile.variables['pr'][(tel_lat>=tel_lat_lo)&(tel_lat<=tel_lat_hi),:][:,(tel_lon>=tel_lon_lo)&(tel_lon<=tel_lon_hi)]
    ncfile.close()


# In[10]:

ncal_latlon = numpy.load('/Volumes/data1/baird/subensembles_data/ncal_latlon_array.npy')
ccal_latlon = numpy.load('/Volumes/data1/baird/subensembles_data/ccal_latlon_array.npy')
scal_latlon = numpy.load('/Volumes/data1/baird/subensembles_data/scal_latlon_array.npy')
coastal_cal_latlon = numpy.load('/Volumes/data1/baird/subensembles_data/coastal_cal_latlon_array.npy')


# In[11]:

# precip indices
pr_indices_lon_ncal = [ numpy.where(numpy.in1d(pr_regional_lon_vals, ncal_latlon[i,0]))[0][0] for i in range(ncal_latlon.shape[0]) ]
pr_indices_lat_ncal = [ numpy.where(numpy.in1d(pr_regional_lat_vals, ncal_latlon[i,1]))[0][0] for i in range(ncal_latlon.shape[0]) ]

pr_indices_lon_ccal = [ numpy.where(numpy.in1d(pr_regional_lon_vals, ccal_latlon[i,0]))[0][0] for i in range(ccal_latlon.shape[0]) ]
pr_indices_lat_ccal = [ numpy.where(numpy.in1d(pr_regional_lat_vals, ccal_latlon[i,1]))[0][0] for i in range(ccal_latlon.shape[0]) ]

pr_indices_lon_scal = [ numpy.where(numpy.in1d(pr_regional_lon_vals, scal_latlon[i,0]))[0][0] for i in range(scal_latlon.shape[0]) ]
pr_indices_lat_scal = [ numpy.where(numpy.in1d(pr_regional_lat_vals, scal_latlon[i,1]))[0][0] for i in range(scal_latlon.shape[0]) ]

pr_indices_lon_coastal_cal = [ numpy.where(numpy.in1d(pr_regional_lon_vals, coastal_cal_latlon[i,0]))[0][0] for i in range(coastal_cal_latlon.shape[0]) ]
pr_indices_lat_coastal_cal = [ numpy.where(numpy.in1d(pr_regional_lat_vals, coastal_cal_latlon[i,1]))[0][0] for i in range(coastal_cal_latlon.shape[0]) ]

# not sure if I'll need teleconnection indices but here they are in case...
tel_indices_lon_ncal = [ numpy.where(numpy.in1d(tel_regional_lon_vals, ncal_latlon[i,0]))[0][0] for i in range(ncal_latlon.shape[0]) ]
tel_indices_lat_ncal = [ numpy.where(numpy.in1d(tel_regional_lat_vals, ncal_latlon[i,1]))[0][0] for i in range(ncal_latlon.shape[0]) ]

tel_indices_lon_ccal = [ numpy.where(numpy.in1d(tel_regional_lon_vals, ccal_latlon[i,0]))[0][0] for i in range(ccal_latlon.shape[0]) ]
tel_indices_lat_ccal = [ numpy.where(numpy.in1d(tel_regional_lat_vals, ccal_latlon[i,1]))[0][0] for i in range(ccal_latlon.shape[0]) ]

tel_indices_lon_scal = [ numpy.where(numpy.in1d(tel_regional_lon_vals, scal_latlon[i,0]))[0][0] for i in range(scal_latlon.shape[0]) ]
tel_indices_lat_scal = [ numpy.where(numpy.in1d(tel_regional_lat_vals, scal_latlon[i,1]))[0][0] for i in range(scal_latlon.shape[0]) ]


# # take all data and ravel

# In[12]:

# NOW TAKE ALL PR DATA AND RAVEL IT
# CALCULATE ENSEMBLE MEAN FOR EOC CONVERGENCE
model_field_mmem_pr = numpy.mean(model_data_eoc_pr, axis=0)
model_field_mmem_ts = numpy.mean(model_data_eoc_ts, axis=0)
model_field_mmem_ua = numpy.mean(model_data_eoc_ua, axis=0)
#model_field_mmem_tel = numpy.mean(model_data_hist_tel, axis=0)

# NOW CALCULATE BIAS AND CONVERGENCE
bias_values_pr = numpy.zeros((nmods))
convergence_values_pr = numpy.zeros((nmods))

bias_values_ts = numpy.zeros((nmods))
convergence_values_ts = numpy.zeros((nmods))

bias_values_ua = numpy.zeros((nmods))
convergence_values_ua = numpy.zeros((nmods))

bias_values_tel = numpy.zeros((nmods))
correlation_vals_tel = numpy.zeros((nmods))

for i in range(nmods):
    hist_field_pr = model_data_hist_pr[i,:,:]
    eoc_field_pr = model_data_eoc_pr[i,:,:]
    
    hist_field_ts = model_data_hist_ts[i,:,:]
    eoc_field_ts = model_data_eoc_ts[i,:,:]

    hist_field_ua = model_data_hist_ua[i,:,:]
    eoc_field_ua = model_data_eoc_ua[i,:,:]

    hist_field_tel = model_data_hist_tel[i,:,:]

    bias_values_pr[i] = numpy.sqrt( numpy.mean((hist_field_pr - obs_field_pr)**2.) )
    convergence_values_pr[i] = numpy.sqrt( numpy.mean((eoc_field_pr - model_field_mmem_pr)**2.) )
    
    bias_values_ts[i] = numpy.sqrt( numpy.mean((hist_field_ts - obs_field_ts)**2.) )
    convergence_values_ts[i] = numpy.sqrt( numpy.mean((eoc_field_ts - model_field_mmem_ts)**2.) )
    
    bias_values_ua[i] = numpy.sqrt( numpy.mean((hist_field_ua - obs_field_ua)**2.) )
    convergence_values_ua[i] = numpy.sqrt( numpy.mean((eoc_field_ua - model_field_mmem_ua)**2.) )
    
    bias_values_tel[i] = numpy.sqrt( numpy.mean((hist_field_tel - obs_field_tel)**2.) )
    correlation_vals_tel[i] = scipy.stats.pearsonr(hist_field_tel.flatten(), obs_field_tel.flatten())[0]

mmem_bias_pr = numpy.sqrt( numpy.mean( (numpy.mean(model_data_hist_pr, axis=0) - obs_field_pr)**2. ))
mmem_bias_ts = numpy.sqrt( numpy.mean( (numpy.mean(model_data_hist_ts, axis=0) - obs_field_ts)**2. ))
mmem_bias_ua = numpy.sqrt( numpy.mean( (numpy.mean(model_data_hist_ua, axis=0) - obs_field_ua)**2. ))
mmem_bias_tel = numpy.sqrt( numpy.mean( (numpy.mean(model_data_hist_tel, axis=0) - obs_field_tel)**2. ))

bias_values_pr_LENS = numpy.zeros((40))
bias_values_ts_LENS = numpy.zeros((40))
bias_values_ua_LENS = numpy.zeros((40))

for i in range(40):
    hist_field_pr = model_data_hist_pr_LENS[i,:,:]
    hist_field_ts = model_data_hist_ts_LENS[i,:,:]
    hist_field_ua = model_data_hist_ua_LENS[i,:,:]
    
    bias_values_pr_LENS[i] = numpy.sqrt( numpy.mean( (hist_field_pr - obs_field_pr)**2.) )
    bias_values_ts_LENS[i] = numpy.sqrt( numpy.mean( (hist_field_ts - obs_field_ts)**2.) )
    bias_values_ua_LENS[i] = numpy.sqrt( numpy.mean( (hist_field_ua - obs_field_ua)**2.) )


# In[13]:

# create dictionaries to be used below
dict_pr = {
'bias_values_mods':bias_values_pr,
'convergence_values_mods':convergence_values_pr,
'bias_values_LENS':bias_values_pr_LENS,
'mmem_bias':mmem_bias_pr,
'nlat':pr_regional_nlat,
'nlon':pr_regional_nlon,
'lats':pr_regional_lat_vals,
'lons':pr_regional_lon_vals,
'fields_hist_mods':model_data_hist_pr,
'fields_hist_mods_LENS':model_data_hist_pr_LENS,
'fields_eoc_mods':model_data_eoc_pr,
'obs_field':obs_field_pr,
'LENS':True
}

dict_ts = {
'bias_values_mods':bias_values_ts,
'convergence_values_mods':convergence_values_ts,
'bias_values_LENS':bias_values_ts_LENS,
'mmem_bias':mmem_bias_ts,
'nlat':ts_regional_nlat,
'nlon':ts_regional_nlon,
'lats':ts_regional_lat_vals,
'lons':ts_regional_lon_vals,
'fields_hist_mods':model_data_hist_ts,
'fields_hist_mods_LENS':model_data_hist_ts_LENS,
'fields_eoc_mods':model_data_eoc_ts,
'obs_field':obs_field_ts,
'LENS':True
}

dict_ua = {
'bias_values_mods':bias_values_ua,
'convergence_values_mods':convergence_values_ua,
'bias_values_LENS':bias_values_ua_LENS,
'mmem_bias':mmem_bias_ua,
'nlat':ua_regional_nlat,
'nlon':ua_regional_nlon,
'lats':ua_regional_lat_vals,
'lons':ua_regional_lon_vals,
'fields_hist_mods':model_data_hist_ua,
'fields_hist_mods_LENS':model_data_hist_ua_LENS,
'fields_eoc_mods':model_data_eoc_ua,
'obs_field':obs_field_ua,
'LENS':True
}

dict_tel = {
'bias_values_mods':bias_values_tel,
'mmem_bias':mmem_bias_tel,
'LENS':False,
'nlat':tel_regional_nlat,
'nlon':tel_regional_nlon,
'fields_hist_mods':model_data_hist_tel,
'obs_field':obs_field_tel
}


# In[14]:

pareto_set_collect_2d_list = []
pareto_set_collect_3d_list = []


# # RASTERIZED VERSION

# # Calculate all biases for CMIP5

# In[ ]:

DATESTRING = datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

dict_x=dict_pr
dict_y=dict_ts
dict_z=dict_ua

N_pareto_loops=5

##################################################
##################################################
##################################################
# do N choose k subensembles
# for each, calculate ensemble mean

k=7

model_numbers = numpy.arange(nmods, dtype=numpy.int)
model_combinations = list(itertools.combinations(model_numbers, k))
random.shuffle(model_combinations)
model_combinations = numpy.array(model_combinations, dtype=numpy.int)

#N_ens = int(2e6) #model_combinations.shape[0]
N_ens = model_combinations.shape[0]
model_combinations = model_combinations[0:N_ens,:]

##################################################
#save these things

save_dict = {}

save_dict['N_ens'] = N_ens
save_dict['model_combinations'] = model_combinations

save_dict['dict_x'] = dict_x
save_dict['dict_y'] = dict_y
save_dict['dict_z'] = dict_z

save_dir = '/Volumes/data1/baird/subensembles_analysis/sampling_k_values/data_files/'
save_filename = 'pareto_front_results_'+DATESTRING+'_NO_LENS_model_combos.npy'
numpy.save(save_dir + save_filename, save_dict)
print('saved model combinations')