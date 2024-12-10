'''
interrogate_obs_files.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 22 Mar 2023

This script's purpose is to examine any given time series of obs for debugging purposes.
'''

import sys
import pathlib
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt

cycle_date_str = '20201125_0600'
sid = 320311005
n_leads = 73

data_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cmaq','raw','airnow','sites_static')

fmt_yyyymmdd = '%Y%m%d'
fmt_yyyymmdd_hhmm = '%Y%m%d_%H%M'
cycle_dt = pd.to_datetime(cycle_date_str, format=fmt_yyyymmdd_hhmm)
cycle_yr = cycle_dt.strftime('%Y')
cycle_mo = cycle_dt.strftime('%m')

data_dir = data_dir_parent.joinpath(cycle_yr,cycle_mo)
data_file_use = data_dir.joinpath('cmaq_airnow_pm2.5_o3_'+cycle_date_str+'_use.nc')
data_file_val = data_dir.joinpath('cmaq_airnow_pm2.5_o3_'+cycle_date_str+'_val.nc')

print('Reading '+str(data_file_use))
use_ds = xr.open_dataset(data_file_use)

obs_o3_sid_use = use_ds.obs_o3_sid_use.values.astype(int)
obs_o3_lat_use = use_ds.obs_o3_lat_use.values
obs_o3_lon_use = use_ds.obs_o3_lon_use.values
obs_o3_con_use = use_ds.obs_o3_con_use.values

print('Reading '+str(data_file_val))
val_ds = xr.open_dataset(data_file_val)

obs_o3_sid_val = val_ds.obs_o3_sid_val.values.astype(int)
obs_o3_lat_val = val_ds.obs_o3_lat_val.values
obs_o3_lon_val = val_ds.obs_o3_lon_val.values
obs_o3_con_val = val_ds.obs_o3_con_val.values

sid_ind_use = np.where(obs_o3_sid_use == sid)[0][0]
#sid_ind_val = np.where(obs_o3_sid_val == sid)[0]

print('O3 concentration observations for station ID '+str(sid)+' ('+str(obs_o3_lon_use[sid_ind_use])+' E, '+str(obs_o3_lat_use[sid_ind_use])+' N):')
for ll in range(n_leads):
	this_dt = cycle_dt+dt.timedelta(hours=ll)
	this_dt_str = this_dt.strftime(fmt_yyyymmdd_hhmm)
	print('valid date = '+this_dt_str+', lead = '+str(ll)+', O3 = '+str(obs_o3_con_use[sid_ind_use,ll]))
