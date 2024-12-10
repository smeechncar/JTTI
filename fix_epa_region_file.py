'''
fix_epa_region_file.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 2 Nov 2023

This script reads in AirNow EPA region files from Ju-Hye and AirNow use/val obs files I created.
The PM2.5 obs are jumbled in Ju-Hye's files somehow, so the EPA region numbers are remapped onto
the use/val obs from the files I created. The O3 obs are fine. New files are then written out.
'''

import sys
import pathlib
import argparse
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
from functions import gen_funcs

def main():
    epa_codes_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI')

    ## Read in the EPA codes files that Ju-Hye made
    ## Note that these data are float32, while my files are float64, which causes inequalities
    epa_codes_file_use = epa_codes_dir.joinpath('EPA_CODE_O3_PM25_USE.nc')
    epa_codes_file_val = epa_codes_dir.joinpath('EPA_CODE_O3_PM25_VAL.nc')

    print('Reading '+str(epa_codes_file_use))
    epa_use_ds = xr.open_dataset(epa_codes_file_use)
    jh_epa_code_o3_use   = epa_use_ds.epa_code_o3.values.astype(int)
    jh_obs_o3_lat_use    = epa_use_ds.obs_o3_lat_use.values
    jh_obs_o3_lon_use    = epa_use_ds.obs_o3_lon_use.values
    jh_epa_code_pm25_use = epa_use_ds.epa_code_pm25.values.astype(int)
    jh_obs_pm25_lat_use  = epa_use_ds.obs_pm25_lat_use.values
    jh_obs_pm25_lon_use  = epa_use_ds.obs_pm25_lon_use.values
    n_jh_obs_o3_use   = len(jh_epa_code_o3_use)
    n_jh_obs_pm25_use = len(jh_epa_code_pm25_use)

    print('Reading '+str(epa_codes_file_val))
    epa_val_ds = xr.open_dataset(epa_codes_file_val)
    jh_epa_code_o3_val   = epa_val_ds.epa_code_o3.values.astype(int)
    jh_obs_o3_lat_val    = epa_val_ds.obs_o3_lat_val.values
    jh_obs_o3_lon_val    = epa_val_ds.obs_o3_lon_val.values
    jh_epa_code_pm25_val = epa_val_ds.epa_code_pm25.values.astype(int)
    jh_obs_pm25_lat_val  = epa_val_ds.obs_pm25_lat_val.values
    jh_obs_pm25_lon_val  = epa_val_ds.obs_pm25_lon_val.values
    n_jh_obs_o3_val   = len(jh_epa_code_o3_val)
    n_jh_obs_pm25_val = len(jh_epa_code_pm25_val)
#    print(n_jh_obs_o3_use)
#    print(n_jh_obs_o3_val)
#    print(n_jh_obs_pm25_use)
#    print(n_jh_obs_pm25_val)

    ## Combine the arrays
    jh_epa_code_o3_all   = np.append(jh_epa_code_o3_use, jh_epa_code_o3_val)
    jh_obs_o3_lat_all    = np.append(jh_obs_o3_lat_use, jh_obs_o3_lat_val)
    jh_obs_o3_lon_all    = np.append(jh_obs_o3_lon_use, jh_obs_o3_lon_val)
    jh_epa_code_pm25_all = np.append(jh_epa_code_pm25_use, jh_epa_code_pm25_val)
    jh_obs_pm25_lat_all  = np.append(jh_obs_pm25_lat_use, jh_obs_pm25_lat_val)
    jh_obs_pm25_lon_all  = np.append(jh_obs_pm25_lon_use, jh_obs_pm25_lon_val)
    n_jh_obs_o3_all   = len(jh_epa_code_o3_all)
    n_jh_obs_pm25_all = len(jh_epa_code_pm25_all)

    ## Read in a sample use/val obs file to get station metadata
    obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static','20200801_00-20211231_23','2020','08')
    obs_file_use = obs_dir.joinpath('airnow_pm2.5_o3_20200801_0000_use.nc')
    obs_file_val = obs_dir.joinpath('airnow_pm2.5_o3_20200801_0000_val.nc')

    print('Reading '+str(obs_file_use))
    obs_use_ds = xr.open_dataset(obs_file_use)
    obs_o3_lat_use = obs_use_ds.o3_lat_use
    obs_o3_lon_use = obs_use_ds.o3_lon_use
    obs_o3_sid_use = obs_use_ds.o3_sid_use
    obs_pm25_lat_use = obs_use_ds.pm25_lat_use
    obs_pm25_lon_use = obs_use_ds.pm25_lon_use
    obs_pm25_sid_use = obs_use_ds.pm25_sid_use

    print('Reading '+str(obs_file_val))
    obs_val_ds = xr.open_dataset(obs_file_val)
    obs_o3_lat_val = obs_val_ds.o3_lat_val
    obs_o3_lon_val = obs_val_ds.o3_lon_val
    obs_o3_sid_val = obs_val_ds.o3_sid_val
    obs_pm25_lat_val = obs_val_ds.pm25_lat_val
    obs_pm25_lon_val = obs_val_ds.pm25_lon_val
    obs_pm25_sid_val = obs_val_ds.pm25_sid_val

    ## Make float32 versions of these variables for comparison with variables from Ju-Hye's files
    obs_o3_lat_use_f32 = obs_o3_lat_use.values.astype(np.float32)
    obs_o3_lat_val_f32 = obs_o3_lat_val.values.astype(np.float32)
    obs_o3_lon_use_f32 = obs_o3_lon_use.values.astype(np.float32)
    obs_o3_lon_val_f32 = obs_o3_lon_val.values.astype(np.float32)
    obs_pm25_lat_use_f32 = obs_pm25_lat_use.values.astype(np.float32)
    obs_pm25_lat_val_f32 = obs_pm25_lat_val.values.astype(np.float32)
    obs_pm25_lon_use_f32 = obs_pm25_lon_use.values.astype(np.float32)
    obs_pm25_lon_val_f32 = obs_pm25_lon_val.values.astype(np.float32)

    n_obs_o3_use = len(obs_o3_lat_use)
    n_obs_o3_val = len(obs_o3_lat_val)
    n_obs_pm25_use = len(obs_pm25_lat_use)
    n_obs_pm25_val = len(obs_pm25_lon_val)
#    print(n_obs_o3_use)
#    print(n_obs_o3_val)
#    print(n_obs_pm25_use)
#    print(n_obs_pm25_val)

    ## Initialize arrays for the reordered EPA region variables
    obs_o3_reg_use = np.full(n_obs_o3_use, 0)
    obs_o3_reg_val = np.full(n_obs_o3_val, 0)
    obs_pm25_reg_use = np.full(n_obs_pm25_use, 0)
    obs_pm25_reg_val = np.full(n_obs_pm25_val, 0)


    '''
    if np.allclose(obs_o3_lat_use, jh_obs_o3_lat_use):
        print('obs_o3_lat_use and jh_obs_o3_lat_use match!')
    else:
        print('mismatch between obs_o3_lat_use and jh_obs_o3_lat_use')
        for nn in range(n_obs_o3_use):
            if np.allclose(obs_o3_lat_use[nn], jh_obs_o3_lat_use[nn]):
                print('nn = '+str(nn)+', equal: lat = '+str(obs_o3_lat_use[nn])+', jh_lat = '+str(jh_obs_o3_lat_use[nn]))
            else:
                print('nn = '+str(nn)+', NOT equal: lat = '+str(obs_o3_lat_use[nn])+', jh_lat = '+str(jh_obs_o3_lat_use[nn]))

    if np.allclose(obs_o3_lat_val, jh_obs_o3_lat_val):
        print('obs_o3_lat_val and jh_obs_o3_lat_val match!')
    else:
        print('mismatch between obs_o3_lat_val and jh_obs_o3_lat_val')
        for nn in range(n_obs_o3_val):
            if np.allclose(obs_o3_lat_val[nn], jh_obs_o3_lat_val[nn]):
                print('nn = '+str(nn)+', equal: lat = '+str(obs_o3_lat_val[nn])+', jh_lat = '+str(jh_obs_o3_lat_val[nn]))
            else:
                print('nn = '+str(nn)+', NOT equal: lat = '+str(obs_o3_lat_val[nn])+', jh_lat = '+str(jh_obs_o3_lat_val[nn]))

    if np.array_equal(obs_pm25_lat_use, jh_obs_pm25_lat_use):
        print('obs_pm25_lat_use and jh_obs_pm25_lat_use match!')
    else:
        print('mismatch between obs_pm25_lat_use and jh_obs_pm25_lat_use')
        for nn in range(n_obs_pm25_use):
            if np.allclose(obs_pm25_lat_use[nn], jh_obs_pm25_lat_use[nn]):
                print('nn = '+str(nn)+', equal: lat = '+str(obs_pm25_lat_use[nn])+', jh_lat = '+str(jh_obs_pm25_lat_use[nn]))
            else:
                print('nn = '+str(nn)+', NOT equal: lat = '+str(obs_pm25_lat_use[nn])+', jh_lat = '+str(jh_obs_pm25_lat_use[nn]))

    if np.array_equal(obs_pm25_lat_val, jh_obs_pm25_lat_val):
        print('obs_pm25_lat_val and jh_obs_pm25_lat_val match!')
    else:
        print('mismatch between obs_pm25_lat_val and jh_obs_pm25_lat_val')
        for nn in range(n_obs_pm25_val):
            if np.allclose(obs_pm25_lat_val[nn], jh_obs_pm25_lat_val[nn]):
                print('nn = '+str(nn)+', equal: lat = '+str(obs_pm25_lat_val[nn])+', jh_lat = '+str(jh_obs_pm25_lat_val[nn]))
            else:
                print('nn = '+str(nn)+', NOT equal: lat = '+str(obs_pm25_lat_val[nn])+', jh_lat = '+str(jh_obs_pm25_lat_val[nn]))
    '''

    ## Get the index value of each use/val obs from the combined set of Ju-Hye's obs/region data
    ## I know from testing that the O3 arrays are ordered just fine, but there's no reason not to double-check
    for nn in range(n_obs_o3_use):
        ind = np.where((obs_o3_lat_use_f32[nn] == jh_obs_o3_lat_all) &
                (obs_o3_lon_use_f32[nn] == jh_obs_o3_lon_all))[0][0]
        obs_o3_reg_use[nn] = jh_epa_code_o3_all[ind]

    for nn in range(n_obs_o3_val):
        ind = np.where((obs_o3_lat_val_f32[nn] == jh_obs_o3_lat_all) &
                (obs_o3_lon_val_f32[nn] == jh_obs_o3_lon_all))[0][0]
        obs_o3_reg_val[nn] = jh_epa_code_o3_all[ind]

    for nn in range(n_obs_pm25_use):
        ind = np.where((obs_pm25_lat_use_f32[nn] == jh_obs_pm25_lat_all) &
                (obs_pm25_lon_use_f32[nn] == jh_obs_pm25_lon_all))[0][0]
        obs_pm25_reg_use[nn] = jh_epa_code_pm25_all[ind]

    for nn in range(n_obs_pm25_val):
        ind = np.where((obs_pm25_lat_val_f32[nn] == jh_obs_pm25_lat_all) &
                (obs_pm25_lon_val_f32[nn] == jh_obs_pm25_lon_all))[0][0]
        obs_pm25_reg_val[nn] = jh_epa_code_pm25_all[ind]

    ## Create xarray DataArrays, then write out the new netcdf files
    n_use_o3 = n_obs_o3_use
    n_val_o3 = n_obs_o3_val
    n_use_pm25 = n_obs_pm25_use
    n_val_pm25 = n_obs_pm25_val

    n_obs_use_o3 = n_obs_o3_use
    n_obs_val_o3 = n_obs_o3_val
    n_obs_use_pm25 = n_obs_pm25_use
    n_obs_val_pm25 = n_obs_pm25_val

    obs_o3_reg_use = xr.DataArray(obs_o3_reg_use,
                        coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
                        attrs={'description':'EPA region (99=outside U.S.) for the use (training) AirNow O3 sites'})
    obs_o3_reg_val = xr.DataArray(obs_o3_reg_val,
                        coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
                        attrs={'description':'EPA region (99=outside U.S.) for the val (validation) AirNow O3 sites'})
    obs_pm25_reg_use = xr.DataArray(obs_pm25_reg_use,
                        coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
                        attrs={'description':'EPA region (99=outside U.S.) for the use (training) AirNow PM2.5 sites'})
    obs_pm25_reg_val = xr.DataArray(obs_pm25_reg_val,
                        coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
                        attrs={'description':'EPA region (99=outside U.S.) for the val (validation) AirNow PM2.5 sites'})

    ## Create new xarray datasets
    ds_use = xr.Dataset(
                data_vars={
                    'obs_o3_lat_use':obs_o3_lat_use, 'obs_o3_lon_use':obs_o3_lon_use,
                    'obs_o3_reg_use':obs_o3_reg_use, 'obs_o3_sid_use':obs_o3_sid_use,
                    'obs_pm25_lat_use':obs_pm25_lat_use, 'obs_pm25_lon_use':obs_pm25_lon_use,
                    'obs_pm25_reg_use':obs_pm25_reg_use, 'obs_pm25_sid_use':obs_pm25_sid_use,  },
                coords={'n_obs_use_o3':n_use_o3, 'n_obs_use_pm25':n_use_pm25},
                attrs={'description':'Metadata for AirNow O3 & PM2.5 stations in the use (training) set'}, )

    ds_val = xr.Dataset(
                data_vars={
                    'obs_o3_lat_val':obs_o3_lat_val, 'obs_o3_lon_val':obs_o3_lon_val,
                    'obs_o3_reg_val':obs_o3_reg_val, 'obs_o3_sid_val':obs_o3_sid_val,
                    'obs_pm25_lat_val':obs_pm25_lat_val, 'obs_pm25_lon_val':obs_pm25_lon_val,
                    'obs_pm25_reg_val':obs_pm25_reg_val, 'obs_pm25_sid_val':obs_pm25_sid_val,  },
                coords={'n_obs_val_o3':n_val_o3, 'n_obs_val_pm25':n_val_pm25},
                attrs={'description':'Metadata for AirNow O3 & PM2.5 stations in the val (validation) set'}, )

    ## Set the filenames
    fname_use = epa_codes_dir.joinpath('epa_code_o3_pm25_use_jared.nc')
    fname_val = epa_codes_dir.joinpath('epa_code_o3_pm25_val_jared.nc')

    ## Write the datasets to NetCDF
    print('Writing '+str(fname_use))
    ds_use.to_netcdf(fname_use)
    fname_use.chmod(0o644)
    print('Writing '+str(fname_val))
    ds_val.to_netcdf(fname_val)
    fname_val.chmod(0o644)


if __name__ == '__main__':
    now_time_beg = dt.datetime.utcnow()
    main()
    now_time_end = dt.datetime.utcnow()
    run_time_tot = now_time_end - now_time_beg
    now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
    now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
    print('\nScript completed successfully.')
    print('   Beg time: '+now_time_beg_str)
    print('   End time: '+now_time_end_str)
    print('   Run time: '+str(run_time_tot)+'\n')
