"""
compare_obs_files.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 9 May 2024

This script reads in AirNow obs from various files to compare them to identify potential bugs and ensure consistency.
"""

import sys
import pathlib
import argparse
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc

def main():
    beg_date_str = '2020-08-01_00'
    end_date_str = '2020-08-01_03'
    stride_h = 1

    obs_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static','20200801_00-20211231_23')
    cams_ra_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4')
    cams_fc_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')
    cmaq_fc_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','cmaq_hourly','sites_static')
    cams_fc_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','sites_static')
    cmaq_fc_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','sites_static')

    fmt_yyyymmdd_hhmm = '%Y%m%d_%H%M'
    fmt_yyyymmdd_hh = '%Y%m%d_%H'
    fmt_yyyymmddhh  = '%Y%m%d%H'
    fmt_yyyymmdd = '%Y%m%d'
    fmt_yyyy = '%Y'
    fmt_mm = '%m'
    fmt_hh = '%H'
    fmt_date = '%Y-%m-%d'
    fmt_date_hh = '%Y-%m-%d_%H'
    fmt_date_file = fmt_yyyymmdd_hhmm
    fmt_date_plot = '%d %b %Y/%H%M'

    ## Create datetime array of valid times
    dt_beg_valid = dt.datetime.strptime(beg_date_str, fmt_date_hh)
    dt_end_valid = dt.datetime.strptime(end_date_str, fmt_date_hh)
    dt_all_valid = pd.date_range(start=dt_beg_valid, end=dt_end_valid, freq=str(stride_h)+'H')
    n_valid_times = len(dt_all_valid)

    ## Loop over valid times
    for tt in range(n_valid_times):
        this_dt = dt_all_valid[tt]
        this_date = this_dt.strftime(fmt_yyyymmdd)
        this_yr = this_dt.strftime(fmt_yyyy)
        this_mo = this_dt.strftime(fmt_mm)
        this_hr = this_dt.strftime(fmt_hh)
        this_date_file = this_dt.strftime(fmt_date_file)

        obs_dir  = obs_dir_parent.joinpath(this_yr, this_mo)
        use_file = obs_dir.joinpath('airnow_pm2.5_o3_'+this_date_file+'_use.nc')
        val_file = obs_dir.joinpath('airnow_pm2.5_o3_'+this_date_file+'_val.nc')

        print('Reading '+str(use_file))
        ds_use = xr.open_dataset(use_file)
        o3_sid_use = ds_use.o3_sid_use
        o3_lon_use = ds_use.o3_lon_use
        o3_lat_use = ds_use.o3_lat_use
        o3_con_use = ds_use.o3_con_use
        pm25_sid_use = ds_use.pm25_sid_use
        pm25_lon_use = ds_use.pm25_lon_use
        pm25_lat_use = ds_use.pm25_lat_use
        pm25_con_use = ds_use.pm25_con_use
        n_use_o3   = len(o3_sid_use)
        n_use_pm25 = len(pm25_sid_use)

        print('Reading '+str(val_file))
        ds_val = xr.open_dataset(val_file)
        o3_sid_val = ds_val.o3_sid_val
        o3_lon_val = ds_val.o3_lon_val
        o3_lat_val = ds_val.o3_lat_val
        o3_con_val = ds_val.o3_con_val
        pm25_sid_val = ds_val.pm25_sid_val
        pm25_lon_val = ds_val.pm25_lon_val
        pm25_lat_val = ds_val.pm25_lat_val
        pm25_con_val = ds_val.pm25_con_val
        n_val_o3   = len(o3_sid_val)
        n_val_pm25 = len(pm25_sid_val)

        ## CAMS RA Raw
        this_dir = cams_ra_raw_dir.joinpath(this_yr, this_mo)
        fname_in_use = this_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_static_'+this_date_file+'_use.nc')
        fname_in_val = this_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_static_'+this_date_file+'_val.nc')

        print('Reading '+str(fname_in_use))
        ds_use = xr.open_dataset(fname_in_use)
        cams_ra_raw_obs_o3_con_use = ds_use.obs_o3_con_use
        cams_ra_raw_obs_o3_lat_use = ds_use.obs_o3_lat_use
        cams_ra_raw_obs_o3_lon_use = ds_use.obs_o3_lon_use
        cams_ra_raw_obs_pm25_con_use = ds_use.obs_pm25_con_use
        cams_ra_raw_obs_pm25_lat_use = ds_use.obs_pm25_lat_use
        cams_ra_raw_obs_pm25_lon_use = ds_use.obs_pm25_lon_use

        ## TODO: Keep building out this script if necessary.
        ##       Scott thinks he found a bug in his scripts, which would make this script unnecessary.


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
