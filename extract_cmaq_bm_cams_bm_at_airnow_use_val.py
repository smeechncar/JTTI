'''
extract_cmaq_bm_cams_bm_at_airnow_use_val.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 9 Mar 2023

This script first reads in AirNow obs NetCDF files that have already been separated into use & val
sets, for use in merging & validation, respectively. Second, it reads in a CMAQ forecast file to
get the grid coordinates of those stations and spatially interpolate the gridded values to the
AirNow sites via nearest-neighbor interpolation. Third, it reads in bias-corrected CMAQ PM2.5 & O3
files that are produced hourly to serve as gridded "analysis"/"truth". Fourth, it reads in bias-
corrected CAMS PM2.5 & O3 hourly "analysis"/"truth" files that have been regridded to the CMAQ grid.
The output files with AirNow obs, BC CMAQ-interpolated values, and BC CAMS-interpolated values,
are in NetCDF.
'''

import os
import sys
import pathlib
import argparse
import datetime as dt
import pandas as pd
import numpy as np
import xarray as xr
import scipy as sp
import scipy.interpolate
import geocat.comp
import pyproj
from functions import gen_funcs

def parse_args():
    ## Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('date_beg', help='beginning date/time to be processed [YYYYMMDD_HH]')
    parser.add_argument('-l', '--date_end', default=None, help='ending date/time to be processed [YYYYMMDD_HH]')

    args = parser.parse_args()
    date_beg = args.date_beg
    date_end = args.date_end

    if len(date_beg) != 11:
        print('ERROR! Incorrect length for positional argument date_beg. Exiting!')
        parser.print_help()
        sys.exit()
    elif date_beg[8] != '_':
        print('ERROR! Incorrect format for positional argument date_beg. Exiting!')
        parser.print_help()
        sys.exit()

    if date_end != None:
        if len(date_end) != 11:
            print('ERROR! Incorrect length for optional argument date_end. Exiting!')
            parser.print_help()
            sys.exit()
        elif date_end[8] != '_':
            print('ERROR! Incorrect format for optional argument date_end. Exiting!')
            parser.print_help()
            sys.exit()
    else:
        date_end = date_beg

    return date_beg, date_end


def main(date_beg, date_end):
    sites_vary = False
    sites_static = True

    ## Set directories
    if sites_vary:
        obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_vary')
    if sites_static:
        obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static',date_beg+'-'+date_end)

    ## PM2.5 and O3
    cmaq_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','BM3','static','cmaq_include0_nn')
    ## PM2.5
    cams_bm1_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','BM1','static','cams_include0')
    ## O3
    cams_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','BM3','static','cams_include0')

    fmt_yyyymmdd_hh = '%Y%m%d_%H'
    fmt_yyyymmddhh  = '%Y%m%d%H'
    fmt_yyyymmdd = '%Y%m%d'
    fmt_yyyy = '%Y'
    fmt_mm = '%m'
    fmt_hh = '%H'
    fmt_date = '%Y-%m-%d'
    fmt_date_hh = '%Y-%m-%d_%H'

    ## Create datetime array
    dt_beg = dt.datetime.strptime(date_beg, fmt_yyyymmdd_hh)
    dt_end = dt.datetime.strptime(date_end, fmt_yyyymmdd_hh)
    dt_array = pd.date_range(start=dt_beg, end=dt_end, freq='1H')
    n_times = len(dt_array)

    get_cmaq_lat_lon = True
    get_cmaq_proj = True

    ## Loop over times
    for tt in range(n_times):
        this_dt = dt_array[tt]
        this_yyyymmdd_hh = this_dt.strftime(fmt_yyyymmdd_hh)
        this_yyyymmddhh  = this_dt.strftime(fmt_yyyymmddhh)
        this_yyyymmdd = this_dt.strftime(fmt_yyyymmdd)
        this_yr = this_dt.strftime(fmt_yyyy)
        this_mo = this_dt.strftime(fmt_mm)
        this_hr = this_dt.strftime(fmt_hh)
        this_date_hh = this_dt.strftime(fmt_date_hh)

        ## First, read in the use & val AirNow PM2.5 & O3 observations
        obs_file_use = obs_dir.joinpath(this_yr,this_mo,'airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_use.nc')
        obs_file_val = obs_dir.joinpath(this_yr,this_mo,'airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_val.nc')

        print('Reading '+str(obs_file_use))
        obs_ds_use = xr.open_dataset(obs_file_use)

        obs_pm25_con_use = obs_ds_use.pm25_con_use
        obs_pm25_sid_use = obs_ds_use.pm25_sid_use
        obs_pm25_lon_use = obs_ds_use.pm25_lon_use
        obs_pm25_lat_use = obs_ds_use.pm25_lat_use
        obs_pm25_cmaq_i_use = obs_ds_use.pm25_cmaq_i_use
        obs_pm25_cmaq_j_use = obs_ds_use.pm25_cmaq_j_use

        obs_o3_con_use = obs_ds_use.o3_con_use
        obs_o3_sid_use = obs_ds_use.o3_sid_use
        obs_o3_lon_use = obs_ds_use.o3_lon_use
        obs_o3_lat_use = obs_ds_use.o3_lat_use
        obs_o3_cmaq_i_use = obs_ds_use.o3_cmaq_i_use
        obs_o3_cmaq_j_use = obs_ds_use.o3_cmaq_j_use

        print('Reading '+str(obs_file_val))
        obs_ds_val = xr.open_dataset(obs_file_val)

        obs_pm25_con_val = obs_ds_val.pm25_con_val
        obs_pm25_sid_val = obs_ds_val.pm25_sid_val
        obs_pm25_lon_val = obs_ds_val.pm25_lon_val
        obs_pm25_lat_val = obs_ds_val.pm25_lat_val
        obs_pm25_cmaq_i_val = obs_ds_val.pm25_cmaq_i_val
        obs_pm25_cmaq_j_val = obs_ds_val.pm25_cmaq_j_val

        obs_o3_con_val = obs_ds_val.o3_con_val
        obs_o3_sid_val = obs_ds_val.o3_sid_val
        obs_o3_lon_val = obs_ds_val.o3_lon_val
        obs_o3_lat_val = obs_ds_val.o3_lat_val
        obs_o3_cmaq_i_val = obs_ds_val.o3_cmaq_i_val
        obs_o3_cmaq_j_val = obs_ds_val.o3_cmaq_j_val

        ## Second, open a CMAQ sample/coordinate file to get grid info. Only need lat/lon, not map projection info.
        cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
        print('Reading CMAQ coordinate data from '+str(cmaq_fname))
        cmaq_ds = xr.open_dataset(cmaq_fname)
        cmaq_lon = cmaq_ds.LON[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})
        cmaq_lat = cmaq_ds.LAT[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})

        ## Get and store the nearest-neighbor grid indices for the AirNow use & val stations for PM2.5 & O3
        ## Only need to do this once
        if tt == 0:
            print('   Finding grid inds for o3_use')
            j_o3_use, i_o3_use = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_o3_lat_use, obs_o3_lon_use)
            print('   Finding grid inds for o3_val')
            j_o3_val, i_o3_val = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_o3_lat_val, obs_o3_lon_val)
            print('   Finding grid inds for pm25_use')
            j_pm25_use, i_pm25_use = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_pm25_lat_use, obs_pm25_lon_use)
            print('   Finding grid inds for pm25_val')
            j_pm25_val, i_pm25_val = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_pm25_lat_val, obs_pm25_lon_val)

        ## Third, read in CMAQ files (Ju-Hye has all CMAQ hourly files in one directory)
#       cmaq_fname_in_this   = cmaq_dir.joinpath('cmaq_airnow_pm2.5_o3_'+this_yyyymmddhh+'.nc')
        cmaq_fname_in_this_o3   = cmaq_bm3_dir.joinpath('cmaq_o3_'+this_yyyymmdd_hh+'00_BM3_static_inc0_nn.nc')
        cmaq_fname_in_this_pm25 = cmaq_bm3_dir.joinpath('cmaq_pm25_'+this_yyyymmdd_hh+'00_BM3_static_inc0_nn.nc')

        if cmaq_fname_in_this_o3.is_file():
            ds_cmaq_o3 = xr.open_dataset(cmaq_fname_in_this_o3)
            cmaq_o3_this = ds_cmaq_o3.cmaq_o3_m
            ## At the first time all values will be 0.0. Make it nan so it doesn't accidentally get used downstream.
            if np.all(cmaq_o3_this == 0.0):
                cmaq_o3_this[:,:] = np.nan
        else:
            print('WARNING: File '+str(cmaq_fname_in_this_o3)+' does not exist!')

        if cmaq_fname_in_this_pm25.is_file():
            ds_cmaq_pm25 = xr.open_dataset(cmaq_fname_in_this_pm25)
            cmaq_pm25_this = ds_cmaq_pm25.cmaq_pm25_m
            if np.all(cmaq_pm25_this == 0.0):
                cmaq_pm25_this[:,:] = np.nan
        else:
            print('WARNING: File '+str(cmaq_fname_in_this_pm25)+' does not exist!')

        ## Fourth, read in CAMS files (Ju-Hye has all CAMS hourly files in one directory as well)
        cams_fname_in_this_o3   = cams_bm3_dir.joinpath('cams_o3_regrid_cmaq_'+this_yyyymmdd_hh+'00_BM3_static_include0.nc')
#       cams_fname_in_this_pm25 = cams_bm1_dir.joinpath('cams_pm25_regrid_cmaq_'+this_yyyymmdd_hh+'00_BM1_static_include0.nc')
        cams_fname_in_this_pm25 = cams_bm3_dir.joinpath('cams_pm25_regrid_cmaq_'+this_yyyymmdd_hh+'00_BM3_static_noQC.nc')

        if cams_fname_in_this_o3.is_file():
            ds_cams_o3 = xr.open_dataset(cams_fname_in_this_o3)
            cams_o3_this = ds_cams_o3.cams_o3_m
            if np.all(cams_o3_this == 0.0):
                cams_o3_this[:,:] = np.nan
        else:
            print('WARNING: File '+str(cams_fname_in_this_o3)+' does not exist!')

        if cams_fname_in_this_pm25.is_file():
            ds_cams_pm25 = xr.open_dataset(cams_fname_in_this_pm25)
            cams_pm25_this = ds_cams_pm25.cams_pm25_m
            if np.all(cams_pm25_this == 0.0):
                cams_pm25_this[:,:] = np.nan
        else:
            print('WARNING: File '+str(cams_fname_in_this_pm25)+' does not exist!')

        ## Create xarray data arrays using this data
        n_use_o3 = len(obs_o3_lat_use)
        n_val_o3 = len(obs_o3_lat_val)
        n_use_pm25 = len(obs_pm25_lat_use)
        n_val_pm25 = len(obs_pm25_lat_val)

        ## Using RectBivariateSpline for interpolation doesn't work on curvilinear grids like CAMS, WRF, etc.
        ## Instead, do nearest-neighbor interpolation to AirNow stations using grid indices found above
        print('Interpolating nearest-neighbor CMAQ values to AirNow stations')
        cmaq_o3_ppb_airnow_use = np.full(n_use_o3, np.nan)
        cmaq_o3_ppb_airnow_val = np.full(n_val_o3, np.nan)
        cmaq_pm25_ug_airnow_use = np.full(n_use_pm25, np.nan)
        cmaq_pm25_ug_airnow_val = np.full(n_val_pm25, np.nan)

        for nn in range(n_use_o3):
            cmaq_o3_ppb_airnow_use[nn] = cmaq_o3_this[j_o3_use[nn], i_o3_use[nn]]
        for nn in range(n_val_o3):
            cmaq_o3_ppb_airnow_val[nn] = cmaq_o3_this[j_o3_val[nn], i_o3_val[nn]]
        for nn in range(n_use_pm25):
            cmaq_pm25_ug_airnow_use[nn] = cmaq_pm25_this[j_pm25_use[nn], i_pm25_use[nn]]
        for nn in range(n_val_pm25):
            cmaq_pm25_ug_airnow_val[nn] = cmaq_pm25_this[j_pm25_val[nn], i_pm25_val[nn]]

        print('Interpolating nearest-neighbor CAMS (regridded onto CMAQ grid) values to AirNow stations')
        cams_o3_ppb_airnow_use = np.full(n_use_o3, np.nan)
        cams_o3_ppb_airnow_val = np.full(n_val_o3, np.nan)
        cams_pm25_ug_airnow_use = np.full(n_use_pm25, np.nan)
        cams_pm25_ug_airnow_val = np.full(n_val_pm25, np.nan)

        for nn in range(n_use_o3):
            cams_o3_ppb_airnow_use[nn] = cams_o3_this[j_o3_use[nn], i_o3_use[nn]]
        for nn in range(n_val_o3):
            cams_o3_ppb_airnow_val[nn] = cams_o3_this[j_o3_val[nn], i_o3_val[nn]]
        for nn in range(n_use_pm25):
            cams_pm25_ug_airnow_use[nn] = cams_pm25_this[j_pm25_use[nn], i_pm25_use[nn]]
        for nn in range(n_val_pm25):
            cams_pm25_ug_airnow_val[nn] = cams_pm25_this[j_pm25_val[nn], i_pm25_val[nn]]

        cmaq_str = 'bias-corrected CMAQ gridded analysis time series (BM3 BC model for both PM2.5 and O3)'
        cmaq_o3_con_use = xr.DataArray(cmaq_o3_ppb_airnow_use,
                                    coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
                                    attrs={'description':cmaq_str+' O3 concentration interpolated to AirNow O3 stations used for merging', 'units':'ppbv'})
        cmaq_o3_con_val = xr.DataArray(cmaq_o3_ppb_airnow_val,
                                    coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
                                    attrs={'description':cmaq_str+' O3 concentration interpolated to AirNow O3 stations withheld for validation', 'units':'ppbv'})

        cmaq_pm25_con_use = xr.DataArray(cmaq_pm25_ug_airnow_use,
                                    coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
                                    attrs={'description':cmaq_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations used for merging', 'units':'ug/m3'})
        cmaq_pm25_con_val = xr.DataArray(cmaq_pm25_ug_airnow_val,
                                    coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
                                    attrs={'description':cmaq_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations withheld for validation', 'units':'ug/m3'})

        cams_str = 'bias-corrected CAMS gridded analysis time series (BM3 BC model for both PM2.5 and O3)'
        cams_o3_con_use = xr.DataArray(cams_o3_ppb_airnow_use,
                                    coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
                                    attrs={'description':cams_str+' O3 concentration interpolated to AirNow O3 stations used for merging', 'units':'ppbv'})
        cams_o3_con_val = xr.DataArray(cams_o3_ppb_airnow_val,
                                    coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
                                    attrs={'description':cams_str+' O3 concentration interpolated to AirNow O3 stations withheld for validation', 'units':'ppbv'})

        cams_pm25_con_use = xr.DataArray(cams_pm25_ug_airnow_use,
                                    coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
                                    attrs={'description':cams_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations used for merging', 'units':'ug/m3'})
        cams_pm25_con_val = xr.DataArray(cams_pm25_ug_airnow_val,
                                    coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
                                    attrs={'description':cams_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations withheld for validation', 'units':'ug/m3'})


        ## Create new xarray datasets
        ds_use = xr.Dataset(
                        data_vars={ 'cmaq_bm3_pm25_con_use':cmaq_pm25_con_use, 'cmaq_bm3_o3_con_use':cmaq_o3_con_use,
                                        'cams_bm3_pm25_con_use':cams_pm25_con_use, 'cams_bm3_o3_con_use':cams_o3_con_use,
                                        'obs_pm25_con_use':obs_pm25_con_use, 'obs_o3_con_use':obs_o3_con_use,
                                        'obs_pm25_sid_use':obs_pm25_sid_use, 'obs_o3_sid_use':obs_o3_sid_use,
                                        'obs_pm25_lon_use':obs_pm25_lon_use, 'obs_o3_lon_use':obs_o3_lon_use,
                                        'obs_pm25_lat_use':obs_pm25_lat_use, 'obs_o3_lat_use':obs_o3_lat_use,
                                        'obs_pm25_cmaq_i_use':obs_pm25_cmaq_i_use, 'obs_o3_cmaq_i_use':obs_o3_cmaq_i_use,
                                        'obs_pm25_cmaq_j_use':obs_pm25_cmaq_j_use, 'obs_o3_cmaq_j_use':obs_o3_cmaq_j_use,
                            },
                        coords={'n_obs_use_pm25':n_use_pm25, 'n_obs_use_o3':n_use_o3},
                        attrs={'description':'Subset of PM2.5 & O3 AirNow observations and interpolated '+cmaq_str+' and '+cams_str+' values retained for use in merging'},
                        )

        ds_val = xr.Dataset(
                        data_vars={ 'cmaq_bm3_pm25_con_val':cmaq_pm25_con_val, 'cmaq_bm3_o3_con_val':cmaq_o3_con_val,
                                        'cams_bm3_pm25_con_val':cams_pm25_con_val, 'cams_bm3_o3_con_val':cams_o3_con_val,
                                        'obs_pm25_con_val':obs_pm25_con_val, 'obs_o3_con_val':obs_o3_con_val,
                                        'obs_pm25_sid_val':obs_pm25_sid_val, 'obs_o3_sid_val':obs_o3_sid_val,
                                        'obs_pm25_lon_val':obs_pm25_lon_val, 'obs_o3_lon_val':obs_o3_lon_val,
                                        'obs_pm25_lat_val':obs_pm25_lat_val, 'obs_o3_lat_val':obs_o3_lat_val,
                                        'obs_pm25_cmaq_i_val':obs_pm25_cmaq_i_val, 'obs_o3_cmaq_i_val':obs_o3_cmaq_i_val,
                                        'obs_pm25_cmaq_j_val':obs_pm25_cmaq_j_val, 'obs_o3_cmaq_j_val':obs_o3_cmaq_j_val,
                            },
                        coords={'n_obs_val_pm25':n_val_pm25, 'n_obs_val_o3':n_val_o3},
                        attrs={'description':'Subset of PM2.5 & O3 AirNow observations and interpolated '+cmaq_str+' and '+cams_str+' values withheld for validation'},
                        )

        ## Set the output paths & filenames
        if sites_vary:
            out_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','sites_vary',this_yr,this_mo)
            out_dir.mkdir(parents=True, exist_ok=True)
            out_dir.chmod(0o755)
            fname_use = out_dir.joinpath('cmaq_bc_cams_bc_airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_use.nc')
            fname_val = out_dir.joinpath('cmaq_bc_cams_bc_airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_val.nc')
        if sites_static:
            out_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','sites_static',this_yr,this_mo)
            out_dir.mkdir(parents=True, exist_ok=True)
            out_dir.chmod(0o755)
            fname_use = out_dir.joinpath('cmaq_bc_cams_bc_airnow_pm2.5_o3_static_'+this_yyyymmdd_hh+'00_use.nc')
            fname_val = out_dir.joinpath('cmaq_bc_cams_bc_airnow_pm2.5_o3_static_'+this_yyyymmdd_hh+'00_val.nc')

        ## Write the datasets to NetCDF
        print('Writing '+str(fname_use))
        ds_use.to_netcdf(fname_use)
        fname_use.chmod(0o644)
        print('Writing '+str(fname_val))
        ds_val.to_netcdf(fname_val)
        fname_val.chmod(0o644)


if __name__ == '__main__':
    now_time_beg = dt.datetime.utcnow()
    date_beg, date_end, = parse_args()
    main(date_beg, date_end)
    now_time_end = dt.datetime.utcnow()
    run_time_tot = now_time_end - now_time_beg
    now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
    now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
    print('\nScript completed successfully.')
    print('   Beg time: '+now_time_beg_str)
    print('   End time: '+now_time_end_str)
    print('   Run time: '+str(run_time_tot)+'\n')
