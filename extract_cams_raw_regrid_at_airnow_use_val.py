'''
extract_cams_raw_regrid_at_airnow_use_val.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 30 Oct 2023

This script first reads in AirNow obs NetCDF files that have already been separated into use & val
sets, for use in merging & validation, respectively. Second, it reads in either a CAMS forecast or
reanalysis (EAC4) file that has already been regridded to the CMAQ grid to get grid coordinates.
Due to RectBivariateSpline not working on curvilinear grids, nearest neighbor interpolation is used
to obtain CAMS data at the AirNow sites. CAMS PM2.5 files are produced hourly, but because CAMS O3
files are only produced 3-hourly, temporal gridded lienar interpolation is also required.
Regridding with regrid_cams_cmaq.py takes care of this linear interpolation, so files are available
every 1 hour, which this script then reads in.
The output files with CAMS-interpolated values at AirNow sites are in NetCDF.
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
from functions import gen_funcs

def parse_args():
    ## Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('date_beg', help='beginning date/time to be processed [YYYYMMDD_HH]')
    parser.add_argument('-l', '--date_end', default=None, help='ending date/time to be processed [YYYYMMDD_HH]')
    parser.add_argument('-f', '--forecast', default=False, help='flag to process CAMS gridded forecasts', action='store_true')
    parser.add_argument('-r', '--reanalysis', default=False, help='flag to process CAMS gridded reanalysis', action='store_true')

    args = parser.parse_args()
    date_beg = args.date_beg
    date_end = args.date_end

    cams_fcst = False
    if args.forecast:
        cams_fcst = True

    cams_eac4 = False
    if args.reanalysis:
        cams_eac4 = True

    if cams_fcst and cams_eac4:
        print('ERROR: Both -f (cams_fcst) and -r (cams_eac4) are used/True. Please only use one of them.')
        sys.exit()
    elif not cams_fcst and not cams_eac4:
        print('ERROR: Neither -f (cams_fcst) nor -r (cams_eac4) are used/True. Please only use one of them.')
        sys.exit()

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

    return date_beg, date_end, cams_fcst, cams_eac4


def main(date_beg, date_end, cams_fcst, cams_eac4):
    sites_vary = False
    sites_static = True

    ## Set directories
    if sites_vary:
        obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_vary')
    if sites_static:
#        obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static',date_beg+'-'+date_end)
        obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static','20200801_00-20211231_23')
    cams_fcst_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')
    cams_eac4_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4')

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

    get_cams_lat_lon = True

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

        ## Second, read in CAMS files
        '''
        ## Start by searching for files at neighboring times if needed for temporal interpolation
        prev2_dt = this_dt - dt.timedelta(hours=2)
        prev1_dt = this_dt - dt.timedelta(hours=1)
        next1_dt = this_dt + dt.timedelta(hours=1)
        next2_dt = this_dt + dt.timedelta(hours=2)

        prev2_yr = prev2_dt.strftime('%Y')
        prev1_yr = prev1_dt.strftime('%Y')
        next1_yr = next1_dt.strftime('%Y')
        next2_yr = next2_dt.strftime('%Y')

        prev2_mo = prev2_dt.strftime('%m')
        prev1_mo = prev1_dt.strftime('%m')
        next1_mo = next1_dt.strftime('%m')
        next2_mo = next2_dt.strftime('%m')

        prev2_hr = prev2_dt.strftime('%H')
        prev1_hr = prev1_dt.strftime('%H')
        next1_hr = next1_dt.strftime('%H')
        next2_hr = next2_dt.strftime('%H')

        prev2_yyyymmdd = prev2_dt.strftime(fmt_yyyymmdd)
        prev1_yyyymmdd = prev1_dt.strftime(fmt_yyyymmdd)
        next1_yyyymmdd = next1_dt.strftime(fmt_yyyymmdd)
        next2_yyyymmdd = next2_dt.strftime(fmt_yyyymmdd)

#       prev2_date_str = prev2_dt.strftime(fmt_date)
#       prev1_date_str = prev1_dt.strftime(fmt_date)
#       next1_date_str = next1_dt.strftime(fmt_date)
#       next2_date_str = next2_dt.strftime(fmt_date)

        prev2_date_hh = prev2_dt.strftime(fmt_date_hh)
        prev1_date_hh = prev1_dt.strftime(fmt_date_hh)
        next1_date_hh = next1_dt.strftime(fmt_date_hh)
        next2_date_hh = next2_dt.strftime(fmt_date_hh)

        if cams_fcst:
            if int(this_hr) < 12:
                this_cycle = this_yyyymmdd+'_00'
            else:
                this_cycle = this_yyyymmdd+'_12'

            this_cams_dir = cams_fcst_dir.joinpath(this_yr,this_mo,this_cycle)
            this_cycle_dt = pd.to_datetime(this_cycle, format=fmt_yyyymmdd_hh)

            if int(prev2_hr) < 12:
                prev2_cycle = prev2_yyyymmdd+'_00'
            else:
                prev2_cycle = prev2_yyyymmdd+'_12'

            if int(prev1_hr) < 12:
                prev1_cycle = prev1_yyyymmdd+'_00'
            else:
                prev1_cycle = prev1_yyyymmdd+'_12'

            if int(next1_hr) < 12:
                next1_cycle = next1_yyyymmdd+'_00'
            else:
                next1_cycle = next1_yyyymmdd+'_12'

            if int(next2_hr) < 12:
                next2_cycle = next2_yyyymmdd+'_00'
            else:
                next2_cycle = next2_yyyymmdd+'_12'

            prev2_cams_dir = cams_fcst_dir.joinpath(prev2_yr,prev2_mo,prev2_cycle)
            prev1_cams_dir = cams_fcst_dir.joinpath(prev1_yr,prev1_mo,prev1_cycle)
            next1_cams_dir = cams_fcst_dir.joinpath(next1_yr,next1_mo,next1_cycle)
            next2_cams_dir = cams_fcst_dir.joinpath(next2_yr,next2_mo,next2_cycle)

        elif cams_eac4:
            this_cams_dir = cams_eac4_dir.joinpath(this_yr,this_mo)

            prev2_cams_dir = cams_eac4_dir.joinpath(prev2_yr,prev2_mo)
            prev1_cams_dir = cams_eac4_dir.joinpath(prev1_yr,prev1_mo)
            next1_cams_dir = cams_eac4_dir.joinpath(next1_yr,next1_mo)
            next2_cams_dir = cams_eac4_dir.joinpath(next2_yr,next2_mo)

        ## Get CAMS file names
        cams_o3_fname_in_this   = this_cams_dir.joinpath('o3sfc_'+this_date_hh+'00.nc')
        cams_pm25_fname_in_this = this_cams_dir.joinpath('pm2.5_'+this_date_hh+'00.nc')

        cams_o3_fname_in_prev2   = prev2_cams_dir.joinpath('o3sfc_'+prev2_date_hh+'00.nc')
        cams_o3_fname_in_prev1   = prev1_cams_dir.joinpath('o3sfc_'+prev1_date_hh+'00.nc')
        cams_o3_fname_in_next1   = next1_cams_dir.joinpath('o3sfc_'+next1_date_hh+'00.nc')
        cams_o3_fname_in_next2   = next2_cams_dir.joinpath('o3sfc_'+next2_date_hh+'00.nc')
        cams_pm25_fname_in_prev2 = prev2_cams_dir.joinpath('pm2.5_'+prev2_date_hh+'00.nc')
        cams_pm25_fname_in_prev1 = prev1_cams_dir.joinpath('pm2.5_'+prev1_date_hh+'00.nc')
        cams_pm25_fname_in_next1 = next1_cams_dir.joinpath('pm2.5_'+next1_date_hh+'00.nc')
        cams_pm25_fname_in_next2 = next2_cams_dir.joinpath('pm2.5_'+next2_date_hh+'00.nc')

        ## Do these files exist? (True/False)
        isfile_this_o3  = cams_o3_fname_in_this.is_file()
        isfile_prev2_o3 = cams_o3_fname_in_prev2.is_file()
        isfile_prev1_o3 = cams_o3_fname_in_prev1.is_file()
        isfile_next1_o3 = cams_o3_fname_in_next1.is_file()
        isfile_next2_o3 = cams_o3_fname_in_next2.is_file()
        isfile_this_pm25  = cams_pm25_fname_in_this.is_file()
        isfile_prev2_pm25 = cams_pm25_fname_in_prev2.is_file()
        isfile_prev1_pm25 = cams_pm25_fname_in_prev1.is_file()
        isfile_next1_pm25 = cams_pm25_fname_in_next1.is_file()
        isfile_next2_pm25 = cams_pm25_fname_in_next2.is_file()

        ## Get CAMS ozone data
        ## Simplest case: The file exists, no interpolation needed
        if isfile_this_o3:
            ds_o3 = xr.open_dataset(cams_o3_fname_in_this)
#           print(ds_o3)

            if get_cams_lat_lon:
                ## Need to reverse CAMS latitude array to be strictly increasing, rather than strictly decreasing
                cams_lat = ds_o3.latitude[::-1]
                ## Need to put CAMS longitude on range (-180, 180] to be consistent with CMAQ
                cams_lon = ds_o3.longitude
                cams_lon_vals = np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon)
                cams_lon = xr.DataArray(np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon),
                    coords={'longitude':cams_lon_vals}, attrs={'units':'degrees_east', 'long_name':'longitude'})
                get_cams_lat_lon = False

            ## Need to reverse the latitude dimension just as cams_lat was reversed
            cams_o3_this = ds_o3.go3[:,::-1,:]

        ## Case 2: This time bounded by files at h-2 and h+1
        elif isfile_prev2_o3 and isfile_next1_o3:
            ds_o3_prev2 = xr.open_dataset(cams_o3_fname_in_prev2)
            ds_o3_next1 = xr.open_dataset(cams_o3_fname_in_next1)
            if get_cams_lat_lon:
                cams_lat = ds_o3_prev2.latitude[::-1]
                cams_lon = ds_o3_prev2.longitude
                cams_lon_vals = np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon)
                cams_lon = xr.DataArray(np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon),
                    coords={'longitude':cams_lon_vals}, attrs={'units':'degrees_east', 'long_name':'longitude'})
                get_cams_lat_lon = False

            ## Reverse the latitude dimension
            cams_o3_prev2 = ds_o3_prev2.go3[:,::-1,:]
            cams_o3_next1 = ds_o3_next1.go3[:,::-1,:]

            ## Temporally interpolate the data
#           cams_o3_this = (1/3)*cams_o3_prev2 + (2/3)*cams_o3_next1
            cams_o3_this_vals = (1/3)*cams_o3_prev2.values + (2/3)*cams_o3_next1.values
            cams_o3_this = cams_o3_prev2.copy(deep=True, data=cams_o3_this_vals)

        ## Case 3: This time bounded by files at h-1 and h+2
        elif isfile_prev1_o3 and isfile_next2_o3:
            ds_o3_prev1 = xr.open_dataset(cams_o3_fname_in_prev1)
            ds_o3_next2 = xr.open_dataset(cams_o3_fname_in_next2)
            if get_cams_lat_lon:
                cams_lat = ds_o3_prev1.latitude[::-1]
                cams_lon = ds_o3_prev1.longitude
                cams_lon_vals = np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon)
                cams_lon = xr.DataArray(np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon),
                    coords={'longitude':cams_lon_vals}, attrs={'units':'degrees_east', 'long_name':'longitude'})
                get_cams_lat_lon = False

            ## Reverse the latitude dimension
            cams_o3_prev1 = ds_o3_prev1.go3[:,::-1,:]
            cams_o3_next2 = ds_o3_next2.go3[:,::-1,:]

            ## Temporally interpolate the data
            ## This will result in a length-0 time dimension, which causes problems when writing it out:
#           cams_o3_this = (2/3)*cams_o3_prev1 + (1/3)*cams_o3_next2
            ## Need to do this in order to preserve a length-1 time dimension:
            cams_o3_this_vals = (2/3)*cams_o3_prev1.values + (1/3)*cams_o3_next2.values
            cams_o3_this = cams_o3_prev1.copy(deep=True, data=cams_o3_this_vals)

        else:
            print('ERROR: This time does not have O3 files at surrounding times for temporal interp. Exiting!')
            sys.exit()

        ## Convert O3 mass mixing ratio (kg/kg) to volume mixing ratio (ppb)
        ## https://confluence.ecmwf.int/pages/viewpage.action?pageId=153391710
        cams_o3_this_ppb = cams_o3_this * 1e9 * (28.9644 / 47.9982)

        ## Get CAMS PM2.5 data
        ## Simplest case: The file exists, no interpolation needed
        if isfile_this_pm25:
            ds_pm25 = xr.open_dataset(cams_pm25_fname_in_this)
#           print(ds_pm25)

            if get_cams_lat_lon:
                ## Need to reverse CAMS latitude array to be strictly increasing, rather than strictly decreasing
                cams_lat = ds_pm25.latitude[::-1]
                ## Need to put CAMS longitude on range (-180, 180] to be consistent with CMAQ
                cams_lon = ds_pm25.longitude
                cams_lon_vals = np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon)
                cams_lon = xr.DataArray(np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon),
                    coords={'longitude':cams_lon_vals}, attrs={'units':'degrees_east', 'long_name':'longitude'})
                get_cams_lat_lon = False

            ## Need to reverse the latitude dimension just as cams_lat was reversed
            cams_pm25_this = ds_pm25.pm2p5[:,::-1,:]

        ## Case 2: This time bounded by files at h-2 and h+1
        elif isfile_prev2_pm25 and isfile_next1_pm25:
            ds_pm25_prev2 = xr.open_dataset(cams_pm25_fname_in_prev2)
            ds_pm25_next1 = xr.open_dataset(cams_pm25_fname_in_next1)
            if get_cams_lat_lon:
                cams_lat = ds_pm25_prev2.latitude[::-1]
                cams_lon = ds_pm25_prev2.longitude
                cams_lon_vals = np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon)
                cams_lon = xr.DataArray(np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon),
                    coords={'longitude':cams_lon_vals}, attrs={'units':'degrees_east', 'long_name':'longitude'})
                get_cams_lat_lon = False
            cams_pm25_prev2 = ds_pm25_prev2.pm2p5[:,::-1,:]
            cams_pm25_next1 = ds_pm25_next1.pm2p5[:,::-1,:]

            ## Temporally interpolate the data
#           cams_pm25_this = (1/3)*cams_pm25_prev2 + (2/3)*cams_pm25_next1
            cams_pm25_this_vals = (1/3)*cams_pm25_prev2.values + (2/3)*cams_pm25_next1.values
            cams_pm25_this = cams_pm25_prev2.copy(deep=True, data=cams_pm25_this_vals)

        ## Case 3: This time bounded by files at h-1 and h+2
        elif isfile_prev1_pm25 and isfile_next2_pm25:
            ds_pm25_prev1 = xr.open_dataset(cams_pm25_fname_in_prev1)
            ds_pm25_next2 = xr.open_dataset(cams_pm25_fname_in_next2)
            if get_cams_lat_lon:
                cams_lat = ds_pm25_prev1.latitude[::-1]
                cams_lon = ds_pm25_prev1.longitude
                cams_lon_vals = np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon)
                cams_lon = xr.DataArray(np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon),
                    coords={'longitude':cams_lon_vals}, attrs={'units':'degrees_east', 'long_name':'longitude'})
                get_cams_lat_lon = False
            cams_pm25_prev1 = ds_pm25_prev1.pm2p5[:,::-1,:]
            cams_pm25_next2 = ds_pm25_next2.pm2p5[:,::-1,:]

            ## Temporally interpolate the data
#           cams_pm25_this = (2/3)*cams_pm25_prev1 + (1/3)*cams_pm25_next2
            cams_pm25_this_vals = (1/3)*cams_pm25_prev1.values + (1/3)*cams_pm25_next2.values
            cams_pm25_this = cams_pm25_prev1.copy(deep=True, data=cams_pm25_this_vals)

        else:
            print('ERROR: This time does not have PM2.5 files at surrounding times for temporal interp. Exiting!')
            sys.exit()

        ## Convert PM2.5 concentration from kg/m3 to ug/m3
        cams_pm25_this_ug = cams_pm25_this * 1e9

        ## Now bilinearly interpolate the CAMS O3 & PM2.5 fields to the AirNow station locations
        ## CAMS is already on a rectangular lat-lon grid, so no transformation is required.
        ## There cannot be any nans in the input 2D array of values, so use array where nan was replaced with _FillValue.
        ## This will necessarily impact values near the boundary of real data, so be careful at the edges.
        ## To keep extracted values of the spline positive, set the degree of the spline to 1 (true bilinear interp).
        ## If the AirNow lat/lon arrays only have a nan (i.e., no valid data), the output is a nan array as well.
        ## 27 Dec 2022: With some spot checks compared to surrounding grid points both near and far from the middle,
        ##   I verified that RectBivariateSpline gives sensible results, as CAMS is a rectangular lat-lon grid.
        ##   Hence there is no need to use nearest-neighbor or some other interpolation scheme on this grid.
        cams_o3_ppb_spline = sp.interpolate.RectBivariateSpline(cams_lat,cams_lon,cams_o3_this_ppb[0,:,:],kx=1,ky=1)
        cams_o3_ppb_airnow_use = cams_o3_ppb_spline.ev(obs_o3_lat_use, obs_o3_lon_use)
        cams_o3_ppb_airnow_val = cams_o3_ppb_spline.ev(obs_o3_lat_val, obs_o3_lon_val)

        cams_pm25_ug_spline = sp.interpolate.RectBivariateSpline(cams_lat,cams_lon,cams_pm25_this_ug[0,:,:],kx=1,ky=1)
        cams_pm25_ug_airnow_use = cams_pm25_ug_spline.ev(obs_pm25_lat_use, obs_pm25_lon_use)
        cams_pm25_ug_airnow_val = cams_pm25_ug_spline.ev(obs_pm25_lat_val, obs_pm25_lon_val)
        '''

        if tt == 0:
            ## Open a CMAQ sample/coordinate file to get grid info. Only need lat/lon, not map projection info.
            cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
            print('Reading CMAQ coordinate data from '+str(cmaq_fname))
            cmaq_ds = xr.open_dataset(cmaq_fname)
            cmaq_lon = cmaq_ds.LON[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})
            cmaq_lat = cmaq_ds.LAT[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})
            cams_lon = cmaq_lon
            cams_lat = cmaq_lat

            ## Get and store the nearest-neighbor grid indices for the AirNow use & val stations for PM2.5 & O3
            ## Only need to do this once
            print('   Finding grid inds for o3_use')
            j_o3_use, i_o3_use = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_o3_lat_use, obs_o3_lon_use)
            print('   Finding grid inds for o3_val')
            j_o3_val, i_o3_val = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_o3_lat_val, obs_o3_lon_val)
            print('   Finding grid inds for pm25_use')
            j_pm25_use, i_pm25_use = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_pm25_lat_use, obs_pm25_lon_use)
            print('   Finding grid inds for pm25_val')
            j_pm25_val, i_pm25_val = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_pm25_lat_val, obs_pm25_lon_val)

        ## Get CAMS file name
        if cams_fcst:
            if int(this_hr) < 12:
                this_cycle = this_yyyymmdd+'_00'
            else:
                this_cycle = this_yyyymmdd+'_12'

            this_cams_dir = cams_fcst_dir.joinpath(this_yr,this_mo,this_cycle)
            this_cycle_dt = pd.to_datetime(this_cycle, format=fmt_yyyymmdd_hh)

        elif cams_eac4:
            this_cams_dir = cams_eac4_dir.joinpath(this_yr,this_mo)

        cams_fname_in_this = this_cams_dir.joinpath('cams_o3_pm25_regrid_cmaq_'+this_yyyymmdd_hh+'00.nc')

        if cams_fname_in_this.is_file():
            ds_cams = xr.open_dataset(cams_fname_in_this)
            cams_o3_this   = ds_cams.cams_o3[0,:,:]
            cams_pm25_this = ds_cams.cams_pm25[0,:,:]
            ## At the first time all values may be 0.0. Make it nan so it doesn't accidentally get used downstream.
            if np.all(cams_o3_this == 0.0):
                cams_o3_this[:,:] = np.nan
            if np.all(cams_pm25_this == 0.0):
                cams_pm25_this[:,:] = np.nan

        else:
            print('WARNING: File '+str(cams_fname_in_this)+' does not exist!')

        ## Create xarray data arrays using this data
        n_use_o3 = len(obs_o3_lat_use)
        n_val_o3 = len(obs_o3_lat_val)
        n_use_pm25 = len(obs_pm25_lat_use)
        n_val_pm25 = len(obs_pm25_lat_val)

        ## Using RectBivariateSpline for interpolation doesn't work on curvilinear grids like CAMS, WRF, etc.
        ## Instead, do nearest-neighbor interpolation to AirNow stations using grid indices found above
        print('Interpolating nearest-neighbor CAMS values (regridded to CMAQ grid) to AirNow stations')
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

        if cams_fcst:
            cams_str = 'CAMS gridded forecast'
        elif cams_eac4:
            cams_str = 'CAMS gridded reanalysis'

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
                        data_vars={ 'cams_pm25_con_use':cams_pm25_con_use, 'cams_o3_con_use':cams_o3_con_use,
                                'obs_pm25_con_use':obs_pm25_con_use, 'obs_o3_con_use':obs_o3_con_use,
                                'obs_pm25_sid_use':obs_pm25_sid_use, 'obs_o3_sid_use':obs_o3_sid_use,
                                'obs_pm25_lon_use':obs_pm25_lon_use, 'obs_o3_lon_use':obs_o3_lon_use,
                                'obs_pm25_lat_use':obs_pm25_lat_use, 'obs_o3_lat_use':obs_o3_lat_use,
                                'obs_pm25_cmaq_i_use':obs_pm25_cmaq_i_use, 'obs_o3_cmaq_i_use':obs_o3_cmaq_i_use,
                                'obs_pm25_cmaq_j_use':obs_pm25_cmaq_j_use, 'obs_o3_cmaq_j_use':obs_o3_cmaq_j_use,
                            },
                        coords={'n_obs_use_pm25':n_use_pm25, 'n_obs_use_o3':n_use_o3},
                        attrs={'description':'Subset of PM2.5 & O3 AirNow observations and interpolated '+cams_str+' values retained for use in merging'},
                        )

        ds_val = xr.Dataset(
                        data_vars={ 'cams_pm25_con_val':cams_pm25_con_val, 'cams_o3_con_val':cams_o3_con_val,
                                'obs_pm25_con_val':obs_pm25_con_val, 'obs_o3_con_val':obs_o3_con_val,
                                'obs_pm25_sid_val':obs_pm25_sid_val, 'obs_o3_sid_val':obs_o3_sid_val,
                                'obs_pm25_lon_val':obs_pm25_lon_val, 'obs_o3_lon_val':obs_o3_lon_val,
                                'obs_pm25_lat_val':obs_pm25_lat_val, 'obs_o3_lat_val':obs_o3_lat_val,
                                'obs_pm25_cmaq_i_val':obs_pm25_cmaq_i_val, 'obs_o3_cmaq_i_val':obs_o3_cmaq_i_val,
                                'obs_pm25_cmaq_j_val':obs_pm25_cmaq_j_val, 'obs_o3_cmaq_j_val':obs_o3_cmaq_j_val,
                            },
                        coords={'n_obs_val_pm25':n_val_pm25, 'n_obs_val_o3':n_val_o3},
                        attrs={'description':'Subset of PM2.5 & O3 AirNow observations and interpolated '+cams_str+' values withheld for validation'},
                        )

        ## Set the output paths & filenames
        out_dir = this_cams_dir
        if sites_vary:
            fname_use = out_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_use.nc')
            fname_val = out_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_val.nc')
        if sites_static:
            fname_use = out_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_static_'+this_yyyymmdd_hh+'00_use.nc')
            fname_val = out_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_static_'+this_yyyymmdd_hh+'00_val.nc')

        ## Write the datasets to NetCDF
        print('Writing file '+str(fname_use))
        ds_use.to_netcdf(fname_use)
        print('Writing file '+str(fname_val))
        ds_val.to_netcdf(fname_val)


if __name__ == '__main__':
    now_time_beg = dt.datetime.utcnow()
    date_beg, date_end, cams_fcst, cams_eac4 = parse_args()
    main(date_beg, date_end, cams_fcst, cams_eac4)
    now_time_end = dt.datetime.utcnow()
    run_time_tot = now_time_end - now_time_beg
    now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
    now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
    print('\nScript completed successfully.')
    print('   Beg time: '+now_time_beg_str)
    print('   End time: '+now_time_end_str)
    print('   Run time: '+str(run_time_tot)+'\n')
