'''
regrid_cams_to_cmaq.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 8 Feb 2022

This script reads in CAMS forecast or reanalysis data files and a single CMAQ forecast/coordinates file,
regrids the CAMS forecast/reanalysis data to the CMAQ grid, and writes out a new NetCDF file.

May add options for plotting the original or regridded output for convenience (and sanity checking).
'''

import os
import sys
import pathlib
import datetime as dt
import pandas as pd
import numpy as np
import xarray as xr
import scipy as sp
import scipy.interpolate
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import matplotlib.ticker as mticker
#import cartopy
#import cartopy.crs as ccrs
#import geocat.comp
from geocat.f2py import rgrid2rcm
#from functions import map_functions

def main():
    beg_date = '2020-08-03_00'
    end_date = '2021-12-31_23'

    cams_fcst = False
    cams_eac4 = True

    ## Plotting options don't yet do anything
    plot_orig   = False
    plot_regrid = False
    plot_o3   = False
    plot_pm25 = False
    plot_type = 'png'

    if cams_fcst and cams_eac4:
        print('ERROR: Both cams_fcst and cams_eac4 set to True. Please only set one of them to True.')
        sys.exit()
    if not cams_fcst and not cams_eac4:
        print('ERROR: Both cams_fcst and cams_eac4 set to False. Please set one of them to True.')
        sys.exit()

    fmt_date = '%Y-%m-%d'
    fmt_date_hh = '%Y-%m-%d_%H'
    fmt_date_wrf = '%Y-%m-%d_%H:%M:%S'
    fmt_yyyymmdd = '%Y%m%d'
    fmt_yyyymmdd_hh = '%Y%m%d_%H'
    fmt_ddmmyyyyhhmm = '%d %b %Y/%H%M'

    mpl_pm25 = '$\mathregular{PM_{2.5}}$'
    mpl_ugm3 = '\u03bcg $\mathregular{m^{-3}}$'

    min_o3 = 0.0
    max_o3 = 100.0
    int_o3 = 5.0

    min_pm25 = 0.0
    max_pm25 = 150.0
    int_pm25 = 5.0

    ## Build array of dates
    beg_dt = pd.to_datetime(beg_date, format=fmt_date_hh)
    end_dt = pd.to_datetime(end_date, format=fmt_date_hh)
    dt_range = pd.date_range(start=beg_dt, end=end_dt, freq='1h')
    n_times  = len(dt_range)

    ## Get Cartopy plotting features if any plots are desired
    if plot_o3 or plot_pm25:
        borders, states, oceans, lakes, rivers, land = map_functions.get_cartopy_features()
        suptitle_y = 0.77

    ## Open the CMAQ sample/coordinate file
    cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
    print('Reading CMAQ coordinate data from '+str(cmaq_fname))
    cmaq_DS = xr.open_dataset(cmaq_fname)
    cmaq_lon = cmaq_DS.LON[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})   # range (-180, 180]
    cmaq_lat = cmaq_DS.LAT[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})
    cmaq_lon.attrs['long_name'] = 'longitude'
    cmaq_lat.attrs['long_name'] = 'latitude'
    cmaq_lon.attrs['units'] = 'degrees_east'
    cmaq_lat.attrs['units'] = 'degrees_north'
    n_cmaq_lon = cmaq_lon.shape[1]
    n_cmaq_lat = cmaq_lat.shape[0]
    truelat1 = cmaq_DS.attrs['P_ALP']
    truelat2 = cmaq_DS.attrs['P_BET']
    stand_lon = cmaq_DS.attrs['P_GAM']
    cen_lat = cmaq_DS.attrs['YCENT']
    cen_lon = cmaq_DS.attrs['XCENT']
    pole_lat = 90.0
    pole_lon = 0.0
    map_proj = 1
    moad_cen_lat = cen_lat
    dx = cmaq_DS.attrs['XCELL']
    dy = cmaq_DS.attrs['YCELL']
    get_cams_lat_lon = True

    ## Loop through dates looking for CAMS files
    for tt in range(n_times):
        this_dt = dt_range[tt]
        this_yr = this_dt.strftime('%Y')
        this_mo = this_dt.strftime('%m')
        this_hr = this_dt.strftime('%H')
        this_yyyymmdd = this_dt.strftime(fmt_yyyymmdd)
        this_date_str = this_dt.strftime(fmt_date)
        this_date_hh_str = this_dt.strftime(fmt_date_hh)
        this_date_plot = this_dt.strftime(fmt_ddmmyyyyhhmm)+' UTC'
        this_date_wrf = this_dt.strftime(fmt_date_wrf)

        print('Processing CAMS data from '+this_date_hh_str)

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

        prev2_date_str = prev2_dt.strftime(fmt_date)
        prev1_date_str = prev1_dt.strftime(fmt_date)
        next1_date_str = next1_dt.strftime(fmt_date)
        next2_date_str = next2_dt.strftime(fmt_date)

        prev2_date_hh_str = prev2_dt.strftime(fmt_date_hh)
        prev1_date_hh_str = prev1_dt.strftime(fmt_date_hh)
        next1_date_hh_str = next1_dt.strftime(fmt_date_hh)
        next2_date_hh_str = next2_dt.strftime(fmt_date_hh)

        if cams_fcst:
            cams_parent_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')
            if int(this_hr) < 12:
                this_cycle = this_yyyymmdd+'_00'
            else:
                this_cycle = this_yyyymmdd+'_12'

            this_dir = cams_parent_dir.joinpath(this_yr,this_mo,this_cycle)
            this_cycle_dt = pd.to_datetime(this_cycle, format=fmt_yyyymmdd_hh)
            this_cycle_plot = this_cycle_dt.strftime(fmt_ddmmyyyyhhmm)+' UTC'

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

            prev2_dir = cams_parent_dir.joinpath(prev2_yr,prev2_mo,prev2_cycle)
            prev1_dir = cams_parent_dir.joinpath(prev1_yr,prev1_mo,prev1_cycle)
            next1_dir = cams_parent_dir.joinpath(next1_yr,next1_mo,next1_cycle)
            next2_dir = cams_parent_dir.joinpath(next2_yr,next2_mo,next2_cycle)

        elif cams_eac4:
            cams_parent_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4')
            this_dir = cams_parent_dir.joinpath(this_yr,this_mo)

            prev2_dir = cams_parent_dir.joinpath(prev2_yr,prev2_mo)
            prev1_dir = cams_parent_dir.joinpath(prev1_yr,prev1_mo)
            next1_dir = cams_parent_dir.joinpath(next1_yr,next1_mo)
            next2_dir = cams_parent_dir.joinpath(next2_yr,next2_mo)

        ## Get file names
        o3_fname_in_this   = this_dir.joinpath('o3sfc_'+this_date_hh_str+'00.nc')
        pm25_fname_in_this = this_dir.joinpath('pm2.5_'+this_date_hh_str+'00.nc')

        o3_fname_in_prev2   = prev2_dir.joinpath('o3sfc_'+prev2_date_hh_str+'00.nc')
        o3_fname_in_prev1   = prev1_dir.joinpath('o3sfc_'+prev1_date_hh_str+'00.nc')
        o3_fname_in_next1   = next1_dir.joinpath('o3sfc_'+next1_date_hh_str+'00.nc')
        o3_fname_in_next2   = next2_dir.joinpath('o3sfc_'+next2_date_hh_str+'00.nc')
        pm25_fname_in_prev2 = prev2_dir.joinpath('pm2.5_'+prev2_date_hh_str+'00.nc')
        pm25_fname_in_prev1 = prev1_dir.joinpath('pm2.5_'+prev1_date_hh_str+'00.nc')
        pm25_fname_in_next1 = next1_dir.joinpath('pm2.5_'+next1_date_hh_str+'00.nc')
        pm25_fname_in_next2 = next2_dir.joinpath('pm2.5_'+next2_date_hh_str+'00.nc')

        ## Do these files exist? (True/False)
        isfile_this_o3  = o3_fname_in_this.is_file()
        isfile_prev2_o3 = o3_fname_in_prev2.is_file()
        isfile_prev1_o3 = o3_fname_in_prev1.is_file()
        isfile_next1_o3 = o3_fname_in_next1.is_file()
        isfile_next2_o3 = o3_fname_in_next2.is_file()
        isfile_this_pm25  = pm25_fname_in_this.is_file()
        isfile_prev2_pm25 = pm25_fname_in_prev2.is_file()
        isfile_prev1_pm25 = pm25_fname_in_prev1.is_file()
        isfile_next1_pm25 = pm25_fname_in_next1.is_file()
        isfile_next2_pm25 = pm25_fname_in_next2.is_file()

        ## Get CAMS ozone data
        ## Simplest case: The file exists, no interpolation needed
        if isfile_this_o3:
            ds_o3 = xr.open_dataset(o3_fname_in_this)
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

#           ds_o3.assign_coords(longitude=(cams_lon_vals))
#           ds_o3.assign_coords(longitude=(((ds_o3.longitude + 180) % 360) -180))
#           print(ds_o3)
#           sys.exit()

            ## Need to reverse the latitude dimension just as cams_lat was reversed
            cams_o3_this = ds_o3.go3[:,::-1,:]
            ## The longitude coordinate array is [0, 360), but I can't figure out how to reassign it.
            ## Examples @ https://xarray.pydata.org/en/stable/generated/xarray.DataArray.assign_coords.html don't work.
#           print(cams_o3_this)
#           print(cams_lon)
#           cams_o3_this.assign_coords(longitude=(((cams_o3_this.longitude + 180) % 360) - 180))
#           cams_o3_this.assign_coords({'longitude': (((cams_o3_this.longitude + 180) % 360) - 180) })
#           cams_o3_this.assign_coords(longitude2=('longitude', cams_lon_vals) )
#           cams_o3_this.assign_coords(longitude=(cams_lon_vals))
#           print(cams_o3_this)
#           sys.exit()

        ## Case 2: This time bounded by files at h-2 and h+1
        elif isfile_prev2_o3 and isfile_next1_o3:
            ds_o3_prev2 = xr.open_dataset(o3_fname_in_prev2)
            ds_o3_next1 = xr.open_dataset(o3_fname_in_next1)
            if get_cams_lat_lon:
                cams_lat = ds_o3_prev2.latitude[::-1]
                cams_lon = ds_o3_prev2.longitude
                cams_lon_vals = np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon)
                cams_lon = xr.DataArray(np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon),
                    coords={'longitude':cams_lon_vals}, attrs={'units':'degrees_east', 'long_name':'longitude'})
                get_cams_lat_lon = False
            cams_o3_prev2 = ds_o3_prev2.go3[:,::-1,:]
            cams_o3_next1 = ds_o3_next1.go3[:,::-1,:]

            ## Temporally interpolate the data
#           cams_o3_this = (1/3)*cams_o3_prev2 + (2/3)*cams_o3_next1
            cams_o3_this_vals = (1/3)*cams_o3_prev2.values + (2/3)*cams_o3_next1.values
            cams_o3_this = cams_o3_prev2.copy(deep=True, data=cams_o3_this_vals)

        ## Case 3: This time bounded by files at h-1 and h+2
        elif isfile_prev1_o3 and isfile_next2_o3:
            ds_o3_prev1 = xr.open_dataset(o3_fname_in_prev1)
            ds_o3_next2 = xr.open_dataset(o3_fname_in_next2)
            if get_cams_lat_lon:
                cams_lat = ds_o3_prev1.latitude[::-1]
                cams_lon = ds_o3_prev1.longitude
                cams_lon_vals = np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon)
                cams_lon = xr.DataArray(np.where(cams_lon > 180.0, cams_lon-360.0, cams_lon),
                    coords={'longitude':cams_lon_vals}, attrs={'units':'degrees_east', 'long_name':'longitude'})
                get_cams_lat_lon = False
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
            ds_pm25 = xr.open_dataset(pm25_fname_in_this)
#           print(ds_pm25)
#           sys.exit()

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
            ds_pm25_prev2 = xr.open_dataset(pm25_fname_in_prev2)
            ds_pm25_next1 = xr.open_dataset(pm25_fname_in_next1)
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
            ds_pm25_prev1 = xr.open_dataset(pm25_fname_in_prev1)
            ds_pm25_next2 = xr.open_dataset(pm25_fname_in_next2)
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
            ## Bug found 9 Oct 2023: First term was multiplied by 1/3 (incorrect) instead of by 2/3 (correct)
            cams_pm25_this_vals = (2/3)*cams_pm25_prev1.values + (1/3)*cams_pm25_next2.values
            cams_pm25_this = cams_pm25_prev1.copy(deep=True, data=cams_pm25_this_vals)

        else:
            print('ERROR: This time does not have PM2.5 files at surrounding times for temporal interp. Exiting!')
            sys.exit()

        ## Convert PM2.5 concentration from kg/m3 to ug/m3
        cams_pm25_this_ug = cams_pm25_this * 1e9

        ## Now regrid both the O3 and PM2.5 fields from CAMS to the CMAQ grid
        print('-- Regridding the CAMS data to the CMAQ grid')
        ## Try regridding with geocat-comp
        ## My support questions have been resolved: https://github.com/NCAR/geocat-comp/issues/200
        #cams_o3_ppb_regrid  = geocat.comp.rgrid2rcm(cams_lat.values, cams_lon.values, cams_o3_this_ppb.values, cmaq_lat, cmaq_lon)
        #cams_pm25_ug_regrid = geocat.comp.rgrid2rcm(cams_lat, cams_lon, cams_pm25_this_ug, cmaq_lat, cmaq_lon)
        cams_o3_ppb_regrid = rgrid2rcm(cams_lat, cams_lon, cams_o3_this_ppb, cmaq_lat, cmaq_lon)
        cams_pm25_ug_regrid = rgrid2rcm(cams_lat, cams_lon, cams_pm25_this_ug, cmaq_lat, cmaq_lon)

#       print('\no3 before regridding')
#       print(cams_o3_this_ppb)
#       print('\no3 after regridding')
#       print(cams_o3_ppb_regrid)
#       print('\npm2.5 before regridding')
#       print(cams_pm25_this_ug)
#       print('\npm2.5 after regridding')
#       print(cams_pm25_ug_regrid)
#       sys.exit()

        ## Write out the regridded data to NetCDF files
        if cams_eac4:
            desc_text = 'CAMS reanalysis data regridded to CMAQ grid'
        elif cams_fcst:
            desc_text = 'CAMS forecast data regridded to CMAQ grid'
        DS_cams_regrid = xr.Dataset(
            data_vars=dict(
                cams_o3=(['time', 'lat', 'lon'], cams_o3_ppb_regrid.data, {'units':'ppb', 'description':'lowest model level ozone concentration', 'coordinates': 'times latitude longitude'}),
                cams_pm25=(['time', 'lat', 'lon'], cams_pm25_ug_regrid.data, {'units':'ug/m3', 'description':'surface-level PM2.5 concentration', 'coordinates': 'times latitude longitude'}),
                ),
            coords=dict(
                latitude=(['lat', 'lon'], cmaq_lat.data, {'units':'degrees_north', 'coordinates':'longitude latitude'}),
                longitude=(['lat', 'lon'], cmaq_lon.data, {'units':'degrees_east', 'coordinates':'longitude latitude'}),
                times=(['time'], np.atleast_1d(this_dt.to_pydatetime()) ),
                ),
            attrs=dict(
                description=desc_text,
                MAP_PROJ=map_proj, CEN_LAT=cen_lat, CEN_LON=cen_lon, TRUELAT1=truelat1, TRUELAT2=truelat2,
                STAND_LON=stand_lon, MOAD_CEN_LAT=moad_cen_lat, POLE_LAT=pole_lat, POLE_LON=pole_lon, DX=dx, DY=dy))

        ## Put the regridded file in the proper directory
        fname = this_dir.joinpath('cams_o3_pm25_regrid_cmaq_'+this_dt.strftime(fmt_yyyymmdd_hh)+'00.nc')
        print('-- Writing '+str(fname))
        DS_cams_regrid.to_netcdf(path=fname, mode='w', format='NETCDF4')


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
