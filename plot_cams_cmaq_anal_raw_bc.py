'''
plot_cams_cmaq_anal_raw_bc.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 16 Oct 2023

This script will plot maps of requested CAMS reanalysis, CAMS fcst, or CMAQ fcst "analysis" files,
either the raw or bias-corrected (BM3) versions.
'''

import sys
import pathlib
import argparse
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import wrf
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
from functions import map_functions

def make_o3_plot(fname, var, suptitle, suptitle_y, title_r,
        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
        mark1_lon=None, mark1_lat=None, mark1_var=None, mark1=None, mark1_col=None,
        mark2_lon=None, mark2_lat=None, mark2_var=None, mark2=None, mark2_col=None,
        mark_fill=False, marksize=None, lg_text=None, lg_loc=None, plot_missing_stations=False):

    mpl_o3   = '$\mathregular{O_3}$'
    min_o3 = 0.0
    max_o3 = 70.1
    int_o3 = 5.0

    var_min = np.nanmin(var)
    var_max = np.nanmax(var)
    var_name = 'Surface-level '+mpl_o3+' Conc.'
    var_cbar = 'Ozone Concentration'
    var_unit = 'ppbv'
    cbar_lab = var_cbar+' ['+var_unit+']'
    extend = 'max'
    #cmap = mpl.cm.get_cmap('rainbow').copy()
    cmap = mpl.colormaps.get_cmap('rainbow').copy()
    if plot_missing_stations:
        cmap.set_bad('white')
    else:
        # Remove stations with missing obs from the map
        inds1 = np.where(~np.isnan(mark1_var))[0]
        inds2 = np.where(~np.isnan(mark2_var))[0]
        mark1_lon = mark1_lon[inds1]
        mark1_lat = mark1_lat[inds1]
        mark1_var = mark1_var[inds1]
        mark2_lon = mark2_lon[inds2]
        mark2_lat = mark2_lat[inds2]
        mark2_var = mark2_var[inds2]
    bounds = np.arange(min_o3, max_o3, int_o3)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)
    title_l = var_name+f'\nMin: {var_min:.1f} '+var_unit+f', Max: {var_max:.1f} '+var_unit

    map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
        cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
        borders=borders, states=states, lakes=lakes, water_color='none',
        title_l=title_l, title_r=title_r, suptitle_y=suptitle_y,
        marker_lon=mark1_lon, marker_lat=mark1_lat, marker_var=mark1_var, marker_val_fill=mark_fill,
        marker=mark1, marker_color=mark1_col, marker_size=marksize,
        marker_lon2=mark2_lon, marker_lat2=mark2_lat, marker_var2=mark2_var, marker_val_fill2=mark_fill,
        marker2=mark2, marker_color2=mark2_col, marker_size2=marksize, marker_zorder=12,
        lg_text=lg_text, lg_loc=lg_loc, border_width=1.0, marker_width=0.75, marker_width2=0.75)
    fname.chmod(0o666)

def make_pm25_plot(fname, var, suptitle, suptitle_y, title_r,
        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
        mark1_lon=None, mark1_lat=None, mark1_var=None, mark1=None, mark1_col=None,
        mark2_lon=None, mark2_lat=None, mark2_var=None, mark2=None, mark2_col=None,
        mark_fill=False, marksize=None, lg_text=None, lg_loc=None, plot_missing_stations=False):

    mpl_pm25 = '$\mathregular{PM_{2.5}}$'
    mpl_ugm3 = '\u03bcg $\mathregular{m^{-3}}$'
    min_pm25 = 0.0
#    max_pm25 = 50.1
    max_pm25 = 30.1
    int_pm25 = 5.0
    cont_pm25 = np.arange(min_pm25, max_pm25, int_pm25)
#    cont_pm25 = np.append(cont_pm25, np.arange(60.0, 100.1, 10.0))
#    cont_pm25 = np.append(cont_pm25, np.arange(125.0, 200.1, 25.0))
    cont_pm25 = np.append(cont_pm25, np.arange(40.0, 50.1, 10.0))
    cont_pm25 = np.append(cont_pm25, np.arange(75.0, 100.1, 25.0))

    var_min = np.nanmin(var)
    var_max = np.nanmax(var)
    var_name = 'Surface-level '+mpl_pm25+' Conc.'
    var_cbar = mpl_pm25+' Concentration'
    var_unit = mpl_ugm3
    cbar_lab = var_cbar+' ['+var_unit+']'
    extend = 'max'
    #cmap = mpl.cm.get_cmap('rainbow').copy()
    cmap = mpl.colormaps.get_cmap('rainbow').copy()
    if plot_missing_stations:
        cmap.set_bad('white')
    else:
        # Remove stations with missing obs from the map
        inds1 = np.where(~np.isnan(mark1_var))[0]
        inds2 = np.where(~np.isnan(mark2_var))[0]
        mark1_lon = mark1_lon[inds1]
        mark1_lat = mark1_lat[inds1]
        mark1_var = mark1_var[inds1]
        mark2_lon = mark2_lon[inds2]
        mark2_lat = mark2_lat[inds2]
        mark2_var = mark2_var[inds2]
    bounds = cont_pm25
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)
    title_l = var_name+f'\nMin: {var_min:.1f} '+var_unit+f', Max: {var_max:.1f} '+var_unit

    map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
        cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
        borders=borders, states=states, lakes=lakes, water_color='none',
        title_l=title_l, title_r=title_r, suptitle_y=suptitle_y,
        marker_lon=mark1_lon, marker_lat=mark1_lat, marker_var=mark1_var, marker_val_fill=mark_fill,
        marker=mark1, marker_color=mark1_col, marker_size=marksize,
        marker_lon2=mark2_lon, marker_lat2=mark2_lat, marker_var2=mark2_var, marker_val_fill2=mark_fill,
        marker2=mark2, marker_color2=mark2_col, marker_size2=marksize, marker_zorder=12,
        lg_text=lg_text, lg_loc=lg_loc, border_width=1.0, marker_width=0.75, marker_width2=0.75)
    fname.chmod(0o666)

def main():
    beg_date_str = '2020-09-14_21'
    end_date_str = '2020-09-14_21'
    plot_stride_h = 3

    plot_type = 'png'
    plot_o3 = True
    plot_pm25 = True
    plot_cams_ra_raw = True
    plot_cams_fc_raw = True
    plot_cams_fc_bm3 = True
    plot_cmaq_fc_raw = True
    plot_cmaq_fc_bm3 = True
    plot_airnow_stns = False
    plot_airnow_vals = True
    plot_use = True
    plot_val = True
    plot_subdomain = False

    if plot_subdomain:
        i_beg, i_end = 20, 60
        j_beg, j_end = 130, 170
        lat_labels = [36, 37, 38, 39, 40, 41]
        lon_labels = [-124, -123, -122, -121, -120, -119, -118]
        suptitle_y = 1.01
        lg_loc = 'lower right'
        marksize = 81
    else:
        i_beg, i_end = 0, -1
        j_beg, j_end = 0, -1
        lat_labels = [25, 30, 35, 40, 45, 50]
        lon_labels = [-130, -120, -110, -110, -100, -90, -80, -70, -60]
        suptitle_y = 0.91
        lg_loc = 'lower left'
        marksize = 49

    cams_ra_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4')
    cams_fc_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')
    cmaq_fc_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','cmaq_hourly')
    cams_fc_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','BM3','static','cams_include0')
    cmaq_fc_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','BM3','static','cmaq_include0_nn')
    obs_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static','20200801_00-20211231_23')
    out_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','plots','anal')

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

    ## Create datetime array of times to plot
    dt_beg_valid = dt.datetime.strptime(beg_date_str, fmt_date_hh)
    dt_end_valid = dt.datetime.strptime(end_date_str, fmt_date_hh)
    dt_all_valid = pd.date_range(start=dt_beg_valid, end=dt_end_valid, freq=str(plot_stride_h)+'h')
    n_valid_times = len(dt_all_valid)

    ## Open a CMAQ sample/coordinate file to get map projection parameters
    cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
    print('Reading CMAQ coordinate data from '+str(cmaq_fname))
    cmaq_ds = xr.open_dataset(cmaq_fname)
    if plot_subdomain:
        cmaq_lon = cmaq_ds.LON[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
        cmaq_lat = cmaq_ds.LAT[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
    else:
        cmaq_lon = cmaq_ds.LON[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})
        cmaq_lat = cmaq_ds.LAT[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})
    cmaq_lon = cmaq_lon.assign_coords(coords={'XLAT':cmaq_lat, 'XLONG':cmaq_lon})
    cmaq_lat = cmaq_lat.assign_coords(coords={'XLAT':cmaq_lat, 'XLONG':cmaq_lon})
    cmaq_lon.attrs['long_name'] = 'longitude'
    cmaq_lat.attrs['long_name'] = 'latitude'
    cmaq_lon.attrs['units'] = 'degrees_east'
    cmaq_lat.attrs['units'] = 'degrees_north'
    n_cmaq_lon = cmaq_lon.shape[1]
    n_cmaq_lat = cmaq_lat.shape[0]
    truelat1 = cmaq_ds.attrs['P_ALP']
    truelat2 = cmaq_ds.attrs['P_BET']
    stand_lon = cmaq_ds.attrs['P_GAM']
    cen_lat = cmaq_ds.attrs['YCENT']
    cen_lon = cmaq_ds.attrs['XCENT']
    pole_lat = 90.0
    pole_lon = 0.0
    map_proj = 1
    moad_cen_lat = cen_lat
    dx = cmaq_ds.attrs['XCELL']
    dy = cmaq_ds.attrs['YCELL']

    ## Create a dictionary of these projection attributes
    dict_proj = {
        'MAP_PROJ':map_proj, 'CEN_LAT':cen_lat, 'CEN_LON':cen_lon,
        'TRUELAT1':truelat1, 'TRUELAT2':truelat2, 'MOAD_CEN_LAT':moad_cen_lat,
        'STAND_LON':stand_lon, 'POLE_LAT':pole_lat, 'POLE_LON':pole_lon, 'DX':dx, 'DY':dy,
        }

    ## Create an object of class wrf.WrfProj, then a cartopy mapping object
    ## This is essentially a manual reproduction of what wrf.get_cartopy() does
    wrf_proj = wrf.util.getproj(**dict_proj)
    proj_obj = wrf_proj.cartopy()
    cart_proj = proj_obj

    ## Now get the cartopy project x and y limits
    cmaq_lat  = cmaq_lat.assign_attrs(projection=wrf_proj)
    cmaq_lon  = cmaq_lon.assign_attrs(projection=wrf_proj)
    cart_xlim = wrf.cartopy_xlim(var=cmaq_lat)
    cart_ylim = wrf.cartopy_ylim(var=cmaq_lat)

    ## Get Cartopy features
    borders, states, oceans, lakes, rivers, land = map_functions.get_cartopy_features()

    ## Loop over valid times
    for tt in range(n_valid_times):
        this_dt = dt_all_valid[tt]
        this_date = this_dt.strftime(fmt_yyyymmdd)
        this_yr = this_dt.strftime(fmt_yyyy)
        this_mo = this_dt.strftime(fmt_mm)
        this_hr = this_dt.strftime(fmt_hh)
        this_date_file = this_dt.strftime(fmt_date_file)
        this_date_plot = this_dt.strftime(fmt_date_plot)

        ## Do we want to add AirNow station locations or values to the maps?
        if plot_airnow_stns or plot_airnow_vals:
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

            lg_text = ['AirNow Training', 'AirNow Validation']
            mark1 = 'v'
            mark2 = 'o'

            if plot_airnow_vals or plot_airnow_stns:
                if plot_use:
                    mark1_lon_o3 = o3_lon_use
                    mark1_lat_o3 = o3_lat_use
                    mark1_lon_pm25 = pm25_lon_use
                    mark1_lat_pm25 = pm25_lat_use
                else:
                    mark1_lon_o3 = None
                    mark1_lat_o3 = None
                    mark1_lon_pm25 = None
                    mark1_lat_pm25 = None
                if plot_val:
                    mark2_lon_o3 = o3_lon_val
                    mark2_lat_o3 = o3_lat_val
                    mark2_lon_pm25 = pm25_lon_val
                    mark2_lat_pm25 = pm25_lat_val
                else:
                    mark2_lon_o3 = None
                    mark2_lat_o3 = None
                    mark2_lon_pm25 = None
                    mark2_lat_pm25 = None

                mark1_col_o3 = 'None'
                mark2_col_o3 = 'None'
                mark1_col_pm25 = 'None'
                mark2_col_pm25 = 'None'

                if plot_airnow_vals:
                    mark_fill = True
                    mark1_var_o3 = o3_con_use
                    mark2_var_o3 = o3_con_val
                    mark1_var_pm25 = pm25_con_use
                    mark2_var_pm25 = pm25_con_val
                else:
                    mark_fill = False
                    mark1_var_o3 = None
                    mark2_var_o3 = None
                    mark1_var_pm25 = None
                    mark2_var_pm25 = None
            else:
                mark1_lon_o3 = None
                mark1_lat_o3 = None
                mark2_lon_o3 = None
                mark2_lat_o3 = None
                mark1_var_o3 = None
                mark2_var_o3 = None
                mark1_col_o3 = None
                mark2_col_o3 = None
                mark1_lon_pm25 = None
                mark1_lat_pm25 = None
                mark2_lon_pm25 = None
                mark2_lat_pm25 = None
                mark1_var_pm25 = None
                mark2_var_pm25 = None
                mark1_col_pm25 = None
                mark2_col_pm25 = None
        else:
            lg_text = None
            mark1 = None
            mark2 = None
            mark1_lon_o3 = None
            mark1_lat_o3 = None
            mark2_lon_o3 = None
            mark2_lat_o3 = None
            mark1_col_o3 = None
            mark2_col_o3 = None
            mark1_var_o3 = None
            mark2_var_o3 = None
            mark1_lon_pm25 = None
            mark1_lat_pm25 = None
            mark2_lon_pm25 = None
            mark2_lat_pm25 = None
            mark1_col_pm25 = None
            mark2_col_pm25 = None
            mark1_var_pm25 = None
            mark2_var_pm25 = None
            mark_fill = False
                    

        if plot_cams_ra_raw:
            this_dir = cams_ra_raw_dir.joinpath(this_yr, this_mo)
            fname_in = this_dir.joinpath('cams_o3_pm25_regrid_cmaq_'+this_date_file+'.nc')

            suptitle = 'CAMS Reanalysis Raw, Regridded to CMAQ Grid'
            title_r  = 'Valid: '+this_date_plot+' UTC'

            if fname_in.is_file():
                print('Reading '+str(fname_in))
                ds_xr = xr.open_dataset(fname_in)
                cams_ra_raw_o3 = ds_xr.cams_o3[0,:,:]
                cams_ra_raw_o3 = cams_ra_raw_o3.rename({'latitude':'XLAT', 'longitude':'XLONG'})
                cams_ra_raw_o3 = cams_ra_raw_o3.assign_attrs(projection=wrf_proj)
                cams_ra_raw_pm25 = ds_xr.cams_pm25[0,:,:]
                cams_ra_raw_pm25 = cams_ra_raw_pm25.rename({'latitude':'XLAT', 'longitude':'XLONG'})
                cams_ra_raw_pm25 = cams_ra_raw_pm25.assign_attrs(projection=wrf_proj)
                var_o3   = cams_ra_raw_o3
                var_pm25 = cams_ra_raw_pm25

                out_dir = out_dir_parent.joinpath('cams_ra_raw', this_yr, this_mo)
                out_dir.mkdir(parents=True, exist_ok=True)
                out_dir.chmod(0o777)
                fname_plot_o3   = out_dir.joinpath('cams_ra_raw_o3_'+this_date_file+'.'+plot_type)
                fname_plot_pm25 = out_dir.joinpath('cams_ra_raw_pm25_'+this_date_file+'.'+plot_type)
                if plot_airnow_vals:
                    fname_plot_o3   = out_dir.joinpath('cams_ra_raw_airnow_vals_o3_'+this_date_file+'.'+plot_type)
                    fname_plot_pm25 = out_dir.joinpath('cams_ra_raw_airnow_vals_pm25_'+this_date_file+'.'+plot_type)
                elif plot_airnow_stns:
                    fname_plot_o3   = out_dir.joinpath('cams_ra_raw_airnow_stns_o3_'+this_date_file+'.'+plot_type)
                    fname_plot_pm25 = out_dir.joinpath('cams_ra_raw_airnow_stns_pm25_'+this_date_file+'.'+plot_type)

                if plot_o3:
                    make_o3_plot(fname_plot_o3, var_o3, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_o3, mark1_lat=mark1_lat_o3, mark1_var=mark1_var_o3,
                        mark1=mark1, mark1_col=mark1_col_o3,
                        mark2_lon=mark2_lon_o3, mark2_lat=mark2_lat_o3, mark2_var=mark2_var_o3,
                        mark2=mark2, mark2_col=mark2_col_o3,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)

                if plot_pm25:
                    make_pm25_plot(fname_plot_pm25, var_pm25, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_pm25, mark1_lat=mark1_lat_pm25, mark1_var=mark1_var_pm25,
                        mark1=mark1, mark1_col=mark1_col_pm25,
                        mark2_lon=mark2_lon_pm25, mark2_lat=mark2_lat_pm25, mark2_var=mark2_var_pm25,
                        mark2=mark2, mark2_col=mark2_col_pm25,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)
            else:
                print('WARNING: File not found: '+str(fname))

        if plot_cams_fc_raw:
            if int(this_hr) < 12:
                this_cycle = this_date+'_00'
            else:
                this_cycle = this_date+'_12'

            this_dir = cams_fc_raw_dir.joinpath(this_yr, this_mo, this_cycle)
            this_cycle_dt = pd.to_datetime(this_cycle, format=fmt_yyyymmdd_hh)
            this_cycle_plot = this_cycle_dt.strftime(fmt_date_plot)+' UTC'
            fname = this_dir.joinpath('cams_o3_pm25_regrid_cmaq_'+this_date_file+'.nc')

            suptitle = 'CAMS Forecast Raw, Regridded to CMAQ Grid'
            title_r  = 'Start: '+this_cycle_plot+'\nValid: '+this_date_plot+' UTC'

            if fname.is_file():
                print('Reading '+str(fname))
                ds_xr = xr.open_dataset(fname)
                cams_fc_raw_o3 = ds_xr.cams_o3[0,:,:]
                cams_fc_raw_o3 = cams_fc_raw_o3.rename({'latitude':'XLAT', 'longitude':'XLONG'})
                cams_fc_raw_o3 = cams_fc_raw_o3.assign_attrs(projection=wrf_proj)
                cams_fc_raw_pm25 = ds_xr.cams_pm25[0,:,:]
                cams_fc_raw_pm25 = cams_fc_raw_pm25.rename({'latitude':'XLAT', 'longitude':'XLONG'})
                cams_fc_raw_pm25 = cams_fc_raw_pm25.assign_attrs(projection=wrf_proj)
                var_o3   = cams_fc_raw_o3
                var_pm25 = cams_fc_raw_pm25

                out_dir = out_dir_parent.joinpath('cams_fc_raw', this_yr, this_mo)
                out_dir.mkdir(parents=True, exist_ok=True)
                out_dir.chmod(0o777)
                fname_plot_o3   = out_dir.joinpath('cams_fc_raw_o3_'+this_date_file+'.'+plot_type)
                fname_plot_pm25 = out_dir.joinpath('cams_fc_raw_pm25_'+this_date_file+'.'+plot_type)
                if plot_airnow_vals:
                    fname_plot_o3   = out_dir.joinpath('cams_fc_raw_airnow_vals_o3_'+this_date_file+'.'+plot_type)
                    fname_plot_pm25 = out_dir.joinpath('cams_fc_raw_airnow_vals_pm25_'+this_date_file+'.'+plot_type)
                elif plot_airnow_stns:
                    fname_plot_o3   = out_dir.joinpath('cams_fc_raw_airnow_stns_o3_'+this_date_file+'.'+plot_type)
                    fname_plot_pm25 = out_dir.joinpath('cams_fc_raw_airnow_stns_pm25_'+this_date_file+'.'+plot_type)

                if plot_o3:
                    make_o3_plot(fname_plot_o3, var_o3, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_o3, mark1_lat=mark1_lat_o3, mark1_var=mark1_var_o3,
                        mark1=mark1, mark1_col=mark1_col_o3,
                        mark2_lon=mark2_lon_o3, mark2_lat=mark2_lat_o3, mark2_var=mark2_var_o3,
                        mark2=mark2, mark2_col=mark2_col_o3,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)

                if plot_pm25:
                    make_pm25_plot(fname_plot_pm25, var_pm25, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_pm25, mark1_lat=mark1_lat_pm25, mark1_var=mark1_var_pm25,
                        mark1=mark1, mark1_col=mark1_col_pm25,
                        mark2_lon=mark2_lon_pm25, mark2_lat=mark2_lat_pm25, mark2_var=mark2_var_pm25,
                        mark2=mark2, mark2_col=mark2_col_pm25,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)
            else:
                print('WARNING: File not found: '+str(fname))

        if plot_cams_fc_bm3:
            if int(this_hr) < 12:
                this_cycle = this_date+'_00'
            else:
                this_cycle = this_date+'_12'

            this_dir = cams_fc_bm3_dir
            this_cycle_dt = pd.to_datetime(this_cycle, format=fmt_yyyymmdd_hh)
            this_cycle_plot = this_cycle_dt.strftime(fmt_date_plot)+' UTC'
            fname_o3 = this_dir.joinpath('cams_o3_regrid_cmaq_'+this_date_file+'_BM3_static_include0.nc')
            fname_pm25 = this_dir.joinpath('cams_pm25_regrid_cmaq_'+this_date_file+'_BM3_static_noQC.nc')

            suptitle = 'CAMS Forecast BC, Regridded to CMAQ Grid'
            title_r  = 'Start: '+this_cycle_plot+'\nValid: '+this_date_plot+' UTC'

            if fname_o3.is_file():
                print('Reading '+str(fname_o3))
                ds_xr = xr.open_dataset(fname_o3)
                cams_fc_bm3_o3 = ds_xr.cams_o3_m
                cams_fc_bm3_o3 = cams_fc_bm3_o3.rename({'y':'XLAT', 'x':'XLONG'})
                cams_fc_bm3_o3 = cams_fc_bm3_o3.assign_attrs(projection=wrf_proj)
                var_o3   = cams_fc_bm3_o3

                out_dir = out_dir_parent.joinpath('cams_fc_bm3', this_yr, this_mo)
                out_dir.mkdir(parents=True, exist_ok=True)
                out_dir.chmod(0o777)
                fname_plot_o3   = out_dir.joinpath('cams_fc_bm3_o3_'+this_date_file+'.'+plot_type)
                if plot_airnow_vals:
                    fname_plot_o3   = out_dir.joinpath('cams_fc_bm3_airnow_vals_o3_'+this_date_file+'.'+plot_type)
                elif plot_airnow_stns:
                    fname_plot_o3   = out_dir.joinpath('cams_fc_bm3_airnow_stns_o3_'+this_date_file+'.'+plot_type)

                if plot_o3:
                    make_o3_plot(fname_plot_o3, var_o3, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_o3, mark1_lat=mark1_lat_o3, mark1_var=mark1_var_o3,
                        mark1=mark1, mark1_col=mark1_col_o3,
                        mark2_lon=mark2_lon_o3, mark2_lat=mark2_lat_o3, mark2_var=mark2_var_o3,
                        mark2=mark2, mark2_col=mark2_col_o3,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)

            else:
                print('WARNING: File not found: '+str(fname_o3))

            if fname_pm25.is_file():
                print('Reading '+str(fname_pm25))
                ds_xr = xr.open_dataset(fname_pm25)
                cams_fc_bm3_pm25 = ds_xr.cams_pm25_m
                cams_fc_bm3_pm25 = cams_fc_bm3_pm25.rename({'y':'XLAT', 'x':'XLONG'})
                cams_fc_bm3_pm25 = cams_fc_bm3_pm25.assign_attrs(projection=wrf_proj)
                var_pm25 = cams_fc_bm3_pm25

                out_dir = out_dir_parent.joinpath('cams_fc_bm3', this_yr, this_mo)
                out_dir.mkdir(parents=True, exist_ok=True)
                out_dir.chmod(0o777)
                fname_plot_pm25 = out_dir.joinpath('cams_fc_bm3_pm25_'+this_date_file+'.'+plot_type)
                if plot_airnow_vals:
                    fname_plot_pm25 = out_dir.joinpath('cams_fc_bm3_airnow_vals_pm25_'+this_date_file+'.'+plot_type)
                elif plot_airnow_stns:
                    fname_plot_pm25 = out_dir.joinpath('cams_fc_bm3_airnow_stns_pm25_'+this_date_file+'.'+plot_type)

                if plot_pm25:
                    make_pm25_plot(fname_plot_pm25, var_pm25, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_pm25, mark1_lat=mark1_lat_pm25, mark1_var=mark1_var_pm25,
                        mark1=mark1, mark1_col=mark1_col_pm25,
                        mark2_lon=mark2_lon_pm25, mark2_lat=mark2_lat_pm25, mark2_var=mark2_var_pm25,
                        mark2=mark2, mark2_col=mark2_col_pm25,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)
            else:
                print('WARNING: File not found: '+str(fname_pm25))

        if plot_cmaq_fc_raw:
            if int(this_hr) >= 0 and int(this_hr) < 6:
                this_cycle = this_date+'_00'
            elif int(this_hr) >= 6 and int(this_hr) < 12:
                this_cycle = this_date+'_06'
            elif int(this_hr) >= 12 and int(this_hr) < 18:
                this_cycle = this_date+'_12'
            elif int(this_hr) >= 18:
                this_cycle = this_date+'_18'

            this_dir = cmaq_fc_raw_dir
            this_cycle_dt = pd.to_datetime(this_cycle, format=fmt_yyyymmdd_hh)
            this_cycle_plot = this_cycle_dt.strftime(fmt_date_plot)+' UTC'
            fname = this_dir.joinpath('cmaq_airnow_pm2.5_o3_'+this_date+this_hr+'.nc')

            suptitle = 'CMAQ Forecast Raw'
            title_r  = 'Start: '+this_cycle_plot+'\nValid: '+this_date_plot+' UTC'

            if fname.is_file():
                print('Reading '+str(fname))
                ds_xr = xr.open_dataset(fname)
                cmaq_fc_raw_o3 = ds_xr.O3
                cmaq_fc_raw_o3 = cmaq_fc_raw_o3.rename({'lat':'XLAT', 'lon':'XLONG'})
                cmaq_fc_raw_o3 = cmaq_fc_raw_o3.assign_attrs(projection=wrf_proj)
                cmaq_fc_raw_pm25 = ds_xr.PM25
                cmaq_fc_raw_pm25 = cmaq_fc_raw_pm25.rename({'lat':'XLAT', 'lon':'XLONG'})
                cmaq_fc_raw_pm25 = cmaq_fc_raw_pm25.assign_attrs(projection=wrf_proj)
                var_o3   = cmaq_fc_raw_o3
                var_pm25 = cmaq_fc_raw_pm25

                out_dir = out_dir_parent.joinpath('cmaq_fc_raw', this_yr, this_mo)
                out_dir.mkdir(parents=True, exist_ok=True)
                out_dir.chmod(0o777)
                fname_plot_o3   = out_dir.joinpath('cmaq_fc_raw_o3_'+this_date_file+'.'+plot_type)
                fname_plot_pm25 = out_dir.joinpath('cmaq_fc_raw_pm25_'+this_date_file+'.'+plot_type)
                if plot_airnow_vals:
                    fname_plot_o3   = out_dir.joinpath('cmaq_fc_raw_airnow_vals_o3_'+this_date_file+'.'+plot_type)
                    fname_plot_pm25 = out_dir.joinpath('cmaq_fc_raw_airnow_vals_pm25_'+this_date_file+'.'+plot_type)
                elif plot_airnow_stns:
                    fname_plot_o3   = out_dir.joinpath('cmaq_fc_raw_airnow_stns_o3_'+this_date_file+'.'+plot_type)
                    fname_plot_pm25 = out_dir.joinpath('cmaq_fc_raw_airnow_stns_pm25_'+this_date_file+'.'+plot_type)

                if plot_o3:
                    make_o3_plot(fname_plot_o3, var_o3, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_o3, mark1_lat=mark1_lat_o3, mark1_var=mark1_var_o3,
                        mark1=mark1, mark1_col=mark1_col_o3,
                        mark2_lon=mark2_lon_o3, mark2_lat=mark2_lat_o3, mark2_var=mark2_var_o3,
                        mark2=mark2, mark2_col=mark2_col_o3,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)

                if plot_pm25:
                    make_pm25_plot(fname_plot_pm25, var_pm25, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_pm25, mark1_lat=mark1_lat_pm25, mark1_var=mark1_var_pm25,
                        mark1=mark1, mark1_col=mark1_col_pm25,
                        mark2_lon=mark2_lon_pm25, mark2_lat=mark2_lat_pm25, mark2_var=mark2_var_pm25,
                        mark2=mark2, mark2_col=mark2_col_pm25,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)
            else:
                print('WARNING: File not found: '+str(fname))

        if plot_cmaq_fc_bm3:
            if int(this_hr) >= 0 and int(this_hr) < 6:
                this_cycle = this_date+'_00'
            elif int(this_hr) >= 6 and int(this_hr) < 12:
                this_cycle = this_date+'_06'
            elif int(this_hr) >= 12 and int(this_hr) < 18:
                this_cycle = this_date+'_12'
            elif int(this_hr) >= 18:
                this_cycle = this_date+'_18'

            this_dir = cmaq_fc_bm3_dir
            this_cycle_dt = pd.to_datetime(this_cycle, format=fmt_yyyymmdd_hh)
            this_cycle_plot = this_cycle_dt.strftime(fmt_date_plot)+' UTC'
            fname_o3 = this_dir.joinpath('cmaq_o3_'+this_date_file+'_BM3_static_inc0_nn.nc')
            fname_pm25 = this_dir.joinpath('cmaq_pm25_'+this_date_file+'_BM3_static_inc0_nn.nc')

            suptitle = 'CMAQ Forecast BC'
            title_r  = 'Start: '+this_cycle_plot+'\nValid: '+this_date_plot+' UTC'

            if fname_o3.is_file():
                print('Reading '+str(fname_o3))
                ds_xr = xr.open_dataset(fname_o3)
                cmaq_fc_bm3_o3 = ds_xr.cmaq_o3_m
                cmaq_fc_bm3_o3 = cmaq_fc_bm3_o3.rename({'y':'XLAT', 'x':'XLONG'})
                cmaq_fc_bm3_o3 = cmaq_fc_bm3_o3.assign_attrs(projection=wrf_proj)
                var_o3   = cmaq_fc_bm3_o3

                out_dir = out_dir_parent.joinpath('cmaq_fc_bm3', this_yr, this_mo)
                out_dir.mkdir(parents=True, exist_ok=True)
                out_dir.chmod(0o777)
                fname_plot_o3   = out_dir.joinpath('cmaq_fc_bm3_o3_'+this_date_file+'.'+plot_type)
                if plot_airnow_vals:
                    fname_plot_o3   = out_dir.joinpath('cmaq_fc_bm3_airnow_vals_o3_'+this_date_file+'.'+plot_type)
                elif plot_airnow_stns:
                    fname_plot_o3   = out_dir.joinpath('cmaq_fc_bm3_airnow_stns_o3_'+this_date_file+'.'+plot_type)

                if plot_o3:
                    make_o3_plot(fname_plot_o3, var_o3, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_o3, mark1_lat=mark1_lat_o3, mark1_var=mark1_var_o3,
                        mark1=mark1, mark1_col=mark1_col_o3,
                        mark2_lon=mark2_lon_o3, mark2_lat=mark2_lat_o3, mark2_var=mark2_var_o3,
                        mark2=mark2, mark2_col=mark2_col_o3,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)

            else:
                print('WARNING: File not found: '+str(fname_o3))

            if fname_pm25.is_file():
                print('Reading '+str(fname_pm25))
                ds_xr = xr.open_dataset(fname_pm25)
                cmaq_fc_bm3_pm25 = ds_xr.cmaq_pm25_m
                cmaq_fc_bm3_pm25 = cmaq_fc_bm3_pm25.rename({'y':'XLAT', 'x':'XLONG'})
                cmaq_fc_bm3_pm25 = cmaq_fc_bm3_pm25.assign_attrs(projection=wrf_proj)
                var_pm25 = cmaq_fc_bm3_pm25

                out_dir = out_dir_parent.joinpath('cmaq_fc_bm3', this_yr, this_mo)
                out_dir.mkdir(parents=True, exist_ok=True)
                out_dir.chmod(0o777)
                fname_plot_pm25 = out_dir.joinpath('cmaq_fc_bm3_pm25_'+this_date_file+'.'+plot_type)
                if plot_airnow_vals:
                    fname_plot_pm25 = out_dir.joinpath('cmaq_fc_bm3_airnow_vals_pm25_'+this_date_file+'.'+plot_type)
                elif plot_airnow_stns:
                    fname_plot_pm25 = out_dir.joinpath('cmaq_fc_bm3_airnow_stns_pm25_'+this_date_file+'.'+plot_type)

                if plot_pm25:
                    make_pm25_plot(fname_plot_pm25, var_pm25, suptitle, suptitle_y, title_r,
                        cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat, borders, states, lakes,
                        mark1_lon=mark1_lon_pm25, mark1_lat=mark1_lat_pm25, mark1_var=mark1_var_pm25,
                        mark1=mark1, mark1_col=mark1_col_pm25,
                        mark2_lon=mark2_lon_pm25, mark2_lat=mark2_lat_pm25, mark2_var=mark2_var_pm25,
                        mark2=mark2, mark2_col=mark2_col_pm25,
                        mark_fill=mark_fill, marksize=marksize, lg_text=lg_text, lg_loc=lg_loc)
            else:
                print('WARNING: File not found: '+str(fname_pm25))




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
