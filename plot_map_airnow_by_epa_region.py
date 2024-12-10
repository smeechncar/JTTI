'''
plot_map_airnow_by_epa_region.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 18 Oct 2023

This script plots AirNow O3 & PM2.5 stations on a map, coding them by EPA region.
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

def plot_stn_map(fname_o3, epa_code_o3, lat_o3, lon_o3, suptitle_o3, fname_pm25, epa_code_pm25, lat_pm25, lon_pm25, suptitle_pm25):
    mpl_o3   = '$\mathregular{O_3}$'
    mpl_pm25 = '$\mathregular{PM_{2.5}}$'
    lat_labels = [25, 30, 35, 40, 45, 50]
    lon_labels = [-130, -120, -110, -110, -100, -90, -80, -70, -60]
    suptitle_y = 0.86
    lg_loc = 'lower right'
    marksize = 49
    fontsize = 14
    lg_fontsize = 10
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
        'tab:brown', 'tab:olive', 'tab:gray', 'tab:pink', 'tab:cyan', 'black']
    markers = ['o', 'v', '^', 's', 'P', 'p', '*', 'H', 'X', 'D', '8']
    msize =   [49,  49,  49,  36,  49,  49,  81,  49,  49,  36,  49]
    lg_text = ['Region 1', 'Region 2', 'Region 3', 'Region 4', 'Region 5',
        'Region 6', 'Region 7', 'Region 8', 'Region 9', 'Region 10', 'Non-U.S.']
#    lg_text = ['Reg. 1', 'Reg. 2', 'Reg. 3', 'Reg. 4', 'Reg. 5', 'Reg. 6', 'Reg. 7', 'Reg. 8', 'Reg. 9', 'Reg. 10', 'Non-U.S.']

    ## Open a CMAQ sample/coordinate file to get map projection parameters
    cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
    print('Reading CMAQ coordinate data from '+str(cmaq_fname))
    cmaq_ds = xr.open_dataset(cmaq_fname)
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

    ## Get indices for obs in each EPA region
    ind01_o3 = np.where(epa_code_o3 == 1)[0]
    ind02_o3 = np.where(epa_code_o3 == 2)[0]
    ind03_o3 = np.where(epa_code_o3 == 3)[0]
    ind04_o3 = np.where(epa_code_o3 == 4)[0]
    ind05_o3 = np.where(epa_code_o3 == 5)[0]
    ind06_o3 = np.where(epa_code_o3 == 6)[0]
    ind07_o3 = np.where(epa_code_o3 == 7)[0]
    ind08_o3 = np.where(epa_code_o3 == 8)[0]
    ind09_o3 = np.where(epa_code_o3 == 9)[0]
    ind10_o3 = np.where(epa_code_o3 == 10)[0]
    ind11_o3 = np.where(epa_code_o3 == 99)[0]

    ind01_pm25 = np.where(epa_code_pm25 == 1)[0]
    ind02_pm25 = np.where(epa_code_pm25 == 2)[0]
    ind03_pm25 = np.where(epa_code_pm25 == 3)[0]
    ind04_pm25 = np.where(epa_code_pm25 == 4)[0]
    ind05_pm25 = np.where(epa_code_pm25 == 5)[0]
    ind06_pm25 = np.where(epa_code_pm25 == 6)[0]
    ind07_pm25 = np.where(epa_code_pm25 == 7)[0]
    ind08_pm25 = np.where(epa_code_pm25 == 8)[0]
    ind09_pm25 = np.where(epa_code_pm25 == 9)[0]
    ind10_pm25 = np.where(epa_code_pm25 == 10)[0]
    ind11_pm25 = np.where(epa_code_pm25 == 99)[0]

    lon01_o3 = lon_o3[ind01_o3]
    lon02_o3 = lon_o3[ind02_o3]
    lon03_o3 = lon_o3[ind03_o3]
    lon04_o3 = lon_o3[ind04_o3]
    lon05_o3 = lon_o3[ind05_o3]
    lon06_o3 = lon_o3[ind06_o3]
    lon07_o3 = lon_o3[ind07_o3]
    lon08_o3 = lon_o3[ind08_o3]
    lon09_o3 = lon_o3[ind09_o3]
    lon10_o3 = lon_o3[ind10_o3]
    lon11_o3 = lon_o3[ind11_o3]

    lat01_o3 = lat_o3[ind01_o3]
    lat02_o3 = lat_o3[ind02_o3]
    lat03_o3 = lat_o3[ind03_o3]
    lat04_o3 = lat_o3[ind04_o3]
    lat05_o3 = lat_o3[ind05_o3]
    lat06_o3 = lat_o3[ind06_o3]
    lat07_o3 = lat_o3[ind07_o3]
    lat08_o3 = lat_o3[ind08_o3]
    lat09_o3 = lat_o3[ind09_o3]
    lat10_o3 = lat_o3[ind10_o3]
    lat11_o3 = lat_o3[ind11_o3]

    lon01_pm = lon_pm25[ind01_pm25]
    lon02_pm = lon_pm25[ind02_pm25]
    lon03_pm = lon_pm25[ind03_pm25]
    lon04_pm = lon_pm25[ind04_pm25]
    lon05_pm = lon_pm25[ind05_pm25]
    lon06_pm = lon_pm25[ind06_pm25]
    lon07_pm = lon_pm25[ind07_pm25]
    lon08_pm = lon_pm25[ind08_pm25]
    lon09_pm = lon_pm25[ind09_pm25]
    lon10_pm = lon_pm25[ind10_pm25]
    lon11_pm = lon_pm25[ind11_pm25]

    lat01_pm = lat_pm25[ind01_pm25]
    lat02_pm = lat_pm25[ind02_pm25]
    lat03_pm = lat_pm25[ind03_pm25]
    lat04_pm = lat_pm25[ind04_pm25]
    lat05_pm = lat_pm25[ind05_pm25]
    lat06_pm = lat_pm25[ind06_pm25]
    lat07_pm = lat_pm25[ind07_pm25]
    lat08_pm = lat_pm25[ind08_pm25]
    lat09_pm = lat_pm25[ind09_pm25]
    lat10_pm = lat_pm25[ind10_pm25]
    lat11_pm = lat_pm25[ind11_pm25]

#    suptitle_o3 = 'AirNow '+mpl_o3+' Stations by EPA Region'
#    suptitle_pm25 = 'AirNow '+mpl_pm25+' Stations by EPA Region'

    map_functions.mpl_map_stations(fname_o3, suptitle_o3, cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat,
        borders=borders, states=states, oceans=oceans, lakes=lakes, water_color='lightblue',
        lat_labels=lat_labels, lon_labels=lon_labels, suptitle_y=suptitle_y,
        lg_text=lg_text, lg_loc=lg_loc, lg_fontsize=lg_fontsize,
        mark1_lon=lon01_o3, mark1_lat=lat01_o3, mark1=markers[0], mark1_size=msize[0], mark1_color=colors[0],
        mark2_lon=lon02_o3, mark2_lat=lat02_o3, mark2=markers[1], mark2_size=msize[1], mark2_color=colors[1],
        mark3_lon=lon03_o3, mark3_lat=lat03_o3, mark3=markers[2], mark3_size=msize[2], mark3_color=colors[2],
        mark4_lon=lon04_o3, mark4_lat=lat04_o3, mark4=markers[3], mark4_size=msize[3], mark4_color=colors[3],
        mark5_lon=lon05_o3, mark5_lat=lat05_o3, mark5=markers[4], mark5_size=msize[4], mark5_color=colors[4],
        mark6_lon=lon06_o3, mark6_lat=lat06_o3, mark6=markers[5], mark6_size=msize[5], mark6_color=colors[5],
        mark7_lon=lon07_o3, mark7_lat=lat07_o3, mark7=markers[6], mark7_size=msize[6], mark7_color=colors[6],
        mark8_lon=lon08_o3, mark8_lat=lat08_o3, mark8=markers[7], mark8_size=msize[7], mark8_color=colors[7],
        mark9_lon=lon09_o3, mark9_lat=lat09_o3, mark9=markers[8], mark9_size=msize[8], mark9_color=colors[8],
        mark10_lon=lon10_o3, mark10_lat=lat10_o3, mark10=markers[9], mark10_size=msize[9], mark10_color=colors[9],
        mark11_lon=lon11_o3, mark11_lat=lat11_o3, mark11=markers[10],mark11_size=msize[10], mark11_color=colors[10],
        )
    fname_o3.chmod(0o666)

    map_functions.mpl_map_stations(fname_pm25, suptitle_pm25, cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat,
        borders=borders, states=states, oceans=oceans, lakes=lakes, water_color='lightblue',
        lat_labels=lat_labels, lon_labels=lon_labels, suptitle_y=suptitle_y,
        lg_text=lg_text, lg_loc=lg_loc, lg_fontsize=lg_fontsize,
        mark1_lon=lon01_pm, mark1_lat=lat01_pm, mark1=markers[0], mark1_size=msize[0], mark1_color=colors[0],
        mark2_lon=lon02_pm, mark2_lat=lat02_pm, mark2=markers[1], mark2_size=msize[1], mark2_color=colors[1],
        mark3_lon=lon03_pm, mark3_lat=lat03_pm, mark3=markers[2], mark3_size=msize[2], mark3_color=colors[2],
        mark4_lon=lon04_pm, mark4_lat=lat04_pm, mark4=markers[3], mark4_size=msize[3], mark4_color=colors[3],
        mark5_lon=lon05_pm, mark5_lat=lat05_pm, mark5=markers[4], mark5_size=msize[4], mark5_color=colors[4],
        mark6_lon=lon06_pm, mark6_lat=lat06_pm, mark6=markers[5], mark6_size=msize[5], mark6_color=colors[5],
        mark7_lon=lon07_pm, mark7_lat=lat07_pm, mark7=markers[6], mark7_size=msize[6], mark7_color=colors[6],
        mark8_lon=lon08_pm, mark8_lat=lat08_pm, mark8=markers[7], mark8_size=msize[7], mark8_color=colors[7],
        mark9_lon=lon09_pm, mark9_lat=lat09_pm, mark9=markers[8], mark9_size=msize[8], mark9_color=colors[8],
        mark10_lon=lon10_pm, mark10_lat=lat10_pm, mark10=markers[9], mark10_size=msize[9], mark10_color=colors[9],
        mark11_lon=lon11_pm, mark11_lat=lat11_pm, mark11=markers[10],mark11_size=msize[10], mark11_color=colors[10],
        )
    fname_pm25.chmod(0o666)

def main():
    plot_type = 'png'

    plot_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','plots')

    mpl_o3   = '$\mathregular{O_3}$'
    mpl_pm25 = '$\mathregular{PM_{2.5}}$'

    ## Read in the EPA codes files that Ju-Hye made (these files have PM2.5 stations misordered/missorted)
    epa_codes_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI')
#    epa_codes_file_use = epa_codes_dir.joinpath('EPA_CODE_O3_PM25_USE.nc')
#    epa_codes_file_val = epa_codes_dir.joinpath('EPA_CODE_O3_PM25_VAL.nc')
    ## Instead read in the corrected EPA codes files that I made
    epa_codes_file_use = epa_codes_dir.joinpath('epa_code_o3_pm25_use_jared.nc')
    epa_codes_file_val = epa_codes_dir.joinpath('epa_code_o3_pm25_val_jared.nc')

    print('Reading '+str(epa_codes_file_use))
    epa_use_ds = xr.open_dataset(epa_codes_file_use)
#    epa_code_o3_use   = epa_use_ds.epa_code_o3.values.astype(int)
    epa_code_o3_use   = epa_use_ds.obs_o3_reg_use.values
    obs_o3_lat_use    = epa_use_ds.obs_o3_lat_use.values
    obs_o3_lon_use    = epa_use_ds.obs_o3_lon_use.values
#    epa_code_pm25_use = epa_use_ds.epa_code_pm25.values.astype(int)
    epa_code_pm25_use = epa_use_ds.obs_pm25_reg_use.values
    obs_pm25_lat_use  = epa_use_ds.obs_pm25_lat_use.values
    obs_pm25_lon_use  = epa_use_ds.obs_pm25_lon_use.values
    n_obs_o3_use   = len(epa_code_o3_use)
    n_obs_pm25_use = len(epa_code_pm25_use)

    print('Reading '+str(epa_codes_file_val))
    epa_val_ds = xr.open_dataset(epa_codes_file_val)
#    epa_code_o3_val   = epa_val_ds.epa_code_o3.values.astype(int)
    epa_code_o3_val   = epa_val_ds.obs_o3_reg_val.values
    obs_o3_lat_val    = epa_val_ds.obs_o3_lat_val.values
    obs_o3_lon_val    = epa_val_ds.obs_o3_lon_val.values
#    epa_code_pm25_val = epa_val_ds.epa_code_pm25.values.astype(int)
    epa_code_pm25_val = epa_val_ds.obs_pm25_reg_val.values
    obs_pm25_lat_val  = epa_val_ds.obs_pm25_lat_val.values
    obs_pm25_lon_val  = epa_val_ds.obs_pm25_lon_val.values
    n_obs_o3_val   = len(epa_code_o3_val)
    n_obs_pm25_val = len(epa_code_pm25_val)

#    for ss in range(n_obs_o3_use):
#        print('EPA Region '+str(epa_code_o3_use[ss])+', Lat = '+str(obs_o3_lat_use[ss])+' N, Lon = '+str(obs_o3_lon_use[ss])+' E')

    ## Combine the arrays
    epa_code_o3_all   = np.append(epa_code_o3_use, epa_code_o3_val)
    obs_o3_lat_all    = np.append(obs_o3_lat_use, obs_o3_lat_val)
    obs_o3_lon_all    = np.append(obs_o3_lon_use, obs_o3_lon_val)
    epa_code_pm25_all = np.append(epa_code_pm25_use, epa_code_pm25_val)
    obs_pm25_lat_all  = np.append(obs_pm25_lat_use, obs_pm25_lat_val)
    obs_pm25_lon_all  = np.append(obs_pm25_lon_use, obs_pm25_lon_val)
    n_obs_o3_all   = len(epa_code_o3_all)
    n_obs_pm25_all = len(epa_code_pm25_all)

    '''
    print('O3 obs')
    for ss in range(n_obs_o3_all):
        if epa_code_o3_all[ss] == 99:
            print('EPA Region '+str(epa_code_o3_all[ss])+', Lat = '+str(obs_o3_lat_all[ss])+' N, Lon = '+str(obs_o3_lon_all[ss])+' E')

    print('')
    print('PM2.5 obs')
    for ss in range(n_obs_pm25_all):
        if epa_code_pm25_all[ss] == 99:
            print('EPA Region '+str(epa_code_pm25_all[ss])+', Lat = '+str(obs_pm25_lat_all[ss])+' N, Lon = '+str(obs_pm25_lon_all[ss])+' E')
    '''

    suptitle_o3 = 'AirNow '+mpl_o3+' All Stations by EPA Region ('+str(n_obs_o3_all)+' Stations)'
    suptitle_pm25 = 'AirNow '+mpl_pm25+' All Stations by EPA Region ('+str(n_obs_pm25_all)+' Stations)'
    fname_o3   = plot_dir.joinpath('map_airnow_stn_reg_o3_all.'+plot_type)
    fname_pm25 = plot_dir.joinpath('map_airnow_stn_reg_pm25_all.'+plot_type)
    plot_stn_map(fname_o3, epa_code_o3_all, obs_o3_lat_all, obs_o3_lon_all, suptitle_o3,
        fname_pm25, epa_code_pm25_all, obs_pm25_lat_all, obs_pm25_lon_all, suptitle_pm25)

    suptitle_o3 = 'AirNow '+mpl_o3+' Training Stations by EPA Region ('+str(n_obs_o3_use)+' Stations)'
    suptitle_pm25 = 'AirNow '+mpl_pm25+' Training Stations by EPA Region ('+str(n_obs_pm25_use)+' Stations)'
    fname_o3   = plot_dir.joinpath('map_airnow_stn_reg_o3_use.'+plot_type)
    fname_pm25 = plot_dir.joinpath('map_airnow_stn_reg_pm25_use.'+plot_type)
    plot_stn_map(fname_o3, epa_code_o3_use, obs_o3_lat_use, obs_o3_lon_use, suptitle_o3,
        fname_pm25, epa_code_pm25_use, obs_pm25_lat_use, obs_pm25_lon_use, suptitle_pm25)

    suptitle_o3 = 'AirNow '+mpl_o3+' Validation Stations by EPA Region ('+str(n_obs_o3_val)+' Stations)'
    suptitle_pm25 = 'AirNow '+mpl_pm25+' Validation Stations by EPA Region ('+str(n_obs_pm25_val)+' Stations)'
    fname_o3   = plot_dir.joinpath('map_airnow_stn_reg_o3_val.'+plot_type)
    fname_pm25 = plot_dir.joinpath('map_airnow_stn_reg_pm25_val.'+plot_type)
    plot_stn_map(fname_o3, epa_code_o3_val, obs_o3_lat_val, obs_o3_lon_val, suptitle_o3,
        fname_pm25, epa_code_pm25_val, obs_pm25_lat_val, obs_pm25_lon_val, suptitle_pm25)

#    plot_stn_map(fname_o3, epa_code_o3_use, obs_o3_lat_use, obs_o3_lon_use,
#        fname_pm25, epa_code_pm25_use, obs_pm25_lat_use, obs_pm25_lon_use)
#    plot_stn_map(fname_o3, epa_code_o3_val, obs_o3_lat_val, obs_o3_lon_val,
#        fname_pm25, epa_code_pm25_val, obs_pm25_lat_val, obs_pm25_lon_val)

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
