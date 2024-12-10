'''
set_anen_box.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 26 Oct 2022

This script reads in the master list of observations (use/val), a CMAQ grid file, defines a small
subregion (box) of the CMAQ domain, and plots a map of the stations contained therein, along with
the CMAQ coordinates defining the subregion.
'''

import sys
import pathlib
import wrf
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from functions import map_functions

def main():
	obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static')
	out_dir = obs_dir

	date_beg = '20200801_00'
	date_end = '20211231_23'
	plot_type = 'png'

	i_beg = 20
	i_end = 60
	j_beg = 130
	j_end = 170

	fmt_yyyymmdd_hh = '%Y%m%d_%H'
	fmt_yyyymmddhh  = '%Y%m%d%H'
	fmt_yyyymmdd = '%Y%m%d'
	fmt_yyyy = '%Y'
	fmt_mm = '%m'
	fmt_hh = '%H'
	fmt_date = '%Y-%m-%d'
	fmt_date_hh = '%Y-%m-%d_%H'

	pm25_str = '$\mathregular{PM_{2.5}}$'
	o3_str = '$\mathregular{O_{3}}$'

	en_dash = u'\u2013'
	em_dash = u'\u2014'

	## Create datetime array
	dt_beg = dt.datetime.strptime(date_beg, fmt_yyyymmdd_hh)
	dt_end = dt.datetime.strptime(date_end, fmt_yyyymmdd_hh)
	dt_array = pd.date_range(start=dt_beg, end=dt_end, freq='1H')
	n_times = len(dt_array)

	## Open single files that contain all times
	all_file = obs_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_all.nc')
	use_file = obs_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_use.nc')
	val_file = obs_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_val.nc')

	print('Reading '+str(all_file))
	ds_all = xr.open_dataset(all_file)
	pm25_sid_all = ds_all.pm25_sid_all
	pm25_lon_all = ds_all.pm25_lon_all
	pm25_lat_all = ds_all.pm25_lat_all
	pm25_pct_all = ds_all.pm25_pct_all
	o3_sid_all = ds_all.o3_sid_all
	o3_lon_all = ds_all.o3_lon_all
	o3_lat_all = ds_all.o3_lat_all
	o3_pct_all = ds_all.o3_pct_all
	n_all_pm25 = len(pm25_sid_all)
	n_all_o3   = len(o3_sid_all)

	print('Reading '+str(use_file))
	ds_use = xr.open_dataset(use_file)
	pm25_sid_use = ds_use.pm25_sid_use
	pm25_lon_use = ds_use.pm25_lon_use
	pm25_lat_use = ds_use.pm25_lat_use
	pm25_pct_use = ds_use.pm25_pct_use
	o3_sid_use = ds_use.o3_sid_use
	o3_lon_use = ds_use.o3_lon_use
	o3_lat_use = ds_use.o3_lat_use
	o3_pct_use = ds_use.o3_pct_use
	n_use_pm25 = len(pm25_sid_use)
	n_use_o3   = len(o3_sid_use)

	print('Reading '+str(val_file))
	ds_val = xr.open_dataset(val_file)
	pm25_sid_val = ds_val.pm25_sid_val
	pm25_lon_val = ds_val.pm25_lon_val
	pm25_lat_val = ds_val.pm25_lat_val
	pm25_pct_val = ds_val.pm25_pct_val
	o3_sid_val = ds_val.o3_sid_val
	o3_lon_val = ds_val.o3_lon_val
	o3_lat_val = ds_val.o3_lat_val
	o3_pct_val = ds_val.o3_pct_val
	n_val_pm25 = len(pm25_sid_val)
	n_val_o3   = len(o3_sid_val)

	## ================
	## PLOTTING SECTION
	## ================

	suptitle_pm25_all = 'AirNow '+pm25_str+' Obs, All Sites'
	suptitle_o3_all   = 'AirNow '+o3_str+' Obs, All Sites'
	suptitle_pm25_use = 'AirNow '+pm25_str+' Obs, Sites Used for Training'
	suptitle_o3_use   = 'AirNow '+o3_str+' Obs, Sites Used for Training'
	suptitle_pm25_val = 'AirNow '+pm25_str+' Obs, Sites Held for Validation'
	suptitle_o3_val   = 'AirNow '+o3_str+' Obs, Sites Held for Validation'

	fontsize = 12
	suptitle_y = 1.00

	## Make maps with percentage valid reports for all stations within the period for the subdomain
	## First, open the CMAQ sample/coordinate file to get map projection parameters
	cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
	print('Reading CMAQ coordinate data from '+str(cmaq_fname))
	cmaq_ds = xr.open_dataset(cmaq_fname)
	cmaq_lon = cmaq_ds.LON[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
	cmaq_lat = cmaq_ds.LAT[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
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

#	print(cmaq_lon)
#	print(cmaq_lat)
#	sys.exit()

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
	data_crs = ccrs.PlateCarree()

	title = date_beg+' to '+date_end+'\nSubdomain CMAQ i: '+str(i_beg)+en_dash+str(i_end)+', CMAQ j: '+str(j_beg)+en_dash+str(j_end)
	cbar_lab = 'Percentage of Hours with Valid Reports'
	extend = 'max'
	cmap = mpl.cm.viridis
	bounds = np.arange(0.0, 95.1, 5.0)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)

	## Make maps showing all stations in the period
	suptitle = suptitle_pm25_all
	fname = out_dir.joinpath('map_sub_pm2.5_valid_obs_'+date_beg+'-'+date_end+'_all.'+plot_type)
	var = pm25_pct_all
	map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
		cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
		borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
		marker_lat=pm25_lat_all, marker_lon=pm25_lon_all, marker='o', marker_val_fill=True, marker_size=64,
		title_c=title, suptitle_y=suptitle_y)

	suptitle = suptitle_o3_all
	fname = out_dir.joinpath('map_sub_o3_valid_obs_'+date_beg+'-'+date_end+'_all.'+plot_type)
	var = o3_pct_all
	map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
		cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
		borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
		marker_lat=o3_lat_all, marker_lon=o3_lon_all, marker='o', marker_val_fill=True, marker_size=64,
		title_c=title, suptitle_y=suptitle_y)

	## Make maps showing use stations in the period
	suptitle = suptitle_pm25_use
	fname = out_dir.joinpath('map_sub_pm2.5_valid_obs_'+date_beg+'-'+date_end+'_use.'+plot_type)
	var = pm25_pct_use
	map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
		cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
		borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
		marker_lat=pm25_lat_use, marker_lon=pm25_lon_use, marker='o', marker_val_fill=True, marker_size=64,
		title_c=title, suptitle_y=suptitle_y)

	suptitle = suptitle_o3_use
	fname = out_dir.joinpath('map_sub_o3_valid_obs_'+date_beg+'-'+date_end+'_use.'+plot_type)
	var = o3_pct_use
	map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
		cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
		borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
		marker_lat=o3_lat_use, marker_lon=o3_lon_use, marker='o', marker_val_fill=True, marker_size=64,
		title_c=title, suptitle_y=suptitle_y)

	## Make maps showing val stations in the period
	suptitle = suptitle_pm25_val
	fname = out_dir.joinpath('map_sub_pm2.5_valid_obs_'+date_beg+'-'+date_end+'_val.'+plot_type)
	var = pm25_pct_val
	map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
		cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
		borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
		marker_lat=pm25_lat_val, marker_lon=pm25_lon_val, marker='o', marker_val_fill=True, marker_size=64,
		title_c=title, suptitle_y=suptitle_y)

	suptitle = suptitle_o3_val
	fname = out_dir.joinpath('map_sub_o3_valid_obs_'+date_beg+'-'+date_end+'_val.'+plot_type)
	var = o3_pct_val
	map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
		cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
		borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
		marker_lat=o3_lat_val, marker_lon=o3_lon_val, marker='o', marker_val_fill=True, marker_size=64,
		title_c=title, suptitle_y=suptitle_y)

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
