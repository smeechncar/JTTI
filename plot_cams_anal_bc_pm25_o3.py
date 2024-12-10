#!/usr/bin/env python
'''
plot_cams_anal_bc_pm25_o3.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 20 Dec 2022

This script will plot maps of requested CAMS bias-corrected analysis fields of PM2.5 and O3.
The CAMS bias-corrected files that are read in and plotted are prepared by Ju-Hye, and have
been regridded to the CMAQ curvilinear 12-km grid.
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

def parse_args():
	## Parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('beg_dt_str', help='beginning time for plotting [YYYYMMDD_HH]')
	parser.add_argument('-e', '--end_dt_str', default=None, help='ending time for plotting [YYYYMMDD_HH]')
	parser.add_argument('-s', '--subdomain', action='store_true', help='optional flag to plot the California subdomain instead of CONUS (default: False)')

	args = parser.parse_args()
	beg_dt_str = args.beg_dt_str
	end_dt_str = args.end_dt_str

	if len(beg_dt_str) != 11:
		print('ERROR! Incorrect length for positional argument beg_dt_str. Exiting!')
		parser.print_help()
		sys.exit()
	elif beg_dt_str[8] != '_':
		print('ERROR! Incorrect format for positional argument beg_dt_str. Exiting!')

	if end_dt_str is not None:
		if len(end_dt_str) != 11:
			print('ERROR! Incorrect length for optional argument end_dt_str. Exiting!')
			parser.print_help()
			sys.exit()
		elif end_dt_str[8] != '_':
			print('ERROR! Incorrect format for optional argument end_dt_str. Exiting!')
	else:
		end_dt_str = beg_dt_str

	plot_subdomain = False
	if args.subdomain:
		plot_subdomain = True

	return beg_dt_str, end_dt_str, plot_subdomain

def main(beg_dt_str, end_dt_str, plot_subdomain):

	plot_type = 'png'
#	plot_subdomain = True
	plot_airnow_stns = False
	plot_airnow_vals = True

	if plot_subdomain:
		i_beg, i_end = 20, 60
		j_beg, j_end = 130, 170
		lat_labels = [36, 37, 38, 39, 40, 41]
		lon_labels = [-124, -123, -122, -121, -120, -119, -118]
		suptitle_y = 1.01
	else:
		i_beg, i_end = 0, -1
		j_beg, j_end = 0, -1
		lat_labels = [25, 30, 35, 40, 45, 50]
		lon_labels = [-130, -120, -110, -110, -100, -90, -80, -70, -60]
		suptitle_y = 0.91

	cams_bm1_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','BM1','static','cams_include0')
	cams_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','BM3','static','cams_include0')
	out_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams_bc')
	obs_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static','20200801_00-20211231_23')

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

	mpl_o3   = '$\mathregular{O_3}$'
	mpl_pm25 = '$\mathregular{PM_{2.5}}$'
	mpl_ugm3 = '\u03bcg $\mathregular{m^{-3}}$'

	min_o3 = 0.0
	max_o3 = 90.1
	int_o3 = 5.0

	min_pm25 = 0.0
	max_pm25 = 50.1
	int_pm25 = 5.0
	cont_pm25 = np.arange(min_pm25, max_pm25, int_pm25)
	cont_pm25 = np.append(cont_pm25, np.arange(60.0, 100.1, 10.0))
	cont_pm25 = np.append(cont_pm25, np.arange(125.0, 200.1, 25.0))

	## Create datetime array of times to read and plot
	beg_dt = dt.datetime.strptime(beg_dt_str, fmt_yyyymmdd_hh)
	end_dt = dt.datetime.strptime(end_dt_str, fmt_yyyymmdd_hh)
	all_dt = pd.date_range(start=beg_dt, end=end_dt, freq='1H')
	n_times = len(all_dt)

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

	## Loop over times
	for tt in range(n_times):
		this_dt = all_dt[tt]
		this_dt_yr = this_dt.strftime(fmt_yyyy)
		this_dt_mo = this_dt.strftime(fmt_mm)
		this_dt_date = this_dt.strftime(fmt_yyyymmdd)
		this_dt_file = this_dt.strftime(fmt_date_file)
		this_dt_plot = this_dt.strftime(fmt_date_plot)

		## Read in AnEn bias-corrected CAMS files (Ju-Hye has all CAMS BM[1-3] hourly files in one directory)
		cams_fname_pm25 = cams_bm1_dir.joinpath('cams_pm25_regrid_cmaq_'+this_dt_file+'_BM1_static_include0.nc')
		cams_fname_o3   = cams_bm3_dir.joinpath('cams_o3_regrid_cmaq_'+this_dt_file+'_BM3_static_include0.nc')

		if cams_fname_pm25.is_file():
			print('Reading '+str(cams_fname_pm25))
			ds_cams_pm25 = xr.open_dataset(cams_fname_pm25)	# Dimensions (latitude 265, longitude 442)
			if plot_subdomain:
				cams_pm25 = ds_cams_pm25.cams_pm25_m[j_beg:j_end,i_beg:i_end]
			else:
				cams_pm25 = ds_cams_pm25.cams_pm25_m[:,:]
			plot_pm25 = True
		else:
			plot_pm25 = False

		if cams_fname_o3.is_file():
			print('Reading '+str(cams_fname_o3))
			ds_cams_o3 = xr.open_dataset(cams_fname_o3)
			if plot_subdomain:
				cams_o3 = ds_cams_o3.cams_o3_m[j_beg:j_end,i_beg:i_end]
			else:
				cams_o3 = ds_cams_o3.cams_o3_m[:,:]
			plot_o3 = True
		else:
			plot_o3 = False

		## Create the output directory if needed, ensuring it has necessary permissions
		out_dir = out_dir_parent.joinpath(this_dt_yr,this_dt_mo)
		out_dir.mkdir(exist_ok=True, parents=True)
		out_dir.chmod(0o777)

		suptitle_o3 = 'CAMS BC Analysis (BM3 Method)'
		suptitle_pm25 = 'CAMS BC Analysis (BM1 Method)'
		if plot_airnow_vals:
			suptitle_o3 = suptitle_o3+' & AirNow Obs'
			suptitle_pm25 = suptitle_pm25+' & AirNow Obs'
		elif plot_airnow_stns:
			suptitle_o3 = suptitle_o3+' & AirNow Sites'
			suptitle_pm25 = suptitle_pm25+' & AirNow Sites'

		title_r = 'Valid: '+this_dt_plot+' UTC'

		## Do we want to add AirNow station locations or values to the maps?
		if plot_airnow_stns or plot_airnow_vals:
			obs_dir = obs_dir_parent.joinpath(this_dt_yr, this_dt_mo)
			use_file = obs_dir.joinpath('airnow_pm2.5_o3_'+this_dt_file+'_use.nc')
			val_file = obs_dir.joinpath('airnow_pm2.5_o3_'+this_dt_file+'_val.nc')

			print('Reading '+str(use_file))
			ds_use = xr.open_dataset(use_file)
			pm25_sid_use = ds_use.pm25_sid_use
			pm25_lon_use = ds_use.pm25_lon_use
			pm25_lat_use = ds_use.pm25_lat_use
			pm25_con_use = ds_use.pm25_con_use
			o3_sid_use = ds_use.o3_sid_use
			o3_lon_use = ds_use.o3_lon_use
			o3_lat_use = ds_use.o3_lat_use
			o3_con_use = ds_use.o3_con_use
			n_use_pm25 = len(pm25_sid_use)
			n_use_o3   = len(o3_sid_use)

			print('Reading '+str(val_file))
			ds_val = xr.open_dataset(val_file)
			pm25_sid_val = ds_val.pm25_sid_val
			pm25_lon_val = ds_val.pm25_lon_val
			pm25_lat_val = ds_val.pm25_lat_val
			pm25_con_val = ds_val.pm25_con_val
			o3_sid_val = ds_val.o3_sid_val
			o3_lon_val = ds_val.o3_lon_val
			o3_lat_val = ds_val.o3_lat_val
			o3_con_val = ds_val.o3_con_val
			n_val_pm25 = len(pm25_sid_val)
			n_val_o3   = len(o3_sid_val)

		if plot_o3:
			if plot_subdomain:
				fname = out_dir.joinpath('map_calif_cams_bc_bm3_o3_'+this_dt_file+'.'+plot_type)
				lg_loc = 'upper left'
				marksize = 81
			else:
				fname = out_dir.joinpath('map_conus_cams_bc_bm3_o3_'+this_dt_file+'.'+plot_type)
				lg_loc = 'lower left'
				marksize = 49

			mark_fill = False
			if plot_airnow_stns or plot_airnow_vals:	
				mark1_lon = o3_lon_use
				mark1_lat = o3_lat_use
				mark2_lon = o3_lon_val
				mark2_lat = o3_lat_val
				mark1_col = 'None'
				mark2_col = 'None'
				mark1 = 'v'
				mark2 = 'o'
				lg_text = ['AirNow Training', 'AirNow Validation']
				if plot_airnow_vals:
					mark_fill = True
					mark1_var = o3_con_use
					mark2_var = o3_con_val
			else:
				mark1_lon = None
				mark1_lat = None
				mark2_lon = None
				mark2_lat = None
				mark1_var = None
				mark2_var = None
				mark1_col = None
				mark2_col = None
				mark1 = None
				mark2 = None
				lg_text = None

			var = cams_o3[:,:]
			var_min = np.nanmin(var)
			var_max = np.nanmax(var)
			var_name = 'Surface-level '+mpl_o3+' Conc.'
			var_cbar = 'Ozone Concentration'
			var_unit = 'ppbv'
			cbar_lab = var_cbar+' ['+var_unit+']'
			extend = 'max'
#			cmap = mpl.cm.rainbow
			cmap = mpl.cm.get_cmap('rainbow').copy()
			cmap.set_bad('white')
			bounds = np.arange(min_o3, max_o3, int_o3)
			norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)
			title_l = var_name+f'\nMin: {var_min:.1f} '+var_unit+f', Max: {var_max:.1f} '+var_unit
			map_functions.mpl_map_plot(fname, var, suptitle_o3, cbar_lab, cart_proj, cart_xlim, cart_ylim,
				cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
				borders=borders, states=states, lakes=lakes, water_color='none',
				title_l=title_l, title_r=title_r, suptitle_y=suptitle_y,
				marker_lon=mark1_lon, marker_lat=mark1_lat, marker_var=mark1_var, marker_val_fill=mark_fill,
				marker=mark1, marker_color=mark1_col, marker_size=marksize,
				marker_lon2=mark2_lon, marker_lat2=mark2_lat, marker_var2=mark2_var,
				marker2=mark2, marker_color2=mark2_col, marker_size2=marksize, marker_zorder=12,
				lg_text=lg_text, lg_loc=lg_loc)
			fname.chmod(0o666)

		if plot_pm25:
			if plot_subdomain:
				fname = out_dir.joinpath('map_calif_cams_bc_bm1_pm25_'+this_dt_file+'.'+plot_type)
				lg_loc = 'lower right'
				marksize = 81
			else:
				fname = out_dir.joinpath('map_conus_cams_bc_bm1_pm25_'+this_dt_file+'.'+plot_type)
				lg_loc = 'lower left'
				marksize = 49

			if plot_airnow_stns or plot_airnow_vals:
				mark1_lon = pm25_lon_use
				mark1_lat = pm25_lat_use
				mark2_lon = pm25_lon_val
				mark2_lat = pm25_lat_val
				mark1_col = 'None'
				mark2_col = 'None'
				mark1 = 'v'
				mark2 = 'o'
				lg_text = ['AirNow Training', 'AirNow Validation']
				if plot_airnow_vals:
					mark1_var = pm25_con_use
					mark2_var = pm25_con_val
			else:
				mark1_lon = None
				mark1_lat = None
				mark2_lon = None
				mark2_lat = None
				mark1_var = None
				mark2_var = None
				mark1_col = None
				mark2_col = None
				mark1 = None
				mark2 = None
				lg_text = None

			var = cams_pm25[:,:]
			var_min = np.nanmin(var)
			var_max = np.nanmax(var)
			var_name = 'Surface-level '+mpl_pm25+' Conc.'
			var_cbar = mpl_pm25+' Concentration'
			var_unit = mpl_ugm3
			cbar_lab = var_cbar+' ['+var_unit+']'
			extend = 'max'
#			cmap = mpl.cm.rainbow
			cmap = mpl.cm.get_cmap('rainbow').copy()
			cmap.set_bad('white')
#			bounds = np.arange(min_pm25, max_pm25, int_pm25)
			bounds = cont_pm25
			norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)
			title_l = var_name+f'\nMin: {var_min:.1f} '+var_unit+f', Max: {var_max:.1f} '+var_unit
			map_functions.mpl_map_plot(fname, var, suptitle_pm25, cbar_lab, cart_proj, cart_xlim, cart_ylim,
				cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
				borders=borders, states=states, lakes=lakes, water_color='none',
				title_l=title_l, title_r=title_r, suptitle_y=suptitle_y,
				marker_lon=mark1_lon, marker_lat=mark1_lat, marker_var=mark1_var, marker_val_fill=mark_fill,
				marker=mark1, marker_color=mark1_col, marker_size=marksize,
				marker_lon2=mark2_lon, marker_lat2=mark2_lat, marker_var2=mark2_var,
				marker2=mark2, marker_color2=mark2_col, marker_size2=marksize, marker_zorder=12,
				lg_text=lg_text, lg_loc=lg_loc)
			fname.chmod(0o666)


if __name__ == '__main__':
	now_time_beg = dt.datetime.utcnow()
	beg_dt_str, end_dt_str, plot_subdomain = parse_args()
	main(beg_dt_str, end_dt_str, plot_subdomain)
	now_time_end = dt.datetime.utcnow()
	run_time_tot = now_time_end - now_time_beg
	now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
	now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
	print('\nScript completed successfully.')
	print('   Beg time: '+now_time_beg_str)
	print('   End time: '+now_time_end_str)
	print('   Run time: '+str(run_time_tot)+'\n')
