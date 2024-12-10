#!/usr/bin/env python
'''
plot_cmaq_ncar_bc_pm25_o3.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 17 Feb 2023

This script will plot maps of requested CMAQ forecast cycles/lead times fields of PM2.5 and O3.
The CMAQ output is pulled from the NCAR AnEn output, and plots the AnEn mean of 10 members.
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
	parser.add_argument('cycle_beg', help='cycle date/time of first CMAQ simulation [YYYYMMDD_HH]')
	parser.add_argument('-l', '--cycle_end', default=None, help='cycle date/time of last CMAQ simulation [YYYYMMDD_HH]')
	parser.add_argument('-b', '--lead_beg', default=0, type=int, help='beginning lead hour for plotting (default: 0)')
	parser.add_argument('-e', '--lead_end', default=72, type=int, help='ending lead hour for plotting (default: 72)')
	parser.add_argument('-s', '--subdomain', action='store_true', help='optional flag to plot the California subdomain instead of CONUS (default: False)')

	args = parser.parse_args()
	cycle_beg = args.cycle_beg
	cycle_end = args.cycle_end
	lead_beg  = args.lead_beg
	lead_end  = args.lead_end

	if len(cycle_beg) != 11:
		print('ERROR! Incorrect length for positional argument cycle_beg. Exiting!')
		parser.print_help()
		sys.exit()
	elif cycle_beg[8] != '_':
		print('ERROR! Incorrect format for positional argument cycle_beg. Exiting!')

	if cycle_end != None:
		if len(cycle_end) != 11:
			print('ERROR! Incorrect length for optional argument cycle_end. Exiting!')
			parser.print_help()
			sys.exit()
		elif cycle_end[8] != '_':
			print('ERROR! Incorrect format for optional argument cycle_end. Exiting!')
	else:
		cycle_end = cycle_beg

	plot_subdomain = False
	if args.subdomain:
		plot_subdomain = True

	return cycle_beg, cycle_end, lead_beg, lead_end, plot_subdomain

def main(cycle_beg, cycle_end, lead_beg, lead_end, plot_subdomain):

	plot_type = 'png'
#	plot_subdomain = True
	plot_o3 = True
	plot_pm25 = True
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

#	cmaq_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bc_cmaq_06z_using_anen')
	cmaq_dir = pathlib.Path('/','glade','u','home','alessand')
	out_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cmaq','ncar_bc','maps')
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

	## Create datetime array of model cycles
	dt_beg_cycle = dt.datetime.strptime(cycle_beg, fmt_yyyymmdd_hh)
	dt_end_cycle = dt.datetime.strptime(cycle_end, fmt_yyyymmdd_hh)
	dt_arr_cycle = pd.date_range(start=dt_beg_cycle, end=dt_end_cycle, freq='24H')
	n_cycle = len(dt_arr_cycle)

	## Create array of lead times
	leads = np.arange(lead_beg, lead_end+1)
	n_leads = len(leads)

	## Open a CMAQ sample/coordinate file to get map projection parameters
	cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
	print('Reading CMAQ coordinate data from '+str(cmaq_fname))
	cmaq_ds = xr.open_dataset(cmaq_fname)
	if plot_subdomain:
		cmaq_lon = cmaq_ds.LON[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
		cmaq_lat = cmaq_ds.LAT[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
	else:
		cmaq_lon = cmaq_ds.LON[0,0,j_beg:,i_beg:].rename({'ROW':'latitude', 'COL':'longitude'})
		cmaq_lat = cmaq_ds.LAT[0,0,j_beg:,i_beg:].rename({'ROW':'latitude', 'COL':'longitude'})
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

	## Read in AnEn output files
	## Stefano has one file for PM2.5, one for O3 on the California subdomain
	## Both AnEn output files cover all cycles (starting with 1 Jan 2021) & lead times
	anen_out_fname_o3   = cmaq_dir.joinpath('anen_output_o3.nc')
	anen_out_fname_pm25 = cmaq_dir.joinpath('anen_output_PM25.nc')

	print('Reading '+str(anen_out_fname_o3))
	ds_cmaq_anen_o3 = xr.open_dataset(anen_out_fname_o3)
	n_anen_days = len(ds_cmaq_anen_o3.NDAYS)
	n_anen_leads = len(ds_cmaq_anen_o3.NLEADTIMES)
	n_anen_stns = len(ds_cmaq_anen_o3.NSTATIONS)		# for the California subdomain, these are the model grid points
	n_anen_members = len(ds_cmaq_anen_o3.NMEMBERS)
	anen_o3_all = ds_cmaq_anen_o3.AnEn	# [NMEMBERS, NLEADTIMES, NDAYS, NSTATIONS]
	
	print('Reading '+str(anen_out_fname_pm25))
	ds_cmaq_anen_pm25 = xr.open_dataset(anen_out_fname_pm25)
	anen_pm25_all = ds_cmaq_anen_pm25.AnEn	# [NMEMBERS, NLEADTIMES, NDAYS, NSTATIONS]

	## Also read in a file with the lat/lon coordinates of all 1600 points in the 1D array in the AnEn output files
	anen_latlon_fname = cmaq_dir.joinpath('InputAnEn.OBS.o3_pm25.1600st.72hr.20200801-20211231.Cali.o3_bm3.pm25_bm1.latlon.include0.nc')
	print('Reading '+str(anen_latlon_fname))
	anen_latlon_ds = xr.open_dataset(anen_latlon_fname)
	anen_lat = anen_latlon_ds.LAT.values
	anen_lon = anen_latlon_ds.LON.values

	## When reshaping these (1600,) 1D arrays into (40,40) 2D arrays, do they match up with the cmaq lat/lon?
	## No, they don't with order='C' (default), which means np.reshape is mismatching the AnEn output geographically.
	## Need to either use order='F' (Fortran-like) or do a transpose of the matrix to get it to match cmaq lat/lon.
	anen_lat_2d = np.reshape(anen_lat, (40,40), order='F')
	anen_lon_2d = np.reshape(anen_lon, (40,40), order='F')
#	print(anen_lat_2d[0,:])
#	print(np.transpose(anen_lat_2d)[0,:])
#	print(cmaq_lat.values[0,:])
#	print('')
#	print(anen_lon_2d)
#	print(cmaq_lon.values)
	## Find closest j,i index of a particular station (use/train pm2.5 station ss = 427, lat = 40.3, lon = -121.24)
	dist = abs(cmaq_lat-40.3)+abs(cmaq_lon-(-121.24))
	cmaq_j,cmaq_i = np.unravel_index(dist.argmin(), dist.shape)
#	print('cmaq_j = '+str(cmaq_j)+', cmaq_i = '+str(cmaq_i))
	## For completeness, find the index value of this same point in the 1D (1600,) array
	ind1d = np.where(anen_lon == cmaq_lon.values[cmaq_j,cmaq_i])[0][0]
#	print(anen_lat[ind1d])
#	print(anen_lon[ind1d])
#	print(cmaq_lat.values[cmaq_j,cmaq_i])
#	print(cmaq_lon.values[cmaq_j,cmaq_i])
#	print(anen_lat_2d[cmaq_j,cmaq_i])
#	print(anen_lon_2d[cmaq_j,cmaq_i])
#	print(np.transpose(anen_lat_2d)[cmaq_j,cmaq_i])
#	print(np.transpose(anen_lon_2d)[cmaq_j,cmaq_i])
#	print(ind1d)
#	sys.exit()

	## Loop over CMAQ forecast cycles
	for cc in range(n_cycle):
		this_cycle_dt = dt_arr_cycle[cc]
		this_cycle_infile = this_cycle_dt.strftime(fmt_yyyymmdd+'.'+fmt_hh+'z')
		this_cycle_file = this_cycle_dt.strftime(fmt_yyyymmdd_hhmm)
		this_cycle_dir = this_cycle_dt.strftime(fmt_yyyymmdd_hh)
		this_cycle_date = this_cycle_dt.strftime(fmt_yyyymmdd)
		this_cycle_yr = this_cycle_dt.strftime(fmt_yyyy)
		this_cycle_mo = this_cycle_dt.strftime(fmt_mm)
		this_cycle_plot = this_cycle_dt.strftime(fmt_date_plot)

		## To use Stefano's AnEn output files, calculate no. of days since first cycle (known to be 1 Jan 2021)
		anen_out_start_str = '20210101_0600'
		time_since_start = this_cycle_dt - dt.datetime.strptime(anen_out_start_str, fmt_yyyymmdd_hhmm)
		days_since_start = time_since_start.days
#		print(anen_o3_all.shape)
#		print('days_since_start = '+str(days_since_start))
		anen_o3_this_cycle = anen_o3_all[:,:,days_since_start,:]
		anen_pm25_this_cycle = anen_pm25_all[:,:,days_since_start,:]

		## Create the output directory if needed, ensuring it has necessary permissions
		out_dir = out_dir_parent.joinpath(this_cycle_dir)
		out_dir.mkdir(exist_ok=True, parents=True)
		out_dir.chmod(0o777)

		suptitle_o3 = 'CMAQ Forecast (NCAR BC AnEn mean)'
		suptitle_pm25 = 'CMAQ Forecast (NCAR BC AnEn mean)'
		if plot_airnow_vals:
			suptitle_o3 = suptitle_o3+' & AirNow Obs'
			suptitle_pm25 = suptitle_pm25+' & AirNow Obs'
		elif plot_airnow_stns:
			suptitle_o3 = suptitle_o3+' & AirNow Sites'
			suptitle_pm25 = suptitle_pm25+' & AirNow Sites'

		## Loop over the lead times
		for ll in leads:
			this_valid_dt = this_cycle_dt + dt.timedelta(hours=int(ll))
			this_valid_file = this_valid_dt.strftime(fmt_date_file)
			this_valid_plot = this_valid_dt.strftime(fmt_date_plot)
			this_lead_str = f'{ll:02d}'

			title_r = 'Start: '+this_cycle_plot+' UTC\nValid: '+this_valid_plot+' UTC'

#			print(anen_o3_this_cycle.shape)
#			print('lead = '+str(ll))

			## Pull out the array slices for this lead time
			anen_o3_this_lead = anen_o3_this_cycle[:,ll,:]
			anen_pm25_this_lead = anen_pm25_this_cycle[:,ll,:]

#			print(anen_o3_this_lead.shape)

			'''
			## Calculate AnEn mean at each station/grid point
			anen_mean_o3_1d = np.nanmean(anen_o3_this_lead, axis=0)
			anen_mean_pm25_1d = np.nanmean(anen_pm25_this_lead, axis=0)

			print(anen_pm25_this_lead.values[:,0])
			print(anen_mean_pm25_1d[0])
			print(anen_mean_o3_1d.shape)

			## Reshape the array from 1D (1600 pts) to 2D (40x40)
			anen_mean_o3_2d = np.reshape(anen_mean_o3_1d, (40,40), order='F')
			anen_mean_pm25_2d = np.reshape(anen_mean_pm25_1d, (40,40), order='F')
			'''

			## Reshape the (10,1600) array to (10,40,40)
			## As above, need to use order='F' or do a transpose to maintain agreement with cmaq lat/lon arrays
			anen_o3_this_lead_2d   = np.reshape(anen_o3_this_lead.values, (10,40,40), order='F')
			anen_pm25_this_lead_2d = np.reshape(anen_pm25_this_lead.values, (10,40,40), order='F')

#			print(anen_pm25_this_lead.values[:,ind1d])
#			print(anen_pm25_this_lead_2d[:,cmaq_j,cmaq_i])

			## Calculate the AnEn mean at each "station"/grid point
			anen_mean_o3_2d   = np.nanmean(anen_o3_this_lead_2d, axis=0)
			anen_mean_pm25_2d = np.nanmean(anen_pm25_this_lead_2d, axis=0)

#			print(anen_mean_pm25_2d.shape)
#			print(anen_mean_pm25_2d[cmaq_j,cmaq_i])

			## Do we want to add AirNow station locations or values to the maps?
			if plot_airnow_stns or plot_airnow_vals:
				obs_dir = obs_dir_parent.joinpath(this_cycle_yr, this_cycle_mo)
				use_file = obs_dir.joinpath('airnow_pm2.5_o3_'+this_valid_file+'_use.nc')
				val_file = obs_dir.joinpath('airnow_pm2.5_o3_'+this_valid_file+'_val.nc')

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
#				for ss in range(len(pm25_sid_use)):
#					if pm25_lat_use[ss] > 40.0 and pm25_lat_use[ss] < 41.0 and pm25_lon_use[ss] < -121.0 and pm25_lon_use[ss] > -122.0:
#						print('ss = '+str(ss)+', lat = '+str(pm25_lat_use.values[ss])+', lon = '+str(pm25_lon_use.values[ss]))
#				sys.exit()

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
					fname = out_dir.joinpath('map_calif_cmaq_ncar_bc_o3_'+this_valid_file+'.'+plot_type)
					lg_loc = 'upper left'
					marksize = 81
				else:
					fname = out_dir.joinpath('map_conus_cmaq_ncar_bc_o3_'+this_valid_file+'.'+plot_type)
					lg_loc = 'lower left'
					marksize = 49

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
						mark_fill = False
						mark1_var = None
						mark2_var = None
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
					mark_fill = False
					lg_text = None

#				var = cmaq_o3[ll,:,:]
				var = anen_mean_o3_2d
				var_min = np.nanmin(var)
				var_max = np.nanmax(var)
				var_name = 'Surface-level '+mpl_o3+' Conc.'
				var_cbar = 'Ozone Concentration'
				var_unit = 'ppbv'
				cbar_lab = var_cbar+' ['+var_unit+']'
				extend = 'max'
#				cmap = mpl.cm.rainbow
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
					fname = out_dir.joinpath('map_calif_cmaq_ncar_bc_pm25_'+this_valid_file+'.'+plot_type)
					lg_loc = 'lower right'
					marksize = 81
				else:
					fname = out_dir.joinpath('map_conus_cmaq_ncar_bc_pm25_'+this_valid_file+'.'+plot_type)
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

#				var = cmaq_pm25[ll,:,:]
				var = anen_mean_pm25_2d
				var_min = np.nanmin(var)
				var_max = np.nanmax(var)
				var_name = 'Surface-level '+mpl_pm25+' Conc.'
				var_cbar = mpl_pm25+' Concentration'
				var_unit = mpl_ugm3
				cbar_lab = var_cbar+' ['+var_unit+']'
				extend = 'max'
#				cmap = mpl.cm.rainbow
				cmap = mpl.cm.get_cmap('rainbow').copy()
				cmap.set_bad('white')
#				bounds = np.arange(min_pm25, max_pm25, int_pm25)
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
	cycle_beg, cycle_end, lead_beg, lead_end, plot_subdomain = parse_args()
	main(cycle_beg, cycle_end, lead_beg, lead_end, plot_subdomain)
	now_time_end = dt.datetime.utcnow()
	run_time_tot = now_time_end - now_time_beg
	now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
	now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
	print('\nScript completed successfully.')
	print('   Beg time: '+now_time_beg_str)
	print('   End time: '+now_time_end_str)
	print('   Run time: '+str(run_time_tot)+'\n')
