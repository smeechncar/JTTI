'''
extract_anen_cmaq_bm_at_airnow_use_val_leadtime.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 14 Dec 2022

This script first reads in AirNow obs NetCDF files that have already been separated into use & val
sets, for use in merging & validation, respectively. Second, it reads in AnEn-corrected CMAQ
forecast files (one file per cycle with all 73 lead times for that cycle) to
get the grid coordinates of those stations and spatially interpolate the gridded values to the
AirNow sites via bilinear interpolation. CMAQ PM2.5 & O3 files are produced hourly. The output
files with CMAQ-interpolated values at AirNow sites are in NetCDF.
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

def parse_args():
	## Parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('date_beg', help='beginning CMAQ cycle date/time to be processed [YYYYMMDD_HH]')
	parser.add_argument('-l', '--date_end', default=None, help='ending CMAQ cycle date/time to be processed [YYYYMMDD_HH]')
	parser.add_argument('-f', '--lead_hrs', default=72, help='forecast lead hours [default: 72, for hours 0-72]')

	args = parser.parse_args()
	date_beg = args.date_beg
	date_end = args.date_end
	lead_hrs = args.lead_hrs

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

	return date_beg, date_end, lead_hrs


def main(date_beg, date_end, lead_hrs):
	sites_vary = False
	sites_static = True

	sub_sites = True	# process sites over a subdomain?

	## Set directories
	if sites_vary:
		obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_vary')
		out_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cmaq','anen_bc','airnow','sites_vary')
		sub_sid_dir = None
	if sites_static:
		obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static','20200801_00-20211231_23')
		out_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cmaq','anen_bc','airnow','sites_static')
		sub_sid_dir = pathlib.Path('/','glade','p','ral','wsap','jaredlee','NOAA_CMAQ_AnEn','calif')
#	cmaq_dir = pathlib.Path('/','glade','scratch','jkim','CMAQ','cmaq_06z_prep_anen')
	cmaq_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bc_cmaq_06z_using_anen')

	fmt_yyyymmdd_hhmm = '%Y%m%d_%H%M'
	fmt_yyyymmdd_hh = '%Y%m%d_%H'
	fmt_yyyymmddhh  = '%Y%m%d%H'
	fmt_yyyymmdd = '%Y%m%d'
	fmt_yyyy = '%Y'
	fmt_mm = '%m'
	fmt_hh = '%H'
	fmt_date = '%Y-%m-%d'
	fmt_date_hh = '%Y-%m-%d_%H'

	## Create datetime array of model cycles
	dt_beg_cycle = dt.datetime.strptime(date_beg, fmt_yyyymmdd_hh)
	dt_end_cycle = dt.datetime.strptime(date_end, fmt_yyyymmdd_hh)
	dt_arr_cycle = pd.date_range(start=dt_beg_cycle, end=dt_end_cycle, freq='24H')
	n_cycle = len(dt_arr_cycle)

	## Create datetime array of valid times
	dt_beg_valid = dt.datetime.strptime(date_beg, fmt_yyyymmdd_hh)
	dt_end_valid = dt.datetime.strptime(date_end, fmt_yyyymmdd_hh)+dt.timedelta(hours=lead_hrs)
	dt_arr_valid = pd.date_range(start=dt_beg_valid, end=dt_end_valid, freq='1H')
	n_valid = len(dt_arr_valid)

	## Create 2D array that holds all valid times for each model cycle
	dt_arr_pairs = np.full([n_cycle, lead_hrs+1], None)
	for cc in range(n_cycle):
		this_cycle = dt_arr_cycle[cc]
		for vv in range(lead_hrs+1):
			this_valid = this_cycle + dt.timedelta(hours=vv)
			dt_arr_pairs[cc,vv] = this_valid

	get_cmaq_lat_lon = True
	get_cmaq_proj = True

	## Loop over valid times to get the observations as a time series
	for vv in range(n_valid):
		this_dt = dt_arr_valid[vv]
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
		obs_ds_use = obs_ds_use.reset_coords(names='n_obs_ise_o3', drop=True)	# Drop unwanted typo coordinate

		if vv == 0:
			obs_pm25_sid_use = obs_ds_use.pm25_sid_use
			obs_pm25_lon_use = obs_ds_use.pm25_lon_use
			obs_pm25_lat_use = obs_ds_use.pm25_lat_use
			obs_pm25_cmaq_i_use = obs_ds_use.pm25_cmaq_i_use
			obs_pm25_cmaq_j_use = obs_ds_use.pm25_cmaq_j_use

			obs_o3_sid_use = obs_ds_use.o3_sid_use
			obs_o3_lon_use = obs_ds_use.o3_lon_use
			obs_o3_lat_use = obs_ds_use.o3_lat_use
			obs_o3_cmaq_i_use = obs_ds_use.o3_cmaq_i_use
			obs_o3_cmaq_j_use = obs_ds_use.o3_cmaq_j_use

			n_obs_pm25_use = obs_pm25_sid_use.shape[0]
			n_obs_o3_use = obs_o3_sid_use.shape[0]

			obs_pm25_con_use_ts = np.full([n_obs_pm25_use, n_valid], np.nan)
			obs_o3_con_use_ts = np.full([n_obs_o3_use, n_valid], np.nan)

		obs_pm25_con_use_ts[:,vv] = obs_ds_use.pm25_con_use
		obs_o3_con_use_ts[:,vv] = obs_ds_use.o3_con_use

		print('Reading '+str(obs_file_val))
		obs_ds_val = xr.open_dataset(obs_file_val)

		if vv == 0:
			obs_pm25_sid_val = obs_ds_val.pm25_sid_val
			obs_pm25_lon_val = obs_ds_val.pm25_lon_val
			obs_pm25_lat_val = obs_ds_val.pm25_lat_val
			obs_pm25_cmaq_i_val = obs_ds_val.pm25_cmaq_i_val
			obs_pm25_cmaq_j_val = obs_ds_val.pm25_cmaq_j_val

			obs_o3_sid_val = obs_ds_val.o3_sid_val
			obs_o3_lon_val = obs_ds_val.o3_lon_val
			obs_o3_lat_val = obs_ds_val.o3_lat_val
			obs_o3_cmaq_i_val = obs_ds_val.o3_cmaq_i_val
			obs_o3_cmaq_j_val = obs_ds_val.o3_cmaq_j_val

			n_obs_pm25_val = obs_pm25_sid_val.shape[0]
			n_obs_o3_val = obs_o3_sid_val.shape[0]

			obs_pm25_con_val_ts = np.full([n_obs_pm25_val, n_valid], np.nan)
			obs_o3_con_val_ts = np.full([n_obs_o3_val, n_valid], np.nan)

		obs_pm25_con_val_ts[:,vv] = obs_ds_val.pm25_con_val
		obs_o3_con_val_ts[:,vv] = obs_ds_val.o3_con_val

	## Create 2D arrays to hold all forecast-observation pairs
	obs_pm25_con_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	obs_pm25_con_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	obs_o3_con_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	obs_o3_con_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_pm25_con_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_pm25_con_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_o3_con_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_o3_con_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	## Put the time series obs into cycle-lead time arrays
	for cc in range(n_cycle):
		for ll in range(lead_hrs+1):
			ind_valid = np.where(dt_arr_valid == dt_arr_pairs[cc,ll])[0][0]
			obs_pm25_con_use_stn[:,cc,ll] = obs_pm25_con_use_ts[:,ind_valid]
			obs_pm25_con_val_stn[:,cc,ll] = obs_pm25_con_val_ts[:,ind_valid]
			obs_o3_con_use_stn[:,cc,ll] = obs_o3_con_use_ts[:,ind_valid]
			obs_o3_con_val_stn[:,cc,ll] = obs_o3_con_val_ts[:,ind_valid]

	## Eventually write output files over a subdomain in addition to CONUS?
	if sub_sites:
		sub_sid_o3_use_file = sub_sid_dir.joinpath('id_calif_stat_train_o3.txt')
		sub_sid_o3_val_file = sub_sid_dir.joinpath('id_calif_stat_verif_o3.txt')
		sub_sid_pm25_use_file = sub_sid_dir.joinpath('id_calif_stat_train_pm25.txt')
		sub_sid_pm25_val_file = sub_sid_dir.joinpath('id_calif_stat_verif_pm25.txt')

		sub_sid_o3_use = pd.read_csv(sub_sid_o3_use_file, sep=' ', header=None).T[0].to_numpy()
		sub_sid_o3_val = pd.read_csv(sub_sid_o3_val_file, sep=' ', header=None).T[0].to_numpy()
		sub_sid_pm25_use = pd.read_csv(sub_sid_pm25_use_file, sep=' ', header=None).T[0].to_numpy()
		sub_sid_pm25_val = pd.read_csv(sub_sid_pm25_val_file, sep=' ', header=None).T[0].to_numpy()

		inds_sub_o3_use = np.in1d(obs_o3_sid_use, sub_sid_o3_use).nonzero()[0]
		inds_sub_o3_val = np.in1d(obs_o3_sid_val, sub_sid_o3_val).nonzero()[0]
		inds_sub_pm25_use = np.in1d(obs_pm25_sid_use, sub_sid_pm25_use).nonzero()[0]
		inds_sub_pm25_val = np.in1d(obs_pm25_sid_val, sub_sid_pm25_val).nonzero()[0]

		n_obs_o3_use_sub = len(inds_sub_o3_use)
		n_obs_o3_val_sub = len(inds_sub_o3_val)
		n_obs_pm25_use_sub = len(inds_sub_pm25_use)
		n_obs_pm25_val_sub = len(inds_sub_pm25_val)

		obs_o3_sid_use_sub = obs_o3_sid_use[inds_sub_o3_use]
		obs_o3_sid_val_sub = obs_o3_sid_val[inds_sub_o3_val]
		obs_o3_lat_use_sub = obs_o3_lat_use[inds_sub_o3_use]
		obs_o3_lat_val_sub = obs_o3_lat_val[inds_sub_o3_val]
		obs_o3_lon_use_sub = obs_o3_lon_use[inds_sub_o3_use]
		obs_o3_lon_val_sub = obs_o3_lon_val[inds_sub_o3_val]
		obs_o3_cmaq_i_use_sub = obs_o3_cmaq_i_use[inds_sub_o3_use]
		obs_o3_cmaq_i_val_sub = obs_o3_cmaq_i_val[inds_sub_o3_val]
		obs_o3_cmaq_j_use_sub = obs_o3_cmaq_j_use[inds_sub_o3_use]
		obs_o3_cmaq_j_val_sub = obs_o3_cmaq_j_val[inds_sub_o3_val]

		obs_pm25_sid_use_sub = obs_pm25_sid_use[inds_sub_pm25_use]
		obs_pm25_sid_val_sub = obs_pm25_sid_val[inds_sub_pm25_val]
		obs_pm25_lat_use_sub = obs_pm25_lat_use[inds_sub_pm25_use]
		obs_pm25_lat_val_sub = obs_pm25_lat_val[inds_sub_pm25_val]
		obs_pm25_lon_use_sub = obs_pm25_lon_use[inds_sub_pm25_use]
		obs_pm25_lon_val_sub = obs_pm25_lon_val[inds_sub_pm25_val]
		obs_pm25_cmaq_i_use_sub = obs_pm25_cmaq_i_use[inds_sub_pm25_use]
		obs_pm25_cmaq_i_val_sub = obs_pm25_cmaq_i_val[inds_sub_pm25_val]
		obs_pm25_cmaq_j_use_sub = obs_pm25_cmaq_j_use[inds_sub_pm25_use]
		obs_pm25_cmaq_j_val_sub = obs_pm25_cmaq_j_val[inds_sub_pm25_val]

#		print(obs_pm25_sid_use_sub[35].values)
#		print(obs_pm25_lat_use_sub[35].values)
#		print(obs_pm25_lon_use_sub[35].values)
#		sys.exit()

		obs_o3_con_use_stn_sub = obs_o3_con_use_stn[inds_sub_o3_use,:,:]
		obs_o3_con_val_stn_sub = obs_o3_con_val_stn[inds_sub_o3_val,:,:]
		obs_pm25_con_use_stn_sub = obs_pm25_con_use_stn[inds_sub_pm25_use,:,:]
		obs_pm25_con_val_stn_sub = obs_pm25_con_val_stn[inds_sub_pm25_val,:,:]

		cmaq_o3_con_use_stn_sub = np.full([n_obs_o3_use_sub, n_cycle, lead_hrs+1], np.nan)
		cmaq_o3_con_val_stn_sub = np.full([n_obs_o3_val_sub, n_cycle, lead_hrs+1], np.nan)
		cmaq_pm25_con_use_stn_sub = np.full([n_obs_pm25_use_sub, n_cycle, lead_hrs+1], np.nan)
		cmaq_pm25_con_val_stn_sub = np.full([n_obs_pm25_val_sub, n_cycle, lead_hrs+1], np.nan)

	## Loop over CMAQ forecast cycles
	for cc in range(n_cycle):
		this_cycle_dt = dt_arr_cycle[cc]
		this_cycle_infile = this_cycle_dt.strftime(fmt_yyyymmdd+'.'+fmt_hh+'z')
		this_cycle_outfile = this_cycle_dt.strftime(fmt_yyyymmdd_hhmm)
		this_cycle_date = this_cycle_dt.strftime(fmt_yyyymmdd)
		this_cycle_yr = this_cycle_dt.strftime(fmt_yyyy)
		this_cycle_mo = this_cycle_dt.strftime(fmt_mm)

		## Second, read in AnEn bias-corrected CMAQ files (Ju-Hye has all CMAQ hourly files in one directory)
		cmaq_fname_in_this = cmaq_dir.joinpath('bc_cmaq_06z_pm2.5_o3_06z_'+this_cycle_date+'.nc')

		if cmaq_fname_in_this.is_file():
			print('Reading '+str(cmaq_fname_in_this))
			ds_cmaq = xr.open_dataset(cmaq_fname_in_this)
			## Dimensions (lead_hours 73, latitude 265, longitude 442) 

			if get_cmaq_lat_lon:
				cmaq_lat = ds_cmaq.latitude
				cmaq_lon = ds_cmaq.longitude
				get_cmaq_lat_lon = False

			cmaq_o3_this     = ds_cmaq.o3_anen_cmaq_bm3
			cmaq_pm25_this   = ds_cmaq.pm25_anen_cmaq_bm1

		else:
			print('WARNING: File '+str(cmaq_fname_in_this)+' does not exist!')

		'''
		## Now bilinearly interpolate the CMAQ O3 & PM2.5 fields to the AirNow station locations
		## CMAQ is on a Lambert conformal grid, so extra steps are required before calculating the spline

		## First, read attributes from a template CMAQ file to get projection parameters
		if get_cmaq_proj:
			cmaq_coord_file = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
			print('Reading CMAQ coordinate data from '+str(cmaq_coord_file))
			cmaq_coord_ds  = xr.open_dataset(cmaq_coord_file)
			cmaq_truelat1  = cmaq_coord_ds.attrs['P_ALP']
			cmaq_truelat2  = cmaq_coord_ds.attrs['P_BET']
			cmaq_stand_lon = cmaq_coord_ds.attrs['P_GAM']
			cmaq_cen_lat   = cmaq_coord_ds.attrs['YCENT']
			cmaq_cen_lon   = cmaq_coord_ds.attrs['XCENT']
			cmaq_dx        = cmaq_coord_ds.attrs['XCELL']
			cmaq_dy        = cmaq_coord_ds.attrs['YCELL']
			cmaq_we_dim    = cmaq_coord_ds.attrs['NCOLS']
			cmaq_ns_dim    = cmaq_coord_ds.attrs['NROWS']
			CMAQ_EARTH_RAD = 6370000
			get_cmaq_proj  = False

			## Second, interpolate AirNow locations into CMAQ (x,y) space
			## Get the projection
			p = pyproj.Proj(proj='lcc', lon_0=cmaq_cen_lon, lat_0=cmaq_cen_lat,
					lat_1=cmaq_truelat1, lat_2=cmaq_truelat2,	R=CMAQ_EARTH_RAD)
			## Get the distance in meters from the projection center point
			x_from_c_o3_use, y_from_c_o3_use = p(obs_o3_lon_use, obs_o3_lat_use)
			x_from_c_o3_val, y_from_c_o3_val = p(obs_o3_lon_val, obs_o3_lat_val)
			x_from_c_pm25_use, y_from_c_pm25_use = p(obs_pm25_lon_use, obs_pm25_lat_use)
			x_from_c_pm25_val, y_from_c_pm25_val = p(obs_pm25_lon_val, obs_pm25_lat_val)
			## Convert to the number of grid points from the center point
			i_from_c_o3_use = x_from_c_o3_use / cmaq_dx
			i_from_c_o3_val = x_from_c_o3_val / cmaq_dx
			j_from_c_o3_use = y_from_c_o3_use / cmaq_dy
			j_from_c_o3_val = y_from_c_o3_val / cmaq_dy
			i_from_c_pm25_use = x_from_c_pm25_use / cmaq_dx
			i_from_c_pm25_val = x_from_c_pm25_val / cmaq_dx
			j_from_c_pm25_use = y_from_c_pm25_use / cmaq_dy
			j_from_c_pm25_val = y_from_c_pm25_val / cmaq_dy
			## Convert to decimal interpolated grid point (start at SW corner with (1,1))
			cmaq_i_o3_use = (cmaq_we_dim // 2) + i_from_c_o3_use
			cmaq_i_o3_val = (cmaq_we_dim // 2) + i_from_c_o3_val
			cmaq_j_o3_use = (cmaq_ns_dim // 2) + j_from_c_o3_use
			cmaq_j_o3_val = (cmaq_ns_dim // 2) + j_from_c_o3_val
			cmaq_i_pm25_use = (cmaq_we_dim // 2) + i_from_c_pm25_use
			cmaq_i_pm25_val = (cmaq_we_dim // 2) + i_from_c_pm25_val
			cmaq_j_pm25_use = (cmaq_ns_dim // 2) + j_from_c_pm25_use
			cmaq_j_pm25_val = (cmaq_ns_dim // 2) + j_from_c_pm25_val
			## Assemble 1-D arrays of CMAQ grid point longitude and latitude indexes (start CMAQ indexing at (1,1))
			cmaq_i_arr = np.asarray(range(1, cmaq_we_dim+1))
			cmaq_j_arr = np.asarray(range(1, cmaq_ns_dim+1))
		'''

#		'''
		## Some useful code stubs if we decide to use nearest-neighbor "interpolation" instead of bilinear interp
		## Find the four closest grid points (i,j) and values to a given lat/lon
		print('Finding CMAQ nearest neighbor values at AirNow station locations')
		for nn in range(n_obs_o3_use):
			## Find closest j,i index
			dist = abs(cmaq_lat-obs_o3_lat_use[nn])+abs(cmaq_lon-obs_o3_lon_use[nn])
			j,i = np.unravel_index(dist.argmin(), dist.shape)
#			print('j = '+str(j)+', i = '+str(i))
#			print(cmaq_o3_this.shape)
#			print(cmaq_o3_con_use_stn.shape)
#			print(cmaq_lat.shape)
#			print('obs lat = '+str(obs_o3_lat_use.values[nn])+', obs lon = '+str(obs_o3_lon_use.values[nn]))
#			print('cmaq nn lat = '+str(cmaq_lat.values[j,i])+', cmaq nn lon = '+str(cmaq_lon.values[j,i]))
#			print(cmaq_o3_this[:,j,i].values)
			cmaq_o3_con_use_stn[nn,cc,:] = cmaq_o3_this[:,j,i]

			'''
			## Find 4 closest lat,lon points (not necessarily the four surrounding points!)
			cmaq_points = np.column_stack([cmaq_lon.values.ravel(), cmaq_lat.values.ravel()])
			tree = sp.spatial.KDTree(cmaq_points)
			print(tree)
			print(tree.data)
			stn_pt = np.column_stack([obs_o3_lon_use[nn], obs_o3_lat_use[nn]])
			dist, idx = tree.query(stn_pt, k=4)
			print(tree.data[idx])			# tuple of 4 closest lat/lon grid points
			print(tree.data[idx][0])		# 2-D array of 4 closest lat/lon grid points
			print(tree.data[idx][0,0])		# array of 1st closest lat/lon grid point
			print(tree.data[idx][0,1])		# array of 2nd closest lat/lon grid point
			print(tree.data[idx][0,0,0])	# longitude of 1st closest grid point
			ind_lon = np.where(cmaq_lon == tree.data[idx][0,0,0])
			ind_lat = np.where(cmaq_lat == tree.data[idx][0,0,1])
			print(ind_lon)
			print(ind_lat)
			print(tree.data[idx][0,:,1])
			for n in range(4):
				a = abs(cmaq_lat-tree.data[idx][0,n,1]) + abs(cmaq_lon-tree.data[idx][0,n,0])
				i,j = np.unravel_index(a.argmin(), a.shape)
				print('i = '+str(i)+', j = '+str(j))
			'''
#			sys.exit()

		for nn in range(n_obs_o3_val):
			## Find closest j,i index
			dist = abs(cmaq_lat-obs_o3_lat_val[nn])+abs(cmaq_lon-obs_o3_lon_val[nn])
			j,i = np.unravel_index(dist.argmin(), dist.shape)
			cmaq_o3_con_val_stn[nn,cc,:] = cmaq_o3_this[:,j,i]

		for nn in range(n_obs_pm25_use):
			dist = abs(cmaq_lat-obs_pm25_lat_use[nn])+abs(cmaq_lon-obs_pm25_lon_use[nn])
			j,i = np.unravel_index(dist.argmin(), dist.shape)
			cmaq_pm25_con_use_stn[nn,cc,:] = cmaq_pm25_this[:,j,i]

		for nn in range(n_obs_pm25_val):
			dist = abs(cmaq_lat-obs_pm25_lat_val[nn])+abs(cmaq_lon-obs_pm25_lon_val[nn])
			j,i = np.unravel_index(dist.argmin(), dist.shape)
			cmaq_pm25_con_val_stn[nn,cc,:] = cmaq_pm25_this[:,j,i]
#		sys.exit()

		'''
		## Bilinearly interpolate the CMAQ fields to AirNow station locations
		print('Interpolating CMAQ O3 data to AirNow training stations')
		cmaq_o3_con_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_o3_this,
			obs_o3_lat_use, obs_o3_lon_use, opt=2)
		print('Interpolating CMAQ O3 data to AirNow validation stations')
		cmaq_o3_con_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_o3_this,
			obs_o3_lat_val, obs_o3_lon_val, opt=2)
		print('Interpolating CMAQ PM2.5 data to AirNow training stations')
		cmaq_pm25_con_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pm25_this,
			obs_pm25_lat_use, obs_pm25_lon_use, opt=2)
		print('Interpolating CMAQ PM2.5 data to AirNow validation stations')
		cmaq_pm25_con_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pm25_this,
			obs_pm25_lat_val, obs_pm25_lon_val, opt=2)
		'''
		'''
		print('Interpolating CMAQ O3 data to AirNow stations')
		obs_o3_lat_combo = np.append(obs_o3_lat_use, obs_o3_lat_val)
		obs_o3_lon_combo = np.append(obs_o3_lon_use, obs_o3_lon_val)
		cmaq_o3_con_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_o3_this,
			obs_o3_lat_combo, obs_o3_lon_combo, opt=2)
		cmaq_o3_con_use_dummy = cmaq_o3_con_dummy[:,0:n_obs_o3_use]
		cmaq_o3_con_val_dummy = cmaq_o3_con_dummy[:,n_obs_o3_use:]

		print('Interpolating CMAQ PM2.5 data to AirNow stations')
		obs_pm25_lat_combo = np.append(obs_pm25_lat_use, obs_pm25_lat_val)
		obs_pm25_lon_combo = np.append(obs_pm25_lon_use, obs_pm25_lon_val)
		cmaq_pm25_con_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pm25_this,
			obs_pm25_lat_combo, obs_pm25_lon_combo, opt=2)
		cmaq_pm25_con_use_dummy = cmaq_pm25_con_dummy[:,0:n_obs_pm25_use]
		cmaq_pm25_con_val_dummy = cmaq_pm25_con_dummy[:,n_obs_pm25_use:]
		'''
			
		'''
		for ll in range(lead_hrs+1):
			cmaq_o3_con_use_stn[:,cc,ll] = cmaq_o3_con_use_dummy[ll,:]
			cmaq_o3_con_val_stn[:,cc,ll] = cmaq_o3_con_val_dummy[ll,:]
			cmaq_pm25_con_use_stn[:,cc,ll] = cmaq_pm25_con_use_dummy[ll,:]
			cmaq_pm25_con_val_stn[:,cc,ll] = cmaq_pm25_con_val_dummy[ll,:]
		'''

		'''
		## This will necessarily impact values near the boundary of real data, so be careful at the edges.
		## To keep extracted values of the spline positive, set the degree of the spline to 1 (true bilinear interp).
		## If the AirNow lat/lon arrays only have a nan (i.e., no valid data), the output is a nan array as well.
		## Loop over lead times, as RectBivariateSpline can only handle a 2D array
#		for ll in range(lead_hrs+1):
#			cmaq_o3_spline = sp.interpolate.RectBivariateSpline(
#											cmaq_j_arr, cmaq_i_arr, cmaq_o3_this[ll,:,:], kx=1, ky=1)
#			cmaq_o3_con_use_stn[:,cc,ll] = cmaq_o3_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
#			cmaq_o3_con_val_stn[:,cc,ll] = cmaq_o3_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)

#			cmaq_pm25_spline = sp.interpolate.RectBivariateSpline(
#										cmaq_j_arr, cmaq_i_arr, cmaq_pm25_this[ll,:,:], kx=1, ky=1)
#			cmaq_pm25_con_use_stn[:,cc,ll] = cmaq_pm25_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
#			cmaq_pm25_con_val_stn[:,cc,ll] = cmaq_pm25_spline.ev(cmaq_j_pm25_val, cmaq_j_pm25_val)


#			print(cmaq_pm25_this[ll,:,:].shape)
#			for nn in range(n_obs_pm25_val):
#				print('j = '+str(cmaq_j_pm25_val[nn])+', i = '+str(cmaq_i_pm25_val[nn]))
#				print('lon = '+str(obs_pm25_lon_val[nn].values)+', lat = '+str(obs_pm25_lat_val[nn].values)+', id = '+str(obs_pm25_sid_val[nn].values))
#			print(cmaq_pm25_spline)
#			sys.exit()

#			print(cmaq_pm25_this[ll,:,:].shape)
#			ig, jg = np.meshgrid(cmaq_i_arr, cmaq_j_arr, indexing='ij')
#			print(ig.shape)
#			print(jg.shape)
#			ig = np.transpose(ig)
#			jg = np.transpose(jg)
#			print(ig.shape)
#			print(jg.shape)
#			sys.exit()
#			cmaq_pm25_interp = sp.interpolate.RegularGridInterpolator(
#										(cmaq_j_arr, cmaq_i_arr), cmaq_pm25_this[ll,:,:], method='linear', bounds_error=False)
#			print('cmaq_j_pm25_val:')
#			print(cmaq_j_pm25_val)
#			print('cmaq_i_pm25_val:')
#			print(cmaq_i_pm25_val)
#			print('combined array:')
#			test_points = np.array([cmaq_j_pm25_val, cmaq_i_pm25_val]).T
#			print(test_points.shape)
#			print(test_points)
#			print(cmaq_pm25_this[ll,:,:].shape)
#			print('min, max of cmaq_j_pm25_val')
#			print(min(cmaq_i_pm25_val))
#			print(max(cmaq_j_pm25_val))
#			print('min, max of cmaq_i_pm25_val')
#			print(min(cmaq_i_pm25_val))
#			print(max(cmaq_i_pm25_val))
#			print(obs_pm25_lat_val[24].values)
#			print(obs_pm25_lon_val[24].values)
			
#			print(cmaq_pm25_interp)
#			print(cmaq_pm25_interp(test_points))
#			cmaq_pm25_interpn = sp.interpolate.interpn(
#										(cmaq_j_arr, cmaq_i_arr), cmaq_pm25_this[ll,:,:], test_points, bounds_error=False)
#			print(cmaq_pm25_interpn)

#			print(cmaq_lat)
#			print(cmaq_lon)
#			print(cmaq_pm25_this[ll,:,:])
#			print(obs_pm25_lat_val)
#			print(obs_pm25_lon_val)
#			cmaq_pm25_interp_val = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pm25_this,
#				obs_pm25_lat_val, obs_pm25_lon_val, opt=2)
#			print(cmaq_pm25_interp_val[0,:].values)
#			print(cmaq_pm25_interp_val.shape)
		'''

		## Create xarray data arrays using the CMAQ and obs data
		n_use_o3 = len(obs_o3_lat_use)
		n_val_o3 = len(obs_o3_lat_val)
		n_use_pm25 = len(obs_pm25_lat_use)
		n_val_pm25 = len(obs_pm25_lat_val)
		n_leads = lead_hrs+1

		if sub_sites:
			n_use_o3_sub = len(obs_o3_lat_use_sub)
			n_val_o3_sub = len(obs_o3_lat_val_sub)
			n_use_pm25_sub = len(obs_pm25_lat_use_sub)
			n_val_pm25_sub = len(obs_pm25_lat_val_sub)

			cmaq_o3_con_use_stn_sub = cmaq_o3_con_use_stn[inds_sub_o3_use,:,:]
			cmaq_o3_con_val_stn_sub = cmaq_o3_con_val_stn[inds_sub_o3_val,:,:]
			cmaq_pm25_con_use_stn_sub = cmaq_pm25_con_use_stn[inds_sub_pm25_use,:,:]
			cmaq_pm25_con_val_stn_sub = cmaq_pm25_con_val_stn[inds_sub_pm25_val,:,:]

#			print(obs_pm25_sid_val_sub[9].values)	# 60670010
#			print(obs_pm25_lat_val_sub[9].values)	# 38.57
#			print(obs_pm25_lon_val_sub[9].values)	# -121.49
#			print('j=142.650, i=46.676:')
#			print(cmaq_pm25_con_val_stn_sub[9,:,:])	# 20210201_06f00: 13.57

#			print('j=142, i=46:')
#			print(cmaq_pm25_this[:,142,46].values)		# 20210201_06f00: 8.52
#			print('j=142, i=47:')
#			print(cmaq_pm25_this[:,142,47].values)		# 20210201_06f00: 9.15
#			print('j=143, i=47:')
#			print(cmaq_pm25_this[:,143,47].values)		# 20210201_06f00: 9.04
#			print('j=143, i=46:')
#			print(cmaq_pm25_this[:,143,46].values)		# 20210201_06f00: 8.57
#			print('')
#			print('point j=152, i=34:')
#			print('lat = '+str(cmaq_lat[152,34].values)+', lon = '+str(cmaq_lon[152,34].values))
#			print(cmaq_pm25_this[:,152,34].values)
#			print('point j=152, i=35:')
#			print('lat = '+str(cmaq_lat[152,35].values)+', lon = '+str(cmaq_lon[152,35].values))
#			print(cmaq_pm25_this[:,152,35].values)
#			print('point j=153, i=35:')
#			print('lat = '+str(cmaq_lat[153,35].values)+', lon = '+str(cmaq_lon[153,35].values))
#			print(cmaq_pm25_this[:,153,35].values)
#			print('point j=153, i=34:')
#			print('lat = '+str(cmaq_lat[153,34].values)+', lon = '+str(cmaq_lon[153,34].values))
#			print(cmaq_pm25_this[:,153,34].values)

			## Find the closest grid point (i,j) to a given lat/lon
#			a = abs(cmaq_lat-obs_pm25_lat_val_sub[9])+abs(cmaq_lon-obs_pm25_lon_val_sub[9])
#			i,j = np.unravel_index(a.argmin(),a.shape)
#			print(i)
#			print(j)

		cmaq_str = 'CMAQ gridded forecast bias-corrected with BM models'

		obs_o3_con_use = xr.DataArray(obs_o3_con_use_stn[:,cc,:],
								coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
								attrs={'description':'AirNow O3 concentration observations at stations used for merging and training', 'units':'ppbv'})
		obs_o3_con_val = xr.DataArray(obs_o3_con_val_stn[:,cc,:],
								coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
								attrs={'description':'AirNow O3 concentration observations at stations withheld for validation', 'units':'ppbv'})
		obs_pm25_con_use = xr.DataArray(obs_pm25_con_use_stn[:,cc,:],
								coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
								attrs={'description':'AirNow PM2.5 concentration observations at stations used for merging and training', 'units':'ug m-3'})
		obs_pm25_con_val = xr.DataArray(obs_pm25_con_val_stn[:,cc,:],
								coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
								attrs={'description':'AirNow PM2.5 concentration observations at stations withheld for validation', 'units':'ug m-3'})

		cmaq_o3_con_use = xr.DataArray(cmaq_o3_con_use_stn[:,cc,:],
									coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
									attrs={'description':cmaq_str+' O3 concentration interpolated to AirNow O3 stations used for merging and training', 'units':'ppbv'})
		cmaq_o3_con_val = xr.DataArray(cmaq_o3_con_val_stn[:,cc,:],
									coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
									attrs={'description':cmaq_str+' O3 concentration interpolated to AirNow O3 stations withheld for validation', 'units':'ppbv'})

		cmaq_pm25_con_use = xr.DataArray(cmaq_pm25_con_use_stn[:,cc,:],
									coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
									attrs={'description':cmaq_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations used for merging and training', 'units':'ug m-3'})
		cmaq_pm25_con_val = xr.DataArray(cmaq_pm25_con_val_stn[:,cc,:],
									coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
									attrs={'description':cmaq_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations withheld for validation', 'units':'ug m-3'})


		if sub_sites:
			n_use_o3_sub = len(obs_o3_lat_use_sub)
			n_val_o3_sub = len(obs_o3_lat_val_sub)
			n_use_pm25_sub = len(obs_pm25_lat_use_sub)
			n_val_pm25_sub = len(obs_pm25_lat_val_sub)

			cmaq_str = 'AnEn bias-corrected CMAQ gridded forecast over a CA/NV subdomain'

			obs_o3_con_use_sub = xr.DataArray(obs_o3_con_use_stn_sub[:,cc,:],
					coords={'n_obs_use_o3_sub':n_use_o3_sub,'lead_hours':n_leads}, dims=['n_use_o3_sub','n_leads'],
					attrs={'description':'AirNow O3 concentration observations at stations used for merging and training', 'units':'ppbv'})
			obs_o3_con_val_sub = xr.DataArray(obs_o3_con_val_stn_sub[:,cc,:],
					coords={'n_obs_val_o3_sub':n_val_o3_sub,'lead_hours':n_leads}, dims=['n_val_o3_sub','n_leads'],
					attrs={'description':'AirNow O3 concentration observations at stations withheld for validation', 'units':'ppbv'})
			obs_pm25_con_use_sub = xr.DataArray(obs_pm25_con_use_stn_sub[:,cc,:],
					coords={'n_obs_use_pm25_sub':n_use_pm25_sub,'lead_hours':n_leads}, dims=['n_use_pm25_sub','n_leads'],
					attrs={'description':'AirNow PM2.5 concentration observations at stations used for merging and training', 'units':'ug m-3'})
			obs_pm25_con_val_sub = xr.DataArray(obs_pm25_con_val_stn_sub[:,cc,:],
					coords={'n_obs_val_pm25_sub':n_val_pm25_sub,'lead_hours':n_leads}, dims=['n_val_pm25_sub','n_leads'],
					attrs={'description':'AirNow PM2.5 concentration observations at stations withheld for validation', 'units':'ug m-3'})

			cmaq_o3_con_use_sub = xr.DataArray(cmaq_o3_con_use_stn_sub[:,cc,:],
					coords={'n_obs_use_o3_sub':n_use_o3_sub,'lead_hours':n_leads}, dims=['n_use_o3_sub','n_leads'],
					attrs={'description':cmaq_str+' O3 concentration interpolated to AirNow O3 stations used for merging and training', 'units':'ppbv'})
			cmaq_o3_con_val_sub = xr.DataArray(cmaq_o3_con_val_stn_sub[:,cc,:],
					coords={'n_obs_val_o3_sub':n_val_o3_sub,'lead_hours':n_leads}, dims=['n_val_o3_sub','n_leads'],
					attrs={'description':cmaq_str+' O3 concentration interpolated to AirNow O3 stations withheld for validation', 'units':'ppbv'})

			cmaq_pm25_con_use_sub = xr.DataArray(cmaq_pm25_con_use_stn_sub[:,cc,:],
					coords={'n_obs_use_pm25_sub':n_use_pm25_sub,'lead_hours':n_leads}, dims=['n_use_pm25_sub','n_leads'],
					attrs={'description':cmaq_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations used for merging and training', 'units':'ug m-3'})
			cmaq_pm25_con_val_sub = xr.DataArray(cmaq_pm25_con_val_stn_sub[:,cc,:],
					coords={'n_obs_val_pm25_sub':n_val_pm25_sub,'lead_hours':n_leads}, dims=['n_val_pm25_sub','n_leads'],
					attrs={'description':cmaq_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations withheld for validation', 'units':'ug m-3'})

		## Create new xarray datasets
		ds_use = xr.Dataset(
						data_vars={ 'cmaq_pm25_con_use':cmaq_pm25_con_use, 'cmaq_o3_con_use':cmaq_o3_con_use,
										'obs_pm25_con_use':obs_pm25_con_use, 'obs_o3_con_use':obs_o3_con_use,
										'obs_pm25_sid_use':obs_pm25_sid_use, 'obs_o3_sid_use':obs_o3_sid_use,
										'obs_pm25_lon_use':obs_pm25_lon_use, 'obs_o3_lon_use':obs_o3_lon_use,
										'obs_pm25_lat_use':obs_pm25_lat_use, 'obs_o3_lat_use':obs_o3_lat_use,
										'obs_pm25_cmaq_i_use':obs_pm25_cmaq_i_use, 'obs_o3_cmaq_i_use':obs_o3_cmaq_i_use,
										'obs_pm25_cmaq_j_use':obs_pm25_cmaq_j_use, 'obs_o3_cmaq_j_use':obs_o3_cmaq_j_use,
							},
						coords={'n_obs_use_pm25':n_use_pm25, 'n_obs_use_o3':n_use_o3, 'lead_hours':n_leads},
						attrs={'description':'Observations and interpolated CMAQ gridded forecast variables at AirNow PM2.5 & O3 stations retained for use in merging and training'},
						)

		ds_val = xr.Dataset(
						data_vars={ 'cmaq_pm25_con_val':cmaq_pm25_con_val, 'cmaq_o3_con_val':cmaq_o3_con_val,
										'obs_pm25_con_val':obs_pm25_con_val, 'obs_o3_con_val':obs_o3_con_val,
										'obs_pm25_sid_val':obs_pm25_sid_val, 'obs_o3_sid_val':obs_o3_sid_val,
										'obs_pm25_lon_val':obs_pm25_lon_val, 'obs_o3_lon_val':obs_o3_lon_val,
										'obs_pm25_lat_val':obs_pm25_lat_val, 'obs_o3_lat_val':obs_o3_lat_val,
										'obs_pm25_cmaq_i_val':obs_pm25_cmaq_i_val, 'obs_o3_cmaq_i_val':obs_o3_cmaq_i_val,
										'obs_pm25_cmaq_j_val':obs_pm25_cmaq_j_val, 'obs_o3_cmaq_j_val':obs_o3_cmaq_j_val,
							},
						coords={'n_obs_val_pm25':n_val_pm25, 'n_obs_val_o3':n_val_o3, 'lead_hours':n_leads},
						attrs={'description':'Observations and interpolated CMAQ gridded forecast variables at AirNow PM2.5 & O3 stations withheld for validation'},
						)

		if sub_sites:
			ds_use_sub = xr.Dataset(
					data_vars={ 'cmaq_pm25_con_use':cmaq_pm25_con_use_sub, 'cmaq_o3_con_use':cmaq_o3_con_use_sub,
									'obs_pm25_con_use':obs_pm25_con_use_sub, 'obs_o3_con_use':obs_o3_con_use_sub,
									'obs_pm25_sid_use':obs_pm25_sid_use_sub, 'obs_o3_sid_use':obs_o3_sid_use_sub,
									'obs_pm25_lon_use':obs_pm25_lon_use_sub, 'obs_o3_lon_use':obs_o3_lon_use_sub,
									'obs_pm25_lat_use':obs_pm25_lat_use_sub, 'obs_o3_lat_use':obs_o3_lat_use_sub,
									'obs_pm25_cmaq_i_use':obs_pm25_cmaq_i_use_sub, 'obs_o3_cmaq_i_use':obs_o3_cmaq_i_use_sub,
									'obs_pm25_cmaq_j_use':obs_pm25_cmaq_j_use_sub, 'obs_o3_cmaq_j_use':obs_o3_cmaq_j_use_sub,
						},
					coords={'n_obs_use_pm25_sub':n_use_pm25_sub, 'n_obs_use_o3_sub':n_use_o3_sub, 'lead_hours':n_leads},
					attrs={'description':'Observations and interpolated CMAQ gridded forecast variables at AirNow PM2.5 & O3 stations retained for use in merging and training (California subdomain)'},
					)

			ds_val_sub = xr.Dataset(
					data_vars={ 'cmaq_pm25_con_val':cmaq_pm25_con_val_sub, 'cmaq_o3_con_val':cmaq_o3_con_val_sub,
									'obs_pm25_con_val':obs_pm25_con_val_sub, 'obs_o3_con_val':obs_o3_con_val_sub,
									'obs_pm25_sid_val':obs_pm25_sid_val_sub, 'obs_o3_sid_val':obs_o3_sid_val_sub,
									'obs_pm25_lon_val':obs_pm25_lon_val_sub, 'obs_o3_lon_val':obs_o3_lon_val_sub,
									'obs_pm25_lat_val':obs_pm25_lat_val_sub, 'obs_o3_lat_val':obs_o3_lat_val_sub,
									'obs_pm25_cmaq_i_val':obs_pm25_cmaq_i_val_sub, 'obs_o3_cmaq_i_val':obs_o3_cmaq_i_val_sub,
									'obs_pm25_cmaq_j_val':obs_pm25_cmaq_j_val_sub, 'obs_o3_cmaq_j_val':obs_o3_cmaq_j_val_sub,
						},
					coords={'n_obs_val_pm25_sub':n_val_pm25_sub, 'n_obs_val_o3_sub':n_val_o3_sub, 'lead_hours':n_leads},
					attrs={'description':'Observations and interpolated CMAQ gridded forecast variables at AirNow PM2.5 & O3 stations withheld for validation (California subdomain)'},
					)

		## Set the output paths & filenames
		out_dir = out_dir_parent.joinpath(this_cycle_yr,this_cycle_mo)
		out_dir.mkdir(parents=True, exist_ok=True)
		out_dir.chmod(0o755)
		fname_use = out_dir.joinpath('conus_anen_bc_cmaq_airnow_pm2.5_o3_'+this_cycle_outfile+'_use.nc')
		fname_val = out_dir.joinpath('conus_anen_bc_cmaq_airnow_pm2.5_o3_'+this_cycle_outfile+'_val.nc')
		fname_use_sub = out_dir.joinpath('calif_anen_bc_cmaq_airnow_pm2.5_o3_'+this_cycle_outfile+'_use.nc')
		fname_val_sub = out_dir.joinpath('calif_anen_bc_cmaq_airnow_pm2.5_o3_'+this_cycle_outfile+'_val.nc')

		## Write the datasets to NetCDF
		print('Writing '+str(fname_use))
		ds_use.to_netcdf(fname_use)
		fname_use.chmod(0o644)
		print('Writing '+str(fname_val))
		ds_val.to_netcdf(fname_val)
		fname_val.chmod(0o644)

		if sub_sites:
			print('Writing '+str(fname_use_sub))
			ds_use_sub.to_netcdf(fname_use_sub)
			fname_use_sub.chmod(0o644)
			print('Writing '+str(fname_val_sub))
			ds_val_sub.to_netcdf(fname_val_sub)
			fname_val_sub.chmod(0o644)


if __name__ == '__main__':
	now_time_beg = dt.datetime.utcnow()
	date_beg, date_end, lead_hrs = parse_args()
	main(date_beg, date_end, lead_hrs)
	now_time_end = dt.datetime.utcnow()
	run_time_tot = now_time_end - now_time_beg
	now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
	now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
	print('\nScript completed successfully.')
	print('   Beg time: '+now_time_beg_str)
	print('   End time: '+now_time_end_str)
	print('   Run time: '+str(run_time_tot)+'\n')
