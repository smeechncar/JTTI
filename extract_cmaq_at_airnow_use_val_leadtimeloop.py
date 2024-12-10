'''
extract_cmaq_at_airnow_use_val_leadtimeloop.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 8 Sep 2022

This script first reads in AirNow obs NetCDF files that have already been separated into use & val
sets, for use in merging & validation, respectively. Second, it reads in a CMAQ forecast file to
get the grid coordinates of those stations and spatially interpolate the gridded values to the
AirNow sites via bilinear interpolation. CMAQ PM2.5 & O3 files are produced hourly. The output
files with CMAQ-interpolated values at AirNow sites are in NetCDF.

23 Dec 2022:
-- Changed the interpolation from scipy.interpolate.RectBivariateSpline to geocat.comp.rcm2points.
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

	## Set directories
	if sites_vary:
		obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_vary')
		out_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cmaq','raw','airnow','sites_vary')
	if sites_static:
		obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static','20200801_00-20211231_23')
		out_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cmaq','raw','airnow','sites_static')
	cmaq_dir = pathlib.Path('/','glade','scratch','jkim','CMAQ','cmaq_06z_prep_anen')

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
#	n_times = len(dt_arr_valid)

	## Create 2D array that holds all valid times for each model cycle
	dt_arr_pairs = np.full([n_cycle, lead_hrs+1], None)
	for cc in range(n_cycle):
		this_cycle = dt_arr_cycle[cc]
		for vv in range(lead_hrs+1):
			this_valid = this_cycle + dt.timedelta(hours=vv)
			dt_arr_pairs[cc,vv] = this_valid
#	print(dt_arr_pairs)
#	sys.exit()

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

	cmaq_t2m_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_t2m_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_t2m_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_t2m_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_ws10m_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_ws10m_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_ws10m_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_ws10m_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_wd10m_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_wd10m_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_wd10m_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_wd10m_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_swdown_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_swdown_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_swdown_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_swdown_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_pbl_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_pbl_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_pbl_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_pbl_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_nox_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_nox_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_nox_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_nox_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_noy_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_noy_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_noy_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_noy_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	## Put the time series obs into cycle-lead time arrays
	for cc in range(n_cycle):
		for ll in range(lead_hrs+1):
			ind_valid = np.where(dt_arr_valid == dt_arr_pairs[cc,ll])[0][0]
			obs_pm25_con_use_stn[:,cc,ll] = obs_pm25_con_use_ts[:,ind_valid]
			obs_pm25_con_val_stn[:,cc,ll] = obs_pm25_con_val_ts[:,ind_valid]
			obs_o3_con_use_stn[:,cc,ll] = obs_o3_con_use_ts[:,ind_valid]
			obs_o3_con_val_stn[:,cc,ll] = obs_o3_con_val_ts[:,ind_valid]

	## Loop over CMAQ forecast cycles
	for cc in range(n_cycle):
		this_cycle_dt = dt_arr_cycle[cc]
		this_cycle_infile = this_cycle_dt.strftime(fmt_yyyymmdd+'.'+fmt_hh+'z')
		this_cycle_outfile = this_cycle_dt.strftime(fmt_yyyymmdd_hhmm)
		this_cycle_yr = this_cycle_dt.strftime(fmt_yyyy)
		this_cycle_mo = this_cycle_dt.strftime(fmt_mm)

		## Second, read in CMAQ files (Ju-Hye has all CMAQ hourly files in one directory)
		cmaq_fname_in_this = cmaq_dir.joinpath('cmaq_9variables_'+this_cycle_infile+'.nc')

		if cmaq_fname_in_this.is_file():
			print('Reading '+str(cmaq_fname_in_this))
			ds_cmaq = xr.open_dataset(cmaq_fname_in_this)
			## Dimensions (lead_hours 73, latitude 265, longitude 442) 

			if get_cmaq_lat_lon:
				cmaq_lat = ds_cmaq.latitude
				cmaq_lon = ds_cmaq.longitude
				get_cmaq_lat_lon = False

			cmaq_o3_this     = ds_cmaq.O3
			cmaq_pm25_this   = ds_cmaq.PM25
			cmaq_t2m_this    = ds_cmaq.TEMP2
			cmaq_ws10m_this  = ds_cmaq.WSPD10
			cmaq_wd10m_this  = ds_cmaq.WDIR10
			cmaq_swdown_this = ds_cmaq.RGRND
			cmaq_pbl_this    = ds_cmaq.PBL
			cmaq_nox_this    = ds_cmaq.NOX
			cmaq_noy_this    = ds_cmaq.NOY

			## Need to calculate U and V components for later interpolation of wind direction
			cmaq_md10m_this = 270 - cmaq_wd10m_this	# get math wind direction from weather wind direction
			cmaq_u10m_this = cmaq_ws10m_this * np.cos(np.deg2rad(cmaq_md10m_this))
			cmaq_v10m_this = cmaq_ws10m_this * np.sin(np.deg2rad(cmaq_md10m_this))

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

		## This will necessarily impact values near the boundary of real data, so be careful at the edges.
		## To keep extracted values of the spline positive, set the degree of the spline to 1 (true bilinear interp).
		## If the AirNow lat/lon arrays only have a nan (i.e., no valid data), the output is a nan array as well.
		## Loop over lead times, as RectBivariateSpline can only handle a 2D array
		for ll in range(lead_hrs+1):
			cmaq_o3_spline = sp.interpolate.RectBivariateSpline(
											cmaq_j_arr, cmaq_i_arr, cmaq_o3_this[ll,:,:], kx=1, ky=1)
			cmaq_o3_con_use_stn[:,cc,ll] = cmaq_o3_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
			cmaq_o3_con_val_stn[:,cc,ll] = cmaq_o3_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)

			cmaq_pm25_spline = sp.interpolate.RectBivariateSpline(
										cmaq_j_arr, cmaq_i_arr, cmaq_pm25_this[ll,:,:], kx=1, ky=1)
			cmaq_pm25_con_use_stn[:,cc,ll] = cmaq_pm25_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
			cmaq_pm25_con_val_stn[:,cc,ll] = cmaq_pm25_spline.ev(cmaq_j_pm25_val, cmaq_j_pm25_val)

			cmaq_t2m_spline = sp.interpolate.RectBivariateSpline(
										cmaq_j_arr, cmaq_i_arr, cmaq_t2m_this[ll,:,:], kx=1, ky=1)
			cmaq_t2m_o3_use_stn[:,cc,ll]   = cmaq_t2m_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
			cmaq_t2m_o3_val_stn[:,cc,ll]   = cmaq_t2m_spline.ev(cmaq_j_o3_val, cmaq_j_o3_val)
			cmaq_t2m_pm25_use_stn[:,cc,ll] = cmaq_t2m_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
			cmaq_t2m_pm25_val_stn[:,cc,ll] = cmaq_t2m_spline.ev(cmaq_j_pm25_val, cmaq_j_pm25_val)

			## Bilinear interpolation of a circular variable like wind direction may lead to unexpected results
			## Instead, interpolate the U & V components to the station locations, then calculate wind speed/direction
			cmaq_u10m_spline = sp.interpolate.RectBivariateSpline(
										cmaq_j_arr, cmaq_i_arr, cmaq_u10m_this[ll,:,:], kx=1, ky=1)
			cmaq_v10m_spline = sp.interpolate.RectBivariateSpline(
										cmaq_j_arr, cmaq_i_arr, cmaq_v10m_this[ll,:,:], kx=1, ky=1)
			cmaq_u10m_o3_use = cmaq_u10m_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
			cmaq_u10m_o3_val = cmaq_u10m_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)
			cmaq_v10m_o3_use = cmaq_v10m_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
			cmaq_v10m_o3_val = cmaq_v10m_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)
			cmaq_u10m_pm25_use = cmaq_u10m_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
			cmaq_u10m_pm25_val = cmaq_u10m_spline.ev(cmaq_j_pm25_val, cmaq_i_pm25_val)
			cmaq_v10m_pm25_use = cmaq_v10m_spline.ev(cmaq_i_pm25_use, cmaq_i_pm25_use)
			cmaq_v10m_pm25_val = cmaq_v10m_spline.ev(cmaq_i_pm25_val, cmaq_i_pm25_val)

#			cmaq_ws10m_spline = sp.interpolate.RectBivariateSpline(
#										cmaq_j_arr, cmaq_i_arr, cmaq_ws10m_this[ll,:,:], kx=1, ky=1)
#			cmaq_ws10m_o3_use[:,cc,ll]   = cmaq_ws10m_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
#			cmaq_ws10m_o3_val[:,cc,ll]   = cmaq_ws10m_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)
#			cmaq_ws10m_pm25_use[:,cc,ll] = cmaq_ws10m_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
#			cmaq_ws10m_pm25_val[:,cc,ll] = cmaq_ws10m_spline.ev(cmaq_j_pm25_val, cmaq_j_pm25_use)

			cmaq_ws10m_o3_use_stn[:,cc,ll]   = np.sqrt(np.square(cmaq_u10m_o3_use) + np.square(cmaq_v10m_o3_use))
			cmaq_ws10m_o3_val_stn[:,cc,ll]   = np.sqrt(np.square(cmaq_u10m_o3_val) + np.square(cmaq_v10m_o3_val))
			cmaq_ws10m_pm25_use_stn[:,cc,ll] = np.sqrt(np.square(cmaq_u10m_pm25_use) + np.square(cmaq_v10m_pm25_use))
			cmaq_ws10m_pm25_val_stn[:,cc,ll] = np.sqrt(np.square(cmaq_u10m_pm25_val) + np.square(cmaq_v10m_pm25_val))

			cmaq_wd10m_o3_use_stn[:,cc,ll]   = 270.0 - np.rad2deg(np.arctan2(cmaq_v10m_o3_use, cmaq_u10m_o3_use))
			cmaq_wd10m_o3_val_stn[:,cc,ll]   = 270.0 - np.rad2deg(np.arctan2(cmaq_v10m_o3_val, cmaq_u10m_o3_val))
			cmaq_wd10m_pm25_use_stn[:,cc,ll] = 270.0 - np.rad2deg(np.arctan2(cmaq_v10m_pm25_use, cmaq_u10m_pm25_use))
			cmaq_wd10m_pm25_val_stn[:,cc,ll] = 270.0 - np.rad2deg(np.arctan2(cmaq_v10m_pm25_val, cmaq_u10m_pm25_val))

			cmaq_swdown_spline = sp.interpolate.RectBivariateSpline(
											cmaq_j_arr, cmaq_i_arr, cmaq_swdown_this[ll,:,:], kx=1, ky=1)
			cmaq_swdown_o3_use_stn[:,cc,ll]   = cmaq_swdown_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
			cmaq_swdown_o3_val_stn[:,cc,ll]   = cmaq_swdown_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)
			cmaq_swdown_pm25_use_stn[:,cc,ll] = cmaq_swdown_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
			cmaq_swdown_pm25_val_stn[:,cc,ll] = cmaq_swdown_spline.ev(cmaq_j_pm25_val, cmaq_i_pm25_val)

			cmaq_pbl_spline = sp.interpolate.RectBivariateSpline(
										cmaq_j_arr, cmaq_i_arr, cmaq_pbl_this[ll,:,:], kx=1, ky=1)
			cmaq_pbl_o3_use_stn[:,cc,ll]   = cmaq_pbl_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
			cmaq_pbl_o3_val_stn[:,cc,ll]   = cmaq_pbl_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)
			cmaq_pbl_pm25_use_stn[:,cc,ll] = cmaq_pbl_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
			cmaq_pbl_pm25_val_stn[:,cc,ll] = cmaq_pbl_spline.ev(cmaq_j_pm25_val, cmaq_i_pm25_val)

			cmaq_nox_spline = sp.interpolate.RectBivariateSpline(
										cmaq_j_arr, cmaq_i_arr, cmaq_nox_this[ll,:,:], kx=1, ky=1)
			cmaq_nox_o3_use_stn[:,cc,ll]   = cmaq_nox_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
			cmaq_nox_o3_val_stn[:,cc,ll]   = cmaq_nox_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)
			cmaq_nox_pm25_use_stn[:,cc,ll] = cmaq_nox_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
			cmaq_nox_pm25_val_stn[:,cc,ll] = cmaq_nox_spline.ev(cmaq_j_pm25_val, cmaq_i_pm25_val)

			cmaq_noy_spline = sp.interpolate.RectBivariateSpline(
										cmaq_j_arr, cmaq_i_arr, cmaq_noy_this[ll,:,:], kx=1, ky=1)
			cmaq_noy_o3_use_stn[:,cc,ll]   = cmaq_noy_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
			cmaq_noy_o3_val_stn[:,cc,ll]   = cmaq_noy_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)
			cmaq_noy_pm25_use_stn[:,cc,ll] = cmaq_noy_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
			cmaq_noy_pm25_val_stn[:,cc,ll] = cmaq_noy_spline.ev(cmaq_j_pm25_val, cmaq_i_pm25_val)

		'''
		'''
		print('Interpolating CMAQ O3 data to AirNow stations')
		cmaq_o3_con_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_o3_this,
			obs_o3_lat_use, obs_o3_lon_use, opt=2)
		cmaq_o3_con_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_o3_this,
			obs_o3_lat_val, obs_o3_lon_val, opt=2)

		print('Interpolating CMAQ PM2.5 data to AirNow stations')
		cmaq_pm25_con_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pm25_this,
			obs_pm25_lat_use, obs_pm25_lon_use, opt=2)
		cmaq_pm25_con_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pm25_this,
			obs_pm25_lat_val, obs_pm25_lon_val, opt=2)

		print('Interpolating CMAQ T2 data to AirNow stations')
		cmaq_t2m_o3_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_t2m_this,
			obs_o3_lat_use, obs_o3_lon_use, opt=2)
		cmaq_t2m_o3_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_t2m_this,
			obs_o3_lat_val, obs_o3_lon_val, opt=2)
		cmaq_t2m_pm25_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_t2m_this,
			obs_pm25_lat_use, obs_pm25_lon_use, opt=2)
		cmaq_t2m_pm25_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_t2m_this,
			obs_pm25_lat_val, obs_pm25_lon_val, opt=2)

		print('Interpolating CMAQ U10 data to AirNow stations')
		cmaq_u10m_o3_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_u10m_this,
			obs_o3_lat_use, obs_o3_lon_use, opt=2)
		cmaq_u10m_o3_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_u10m_this,
			obs_o3_lat_val, obs_o3_lon_val, opt=2)
		cmaq_u10m_pm25_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_u10m_this,
			obs_pm25_lat_use, obs_pm25_lon_use, opt=2)
		cmaq_u10m_pm25_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_u10m_this,
			obs_pm25_lat_val, obs_pm25_lon_val, opt=2)

		print('Interpolating CMAQ V10 data to AirNow stations')
		cmaq_v10m_o3_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_u10m_this,
			obs_o3_lat_use, obs_o3_lon_use, opt=2)
		cmaq_v10m_o3_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_u10m_this,
			obs_o3_lat_val, obs_o3_lon_val, opt=2)
		cmaq_v10m_pm25_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_u10m_this,
			obs_pm25_lat_use, obs_pm25_lon_use, opt=2)
		cmaq_v10m_pm25_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_u10m_this,
			obs_pm25_lat_val, obs_pm25_lon_val, opt=2)

		cmaq_ws10m_o3_use_dummy = np.sqrt(np.square(cmaq_u10m_o3_use_dummy)+np.square(cmaq_v10m_o3_use_dummy))
		cmaq_ws10m_o3_val_dummy = np.sqrt(np.square(cmaq_u10m_o3_val_dummy)+np.square(cmaq_v10m_o3_val_dummy))
		cmaq_ws10m_pm25_use_dummy = np.sqrt(np.square(cmaq_u10m_pm25_use_dummy)+np.square(cmaq_v10m_pm25_use_dummy))
		cmaq_ws10m_pm25_val_dummy = np.sqrt(np.square(cmaq_v10m_pm25_val_dummy)+np.square(cmaq_v10m_pm25_val_dummy))

		cmaq_wd10m_o3_use_dummy = 270.0 - np.rad2deg(np.arctan2(cmaq_v10m_o3_use_dummy, cmaq_u10m_o3_use_dummy))
		cmaq_wd10m_o3_val_dummy = 270.0 - np.rad2deg(np.arctan2(cmaq_v10m_o3_val_dummy, cmaq_u10m_o3_val_dummy))
		cmaq_wd10m_pm25_use_dummy = 270.0 - np.rad2deg(np.arctan2(cmaq_v10m_pm25_use_dummy, cmaq_u10m_pm25_use_dummy))
		cmaq_wd10m_pm25_val_dummy = 270.0 - np.rad2deg(np.arctan2(cmaq_v10m_pm25_val_dummy, cmaq_u10m_pm25_val_dummy))

		print('Interpolating CMAQ SWDOWN data to AirNow stations')
		cmaq_swdown_o3_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_swdown_this,
			obs_o3_lat_use, obs_o3_lon_use, opt=2)
		cmaq_swdown_o3_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_swdown_this,
			obs_o3_lat_val, obs_o3_lon_val, opt=2)
		cmaq_swdown_pm25_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_swdown_this,
			obs_pm25_lat_use, obs_pm25_lon_use, opt=2)
		cmaq_swdown_pm25_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_swdown_this,
			obs_pm25_lat_val, obs_pm25_lon_val, opt=2)

		print('Interpolating CMAQ PBL data to AirNow stations')
		cmaq_pbl_o3_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pbl_this,
			obs_o3_lat_use, obs_o3_lon_use, opt=2)
		cmaq_pbl_o3_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pbl_this,
			obs_o3_lat_val, obs_o3_lon_val, opt=2)
		cmaq_pbl_pm25_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pbl_this,
			obs_pm25_lat_use, obs_pm25_lon_use, opt=2)
		cmaq_pbl_pm25_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_pbl_this,
			obs_pm25_lat_val, obs_pm25_lon_val, opt=2)

		print('Interpolating CMAQ NOX data to AirNow stations')
		cmaq_nox_o3_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_nox_this,
			obs_o3_lat_use, obs_o3_lon_use, opt=2)
		cmaq_nox_o3_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_nox_this,
			obs_o3_lat_val, obs_o3_lon_val, opt=2)
		cmaq_nox_pm25_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_nox_this,
			obs_pm25_lat_use, obs_pm25_lon_use, opt=2)
		cmaq_nox_pm25_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_nox_this,
			obs_pm25_lat_val, obs_pm25_lon_val, opt=2)

		print('Interpolating CMAQ NOY data to AirNow stations')
		cmaq_noy_o3_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_noy_this,
			obs_o3_lat_use, obs_o3_lon_use, opt=2)
		cmaq_noy_o3_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_noy_this,
			obs_o3_lat_val, obs_o3_lon_val, opt=2)
		cmaq_noy_pm25_use_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_noy_this,
			obs_pm25_lat_use, obs_pm25_lon_use, opt=2)
		cmaq_noy_pm25_val_dummy = geocat.comp.rcm2points(cmaq_lat, cmaq_lon, cmaq_noy_this,
			obs_pm25_lat_val, obs_pm25_lon_val, opt=2)

		for ll in range(lead_hrs+1):
			cmaq_o3_con_use_stn[:,cc,ll] = cmaq_o3_con_use_dummy[ll,:]
			cmaq_o3_con_val_stn[:,cc,ll] = cmaq_o3_con_val_dummy[ll,:]
			cmaq_pm25_con_use_stn[:,cc,ll] = cmaq_pm25_con_use_dummy[ll,:]
			cmaq_pm25_con_val_stn[:,cc,ll] = cmaq_pm25_con_val_dummy[ll,:]

			cmaq_t2m_o3_use_stn[:,cc,ll] = cmaq_t2m_o3_use_dummy[ll,:]
			cmaq_t2m_o3_val_stn[:,cc,ll] = cmaq_t2m_o3_val_dummy[ll,:]
			cmaq_t2m_pm25_use_stn[:,cc,ll] = cmaq_t2m_pm25_use_dummy[ll,:]
			cmaq_t2m_pm25_val_stn[:,cc,ll] = cmaq_t2m_pm25_val_dummy[ll,:]

			cmaq_ws10m_o3_use_stn[:,cc,ll] = cmaq_ws10m_o3_use_dummy[ll,:]
			cmaq_ws10m_o3_val_stn[:,cc,ll] = cmaq_ws10m_o3_val_dummy[ll,:]
			cmaq_ws10m_pm25_use_stn[:,cc,ll] = cmaq_ws10m_pm25_use_dummy[ll,:]
			cmaq_ws10m_pm25_val_stn[:,cc,ll] = cmaq_ws10m_pm25_val_dummy[ll,:]

			cmaq_wd10m_o3_use_stn[:,cc,ll] = cmaq_wd10m_o3_use_dummy[ll,:]
			cmaq_wd10m_o3_val_stn[:,cc,ll] = cmaq_wd10m_o3_val_dummy[ll,:]
			cmaq_wd10m_pm25_use_stn[:,cc,ll] = cmaq_wd10m_pm25_use_dummy[ll,:]
			cmaq_wd10m_pm25_val_stn[:,cc,ll] = cmaq_wd10m_pm25_val_dummy[ll,:]

			cmaq_swdown_o3_use_stn[:,cc,ll] = cmaq_swdown_o3_use_dummy[ll,:]
			cmaq_swdown_o3_val_stn[:,cc,ll] = cmaq_swdown_o3_val_dummy[ll,:]
			cmaq_swdown_pm25_use_stn[:,cc,ll] = cmaq_swdown_pm25_use_dummy[ll,:]
			cmaq_swdown_pm25_val_stn[:,cc,ll] = cmaq_swdown_pm25_val_dummy[ll,:]

			cmaq_pbl_o3_use_stn[:,cc,ll] = cmaq_pbl_o3_use_dummy[ll,:]
			cmaq_pbl_o3_val_stn[:,cc,ll] = cmaq_pbl_o3_val_dummy[ll,:]
			cmaq_pbl_pm25_use_stn[:,cc,ll] = cmaq_pbl_pm25_use_dummy[ll,:]
			cmaq_pbl_pm25_val_stn[:,cc,ll] = cmaq_pbl_pm25_val_dummy[ll,:]

			cmaq_nox_o3_use_stn[:,cc,ll] = cmaq_nox_o3_use_dummy[ll,:]
			cmaq_nox_o3_val_stn[:,cc,ll] = cmaq_nox_o3_val_dummy[ll,:]
			cmaq_nox_pm25_use_stn[:,cc,ll] = cmaq_nox_pm25_use_dummy[ll,:]
			cmaq_nox_pm25_val_stn[:,cc,ll] = cmaq_nox_pm25_val_dummy[ll,:]

			cmaq_noy_o3_use_stn[:,cc,ll] = cmaq_noy_o3_use_dummy[ll,:]
			cmaq_noy_o3_val_stn[:,cc,ll] = cmaq_noy_o3_val_dummy[ll,:]
			cmaq_noy_pm25_use_stn[:,cc,ll] = cmaq_noy_pm25_use_dummy[ll,:]
			cmaq_noy_pm25_val_stn[:,cc,ll] = cmaq_noy_pm25_val_dummy[ll,:]
		'''

		## Using RectBivariateSpline for interpolation doesn't work on curvilinear grids like CAMS, WRF, etc.
		## Instead, do nearest-neighbor "interpolation" to AirNow stations
		## First, get and store the nearest-neighbor grid indices for the AirNow use & val stations for PM2.5 & O3
		## Only need to do this once
		if cc == 0:
			print('   Finding grid inds for o3_use')
			j_o3_use, i_o3_use = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_o3_lat_use, obs_o3_lon_use)
			print('   Finding grid inds for o3_val')
			j_o3_val, i_o3_val = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_o3_lat_val, obs_o3_lat_val)
			print('   Finding grid inds for pm25_use')
			j_pm25_use, i_pm25_use = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_pm25_lat_use, obs_pm25_lon_use)
			print('   Finding grid inds for pm25_val')
			j_pm25_val, i_pm25_val = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_pm25_lat_val, obs_pm25_lon_val)

		print('Interpolating nearest-neighbor CMAQ values to AirNow stations')
		for nn in range(n_obs_o3_use):
#			## Find closest j,i index
#			dist = abs(cmaq_lat-obs_o3_lat_use[nn])+abs(cmaq_lon-obs_o3_lon_use[nn])
#			j, i = np.unravel_index(dist.argmin(), dist.shape)
#			cmaq_o3_con_use_stn[nn,cc,:] = cmaq_o3_this[:,j,i]
#			cmaq_t2m_o3_use_stn[nn,cc,:] = cmaq_t2m_this[:,j,i]
#			cmaq_ws10m_o3_use_stn[nn,cc,:] = cmaq_ws10m_this[:,j,i]
#			cmaq_wd10m_o3_use_stn[nn,cc,:] = cmaq_wd10m_this[:,j,i]
#			cmaq_swdown_o3_use_stn[nn,cc,:] = cmaq_swdown_this[:,j,i]
#			cmaq_pbl_o3_use_stn[nn,cc,:] = cmaq_pbl_this[:,j,i]
#			cmaq_nox_o3_use_stn[nn,cc,:] = cmaq_nox_this[:,j,i]
#			cmaq_noy_o3_use_stn[nn,cc,:] = cmaq_noy_this[:,j,i]
			cmaq_o3_con_use_stn[nn,cc,:] = cmaq_o3_this[:, j_o3_use[nn], i_o3_use[nn]]
			cmaq_t2m_o3_use_stn[nn,cc,:] = cmaq_t2m_this[:, j_o3_use[nn], i_o3_use[nn]]
			cmaq_ws10m_o3_use_stn[nn,cc,:] = cmaq_ws10m_this[:, j_o3_use[nn], i_o3_use[nn]]
			cmaq_wd10m_o3_use_stn[nn,cc,:] = cmaq_wd10m_this[:, j_o3_use[nn], i_o3_use[nn]]
			cmaq_swdown_o3_use_stn[nn,cc,:] = cmaq_swdown_this[:, j_o3_use[nn], i_o3_use[nn]]
			cmaq_pbl_o3_use_stn[nn,cc,:] = cmaq_pbl_this[:, j_o3_use[nn], i_o3_use[nn]]
			cmaq_nox_o3_use_stn[nn,cc,:] = cmaq_nox_this[:, j_o3_use[nn], i_o3_use[nn]]
			cmaq_noy_o3_use_stn[nn,cc,:] = cmaq_noy_this[:, j_o3_use[nn], i_o3_use[nn]]

		for nn in range(n_obs_o3_val):
#			## Find closest j,i index
#			dist = abs(cmaq_lat-obs_o3_lat_val[nn])+abs(cmaq_lon-obs_o3_lon_val[nn])
#			j, i = np.unravel_index(dist.argmin(), dist.shape)
#			cmaq_o3_con_val_stn[nn,cc,:] = cmaq_o3_this[:,j,i]
#			cmaq_t2m_o3_val_stn[nn,cc,:] = cmaq_t2m_this[:,j,i]
#			cmaq_ws10m_o3_val_stn[nn,cc,:] = cmaq_ws10m_this[:,j,i]
#			cmaq_wd10m_o3_val_stn[nn,cc,:] = cmaq_wd10m_this[:,j,i]
#			cmaq_swdown_o3_val_stn[nn,cc,:] = cmaq_swdown_this[:,j,i]
#			cmaq_pbl_o3_val_stn[nn,cc,:] = cmaq_pbl_this[:,j,i]
#			cmaq_nox_o3_val_stn[nn,cc,:] = cmaq_nox_this[:,j,i]
#			cmaq_noy_o3_val_stn[nn,cc,:] = cmaq_noy_this[:,j,i]
			cmaq_o3_con_val_stn[nn,cc,:] = cmaq_o3_this[:, j_o3_val[nn], i_o3_val[nn]]
			cmaq_t2m_o3_val_stn[nn,cc,:] = cmaq_t2m_this[:, j_o3_val[nn], i_o3_val[nn]]
			cmaq_ws10m_o3_val_stn[nn,cc,:] = cmaq_ws10m_this[:, j_o3_val[nn], i_o3_val[nn]]
			cmaq_wd10m_o3_val_stn[nn,cc,:] = cmaq_wd10m_this[:, j_o3_val[nn], i_o3_val[nn]]
			cmaq_swdown_o3_val_stn[nn,cc,:] = cmaq_swdown_this[:, j_o3_val[nn], i_o3_val[nn]]
			cmaq_pbl_o3_val_stn[nn,cc,:] = cmaq_pbl_this[:, j_o3_val[nn], i_o3_val[nn]]
			cmaq_nox_o3_val_stn[nn,cc,:] = cmaq_nox_this[:, j_o3_val[nn], i_o3_val[nn]]
			cmaq_noy_o3_val_stn[nn,cc,:] = cmaq_noy_this[:, j_o3_val[nn], i_o3_val[nn]]

		for nn in range(n_obs_pm25_use):
#			## Find closest j,i index
#			dist = abs(cmaq_lat-obs_pm25_lat_use[nn])+abs(cmaq_lon-obs_pm25_lon_use[nn])
#			j, i = np.unravel_index(dist.argmin(), dist.shape)
#			cmaq_pm25_con_use_stn[nn,cc,:] = cmaq_pm25_this[:,j,i]
#			cmaq_t2m_pm25_use_stn[nn,cc,:] = cmaq_t2m_this[:,j,i]
#			cmaq_ws10m_pm25_use_stn[nn,cc,:] = cmaq_ws10m_this[:,j,i]
#			cmaq_wd10m_pm25_use_stn[nn,cc,:] = cmaq_wd10m_this[:,j,i]
#			cmaq_swdown_pm25_use_stn[nn,cc,:] = cmaq_swdown_this[:,j,i]
#			cmaq_pbl_pm25_use_stn[nn,cc,:] = cmaq_pbl_this[:,j,i]
#			cmaq_nox_pm25_use_stn[nn,cc,:] = cmaq_nox_this[:,j,i]
#			cmaq_noy_pm25_use_stn[nn,cc,:] = cmaq_noy_this[:,j,i]
			cmaq_pm25_con_use_stn[nn,cc,:] = cmaq_pm25_this[:, j_pm25_use[nn], i_pm25_use[nn]]
			cmaq_t2m_pm25_use_stn[nn,cc,:] = cmaq_t2m_this[:, j_pm25_use[nn], i_pm25_use[nn]]
			cmaq_ws10m_pm25_use_stn[nn,cc,:] = cmaq_ws10m_this[:, j_pm25_use[nn], i_pm25_use[nn]]
			cmaq_wd10m_pm25_use_stn[nn,cc,:] = cmaq_wd10m_this[:, j_pm25_use[nn], i_pm25_use[nn]]
			cmaq_swdown_pm25_use_stn[nn,cc,:] = cmaq_swdown_this[:, j_pm25_use[nn], i_pm25_use[nn]]
			cmaq_pbl_pm25_use_stn[nn,cc,:] = cmaq_pbl_this[:, j_pm25_use[nn], i_pm25_use[nn]]
			cmaq_nox_pm25_use_stn[nn,cc,:] = cmaq_nox_this[:, j_pm25_use[nn], i_pm25_use[nn]]
			cmaq_noy_pm25_use_stn[nn,cc,:] = cmaq_noy_this[:, j_pm25_use[nn], i_pm25_use[nn]]

		for nn in range(n_obs_pm25_val):
#			## Find closest j,i index
#			dist = abs(cmaq_lat-obs_pm25_lat_val[nn])+abs(cmaq_lon-obs_pm25_lon_val[nn])
#			j, i = np.unravel_index(dist.argmin(), dist.shape)
#			cmaq_pm25_con_val_stn[nn,cc,:] = cmaq_pm25_this[:,j,i]
#			cmaq_t2m_pm25_val_stn[nn,cc,:] = cmaq_t2m_this[:,j,i]
#			cmaq_ws10m_pm25_val_stn[nn,cc,:] = cmaq_ws10m_this[:,j,i]
#			cmaq_wd10m_pm25_val_stn[nn,cc,:] = cmaq_wd10m_this[:,j,i]
#			cmaq_swdown_pm25_val_stn[nn,cc,:] = cmaq_swdown_this[:,j,i]
#			cmaq_pbl_pm25_val_stn[nn,cc,:] = cmaq_pbl_this[:,j,i]
#			cmaq_nox_pm25_val_stn[nn,cc,:] = cmaq_nox_this[:,j,i]
#			cmaq_noy_pm25_val_stn[nn,cc,:] = cmaq_noy_this[:,j,i]
			cmaq_pm25_con_val_stn[nn,cc,:] = cmaq_pm25_this[:, j_pm25_val[nn], i_pm25_val[nn]]
			cmaq_t2m_pm25_val_stn[nn,cc,:] = cmaq_t2m_this[:, j_pm25_val[nn], i_pm25_val[nn]]
			cmaq_ws10m_pm25_val_stn[nn,cc,:] = cmaq_ws10m_this[:, j_pm25_val[nn], i_pm25_val[nn]]
			cmaq_wd10m_pm25_val_stn[nn,cc,:] = cmaq_wd10m_this[:, j_pm25_val[nn], i_pm25_val[nn]]
			cmaq_swdown_pm25_val_stn[nn,cc,:] = cmaq_swdown_this[:, j_pm25_val[nn], i_pm25_val[nn]]
			cmaq_pbl_pm25_val_stn[nn,cc,:] = cmaq_pbl_this[:, j_pm25_val[nn], i_pm25_val[nn]]
			cmaq_nox_pm25_val_stn[nn,cc,:] = cmaq_nox_this[:, j_pm25_val[nn], i_pm25_val[nn]]
			cmaq_noy_pm25_val_stn[nn,cc,:] = cmaq_noy_this[:, j_pm25_val[nn], i_pm25_val[nn]]
		
		## Create xarray data arrays using the CMAQ and obs data
		n_use_o3 = len(obs_o3_lat_use)
		n_val_o3 = len(obs_o3_lat_val)
		n_use_pm25 = len(obs_pm25_lat_use)
		n_val_pm25 = len(obs_pm25_lat_val)
		n_leads = lead_hrs+1

		cmaq_str = 'CMAQ raw gridded forecast'

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

		cmaq_t2m_o3_use = xr.DataArray(cmaq_t2m_o3_use_stn[:,cc,:],
									coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
									attrs={'description':cmaq_str+' 2-m temperature interpolated to AirNow O3 stations used for merging and training', 'units':'K'})
		cmaq_t2m_o3_val = xr.DataArray(cmaq_t2m_o3_val_stn[:,cc,:],
									coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
									attrs={'description':cmaq_str+' 2-m temperature interpolated to AirNow O3 stations withheld for validation', 'units':'K'})
		cmaq_t2m_pm25_use = xr.DataArray(cmaq_t2m_pm25_use_stn[:,cc,:],
									coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
									attrs={'description':cmaq_str+' 2-m temperature interpolated to AirNow PM2.5 stations used for merging and training', 'units':'K'})
		cmaq_t2m_pm25_val = xr.DataArray(cmaq_t2m_pm25_val_stn[:,cc,:],
									coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
									attrs={'description':cmaq_str+' 2-m temperature interpolated to AirNow PM2.5 stations withheld for validation', 'units':'K'})


		cmaq_ws10m_o3_use = xr.DataArray(cmaq_ws10m_o3_use_stn[:,cc,:],
									coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
									attrs={'description':cmaq_str+' 10-m wind speed interpolated to AirNow O3 stations used for merging and training', 'units':'m s-1'})
		cmaq_ws10m_o3_val = xr.DataArray(cmaq_ws10m_o3_val_stn[:,cc,:],
									coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
									attrs={'description':cmaq_str+' 10-m wind speed interpolated to AirNow O3 stations withheld for validation', 'units':'m s-1'})
		cmaq_ws10m_pm25_use = xr.DataArray(cmaq_ws10m_pm25_use_stn[:,cc,:],
									coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
									attrs={'description':cmaq_str+' 10-m wind speed interpolated to AirNow PM2.5 stations used for merging and training', 'units':'m s-1'})
		cmaq_ws10m_pm25_val = xr.DataArray(cmaq_ws10m_pm25_val_stn[:,cc,:],
									coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
									attrs={'description':cmaq_str+' 10-m wind speed interpolated to AirNow PM2.5 stations withheld for validation', 'units':'m s-1'})

		cmaq_wd10m_o3_use = xr.DataArray(cmaq_wd10m_o3_use_stn[:,cc,:],
									coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
									attrs={'description':cmaq_str+' 10-m wind direction interpolated to AirNow O3 stations used for merging and training', 'units':'deg'})
		cmaq_wd10m_o3_val = xr.DataArray(cmaq_wd10m_o3_val_stn[:,cc,:],
									coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
									attrs={'description':cmaq_str+' 10-m wind direction interpolated to AirNow O3 stations withheld for validation', 'units':'deg'})
		cmaq_wd10m_pm25_use = xr.DataArray(cmaq_wd10m_pm25_use_stn[:,cc,:],
									coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
									attrs={'description':cmaq_str+' 10-m wind direction interpolated to AirNow PM2.5 stations used for merging and training', 'units':'deg'})
		cmaq_wd10m_pm25_val = xr.DataArray(cmaq_wd10m_pm25_val_stn[:,cc,:],
									coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
									attrs={'description':cmaq_str+' 10-m wind direction interpolated to AirNow PM2.5 stations withheld for validation', 'units':'deg'})

		cmaq_swdown_o3_use = xr.DataArray(cmaq_swdown_o3_use_stn[:,cc,:],
									coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
									attrs={'description':cmaq_str+' GHI interpolated to AirNow O3 stations used for merging and training', 'units':'W m-2'})
		cmaq_swdown_o3_val = xr.DataArray(cmaq_swdown_o3_val_stn[:,cc,:],
									coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
									attrs={'description':cmaq_str+' GHI interpolated to AirNow O3 stations withheld for validation', 'units':'W m-2'})
		cmaq_swdown_pm25_use = xr.DataArray(cmaq_swdown_pm25_use_stn[:,cc,:],
									coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
									attrs={'description':cmaq_str+' GHI interpolated to AirNow PM2.5 stations used for merging and training', 'units':'W m-2'})
		cmaq_swdown_pm25_val = xr.DataArray(cmaq_swdown_pm25_val_stn[:,cc,:],
									coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
									attrs={'description':cmaq_str+' GHI interpolated to AirNow PM2.5 stations withheld for validation', 'units':'W m-2'})

		cmaq_pbl_o3_use = xr.DataArray(cmaq_pbl_o3_use_stn[:,cc,:],
									coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
									attrs={'description':cmaq_str+' PBL depth interpolated to AirNow O3 stations used for merging and training', 'units':'m'})
		cmaq_pbl_o3_val = xr.DataArray(cmaq_pbl_o3_val_stn[:,cc,:],
									coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
									attrs={'description':cmaq_str+' PBL depth interpolated to AirNow O3 stations withheld for validation', 'units':'m'})
		cmaq_pbl_pm25_use = xr.DataArray(cmaq_pbl_pm25_use_stn[:,cc,:],
									coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
									attrs={'description':cmaq_str+' PBL depth interpolated to AirNow PM2.5 stations used for merging and training', 'units':'m'})
		cmaq_pbl_pm25_val = xr.DataArray(cmaq_pbl_pm25_val_stn[:,cc,:],
									coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
									attrs={'description':cmaq_str+' PBL depth interpolated to AirNow PM2.5 stations withheld for validation', 'units':'m'})

		cmaq_nox_o3_use = xr.DataArray(cmaq_nox_o3_use_stn[:,cc,:],
									coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
									attrs={'description':cmaq_str+' NOx interpolated to AirNow O3 stations used for merging and training', 'units':'ppmv'})
		cmaq_nox_o3_val = xr.DataArray(cmaq_nox_o3_val_stn[:,cc,:],
									coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
									attrs={'description':cmaq_str+' NOx interpolated to AirNow O3 stations withheld for validation', 'units':'ppmv'})
		cmaq_nox_pm25_use = xr.DataArray(cmaq_nox_pm25_use_stn[:,cc,:],
									coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
									attrs={'description':cmaq_str+' NOx interpolated to AirNow PM2.5 stations used for merging and training', 'units':'ppmv'})
		cmaq_nox_pm25_val = xr.DataArray(cmaq_nox_pm25_val_stn[:,cc,:],
									coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
									attrs={'description':cmaq_str+' NOx interpolated to AirNow PM2.5 stations withheld for validation', 'units':'ppmv'})

		cmaq_noy_o3_use = xr.DataArray(cmaq_noy_o3_use_stn[:,cc,:],
									coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
									attrs={'description':cmaq_str+' NOy interpolated to AirNow O3 stations used for merging and training', 'units':'ppmv'})
		cmaq_noy_o3_val = xr.DataArray(cmaq_noy_o3_val_stn[:,cc,:],
									coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
									attrs={'description':cmaq_str+' NOy interpolated to AirNow O3 stations withheld for validation', 'units':'ppmv'})
		cmaq_noy_pm25_use = xr.DataArray(cmaq_noy_pm25_use_stn[:,cc,:],
									coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
									attrs={'description':cmaq_str+' NOy interpolated to AirNow PM2.5 stations used for merging and training', 'units':'ppmv'})
		cmaq_noy_pm25_val = xr.DataArray(cmaq_noy_pm25_val_stn[:,cc,:],
									coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
									attrs={'description':cmaq_str+' NOy interpolated to AirNow PM2.5 stations withheld for validation', 'units':'ppmv'})


		## Create new xarray datasets
		ds_use = xr.Dataset(
						data_vars={ 'cmaq_pm25_con_use':cmaq_pm25_con_use, 'cmaq_o3_con_use':cmaq_o3_con_use,
										'cmaq_t2m_pm25_use':cmaq_t2m_pm25_use, 'cmaq_t2m_o3_use':cmaq_t2m_o3_use,
										'cmaq_ws10m_pm25_use':cmaq_ws10m_pm25_use, 'cmaq_ws10m_o3_use':cmaq_ws10m_o3_use,
										'cmaq_wd10m_pm25_use':cmaq_wd10m_pm25_use, 'cmaq_wd10m_o3_use':cmaq_wd10m_o3_use,
										'cmaq_swdown_pm25_use':cmaq_swdown_pm25_use, 'cmaq_swdown_o3_use':cmaq_swdown_o3_use,
										'cmaq_pbl_pm25_use':cmaq_pbl_pm25_use, 'cmaq_pbl_o3_use':cmaq_pbl_o3_use,
										'cmaq_nox_pm25_use':cmaq_nox_pm25_use, 'cmaq_nox_o3_use':cmaq_nox_o3_use,
										'cmaq_noy_pm25_use':cmaq_noy_pm25_use, 'cmaq_noy_o3_use':cmaq_noy_o3_use,
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
										'cmaq_t2m_pm25_val':cmaq_t2m_pm25_val, 'cmaq_t2m_o3_val':cmaq_t2m_o3_val,
										'cmaq_ws10m_pm25_val':cmaq_ws10m_pm25_val, 'cmaq_ws10m_o3_val':cmaq_ws10m_o3_val,
										'cmaq_wd10m_pm25_val':cmaq_wd10m_pm25_val, 'cmaq_wd10m_o3_val':cmaq_wd10m_o3_val,
										'cmaq_swdown_pm25_val':cmaq_swdown_pm25_val, 'cmaq_swdown_o3_val':cmaq_swdown_o3_val,
										'cmaq_pbl_pm25_val':cmaq_pbl_pm25_val, 'cmaq_pbl_o3_val':cmaq_pbl_o3_val,
										'cmaq_nox_pm25_val':cmaq_nox_pm25_val, 'cmaq_nox_o3_val':cmaq_nox_o3_val,
										'cmaq_noy_pm25_val':cmaq_noy_pm25_val, 'cmaq_noy_o3_val':cmaq_noy_o3_val,
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

		## Set the output paths & filenames
		out_dir = out_dir_parent.joinpath(this_cycle_yr,this_cycle_mo)
		out_dir.mkdir(parents=True, exist_ok=True)
		out_dir.chmod(0o755)
		fname_use = out_dir.joinpath('cmaq_airnow_pm2.5_o3_'+this_cycle_outfile+'_use.nc')
		fname_val = out_dir.joinpath('cmaq_airnow_pm2.5_o3_'+this_cycle_outfile+'_val.nc')

#		if sites_vary:
#			fname_use = '/glade/scratch/jkim/CMAQ/airnow_cmaq_use_coll/cmaq_airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_use.nc'
#			fname_val = '/glade/scratch/jkim/CMAQ/airnow_cmaq_use_coll/cmaq_airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_val.nc'
#		if sites_static:
#			fname_use = '/glade/scratch/jkim/CMAQ/airnow_cmaq_use_coll/cmaq_airnow_pm2.5_o3_static_'+this_yyyymmdd_hh+'00_use.nc'
#			fname_val = '/glade/scratch/jkim/CMAQ/airnow_cmaq_use_coll/cmaq_airnow_pm2.5_o3_static_'+this_yyyymmdd_hh+'00_val.nc'
#			fname_use = out_dir.joinpath('cmaq_airnow_pm2.5_o3_'+this_cycle_outfile+'_use.nc'
#			fname_val = out_dir.joinpath('cmaq_airnow_pm2.5_o3_'+this_cycle_outfile+'_val.nc'

		## Write the datasets to NetCDF
		print('Writing '+str(fname_use))
		ds_use.to_netcdf(fname_use)
		fname_use.chmod(0o644)
		print('Writing '+str(fname_val))
		ds_val.to_netcdf(fname_val)
		fname_val.chmod(0o644)


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
