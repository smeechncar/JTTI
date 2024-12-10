'''
extract_cmaq_raw_bm_cams_bm_at_airnow_use_val_leadtimeloop.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 15 Mar 2023

This script:
1. Reads in AirNow obs NetCDF files that have already been separated into use & val sets, for use
in merging & validation, respectively.
2. Reads in a CMAQ sample gridded file to get the grid coordinates of those AirNow stations and
finds the grid indices of the nearest-neighbor grid points to each station.
3. Reads in bias-corrected CMAQ PM2.5 & O3 files that are produced hourly to serve as gridded
"analysis"/"truth".
4. Reads in bias-corrected CAMS PM2.5 & O3 files that are produced hourly and have been regridded
to the CMAQ grid.
5. Reads in raw CMAQ forecast output to get the 9 variables needed for AnEn predictors for all 72
hours in each 06z forecast cycle.
6. Writes NetCDF output files with AirNow obs, BC CMAQ-interpolated values, BC CAMS-interpolated values,
and raw CMAQ variables interpolated to AirNow sites for each forecast cycle.
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
	parser.add_argument('-f', '--lead_hrs', type=int, default=72, help='forecast lead hours [default: 72, for hours 0-72]')

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

	## All 9 AnEn predictors
	cmaq_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','06z_cmaq_predictors')
	## PM2.5 and O3
	cmaq_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','BM3','static','cmaq_include0_nn')
	## PM2.5 and O3
	cams_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','BM3','static','cams_include0')

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
	dt_end_valid = dt.datetime.strptime(date_end, fmt_yyyymmdd_hh) + dt.timedelta(hours=lead_hrs)
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

	## 1. Loop over valid times to get the observations as a time series
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

	## 2. Open a CMAQ sample/coordinate file to get grid info. Only need lat/lon, not map projection info.
	cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
	print('Reading CMAQ coordinate data from '+str(cmaq_fname))
	cmaq_ds = xr.open_dataset(cmaq_fname)
	cmaq_lon = cmaq_ds.LON[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})
	cmaq_lat = cmaq_ds.LAT[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})

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

	## Create valid-time arrays to hold time series of CMAQ BM & CAMS BM PM2.5 & O3
	cmaq_bm3_o3_con_use_ts = np.full([n_obs_o3_use, n_valid], np.nan)
	cmaq_bm3_o3_con_val_ts = np.full([n_obs_o3_val, n_valid], np.nan)
	cmaq_bm3_pm25_con_use_ts = np.full([n_obs_pm25_use, n_valid], np.nan)
	cmaq_bm3_pm25_con_val_ts = np.full([n_obs_pm25_val, n_valid], np.nan)

	cams_bm3_o3_con_use_ts = np.full([n_obs_o3_use, n_valid], np.nan)
	cams_bm3_o3_con_val_ts = np.full([n_obs_o3_val, n_valid], np.nan)
	cams_bm3_pm25_con_use_ts = np.full([n_obs_pm25_use, n_valid], np.nan)
	cams_bm3_pm25_con_val_ts = np.full([n_obs_pm25_val, n_valid], np.nan)

	n_use_o3 = len(obs_o3_lat_use)
	n_val_o3 = len(obs_o3_lat_val)
	n_use_pm25 = len(obs_pm25_lat_use)
	n_val_pm25 = len(obs_pm25_lat_val)
	n_leads = lead_hrs+1

	## Loop over valid times
	for vv in range(n_valid):
		this_dt = dt_arr_valid[vv]
		this_yyyymmdd_hh = this_dt.strftime(fmt_yyyymmdd_hh)
		this_yyyymmddhh  = this_dt.strftime(fmt_yyyymmddhh)
		this_yyyymmdd = this_dt.strftime(fmt_yyyymmdd)
		this_yr = this_dt.strftime(fmt_yyyy)
		this_mo = this_dt.strftime(fmt_mm)
		this_hr = this_dt.strftime(fmt_hh)
		this_date_hh = this_dt.strftime(fmt_date_hh)

		## 3. Read in CMAQ BM3 "truth" files (Ju-Hye has all CMAQ BM3 hourly files in one directory)
		cmaq_fname_in_this_o3   = cmaq_bm3_dir.joinpath('cmaq_o3_'+this_yyyymmdd_hh+'00_BM3_static_inc0_nn.nc')
		cmaq_fname_in_this_pm25 = cmaq_bm3_dir.joinpath('cmaq_pm25_'+this_yyyymmdd_hh+'00_BM3_static_noQC.nc')

		if cmaq_fname_in_this_o3.is_file():
			print('Reading '+str(cmaq_fname_in_this_o3))
			ds_cmaq_o3 = xr.open_dataset(cmaq_fname_in_this_o3)
			cmaq_o3_this = ds_cmaq_o3.cmaq_o3_m
			## At the first time all values will be 0.0. Make it nan so it doesn't accidentally get used downstream.
			if np.all(cmaq_o3_this == 0.0):
				cmaq_o3_this[:,:] = np.nan

			## Extract nearest-neighbor values at AirNow stations
			for nn in range(n_use_o3):
				cmaq_bm3_o3_con_use_ts[nn,vv] = cmaq_o3_this[j_o3_use[nn], i_o3_use[nn]]
			for nn in range(n_val_o3):
				cmaq_bm3_o3_con_val_ts[nn,vv] = cmaq_o3_this[j_o3_val[nn], i_o3_val[nn]]
		else:
			print('WARNING: File '+str(cmaq_fname_in_this_o3)+' does not exist!')

		if cmaq_fname_in_this_pm25.is_file():
			print('Reading '+str(cmaq_fname_in_this_pm25))
			ds_cmaq_pm25 = xr.open_dataset(cmaq_fname_in_this_pm25)
			cmaq_pm25_this = ds_cmaq_pm25.cmaq_pm25_m
			if np.all(cmaq_pm25_this == 0.0):
				cmaq_pm25_this[:,:] = np.nan

			## Extract nearest-neighbor values at AirNow stations
			for nn in range(n_use_pm25):
				cmaq_bm3_pm25_con_use_ts[nn,vv] = cmaq_pm25_this[j_pm25_use[nn], i_pm25_use[nn]]
			for nn in range(n_val_pm25):
				cmaq_bm3_pm25_con_val_ts[nn,vv] = cmaq_pm25_this[j_pm25_val[nn], i_pm25_val[nn]]
		else:
			print('WARNING: File '+str(cmaq_fname_in_this_pm25)+' does not exist!')

		## 4. Read in CAMS BM3 "truth" files (Ju-Hye has all CAMS BM files in one directory as well)
		cams_fname_in_this_o3   = cams_bm3_dir.joinpath('cams_o3_regrid_cmaq_'+this_yyyymmdd_hh+'00_BM3_static_include0.nc')
		cams_fname_in_this_pm25 = cams_bm3_dir.joinpath('cams_pm25_regrid_cmaq_'+this_yyyymmdd_hh+'00_BM3_static_noQC.nc')

		if cams_fname_in_this_o3.is_file():
			print('Reading '+str(cams_fname_in_this_o3))
			ds_cams_o3 = xr.open_dataset(cams_fname_in_this_o3)
			cams_o3_this = ds_cams_o3.cams_o3_m
			if np.all(cams_o3_this == 0.0):
				cams_o3_this[:,:] = np.nan

			## Extract nearest-neighbor values at AirNow stations
			for nn in range(n_use_o3):
				cams_bm3_o3_con_use_ts[nn,vv] = cams_o3_this[j_o3_use[nn], i_o3_use[nn]]
			for nn in range(n_val_o3):
				cams_bm3_o3_con_val_ts[nn,vv] = cams_o3_this[j_o3_val[nn], i_o3_val[nn]]
		else:
			print('WARNING: File '+str(cams_fname_in_this_o3)+' does not exist!')

		if cams_fname_in_this_pm25.is_file():
			print('Reading '+str(cams_fname_in_this_pm25))
			ds_cams_pm25 = xr.open_dataset(cams_fname_in_this_pm25)
			cams_pm25_this = ds_cams_pm25.cams_pm25_m
			if np.all(cams_pm25_this == 0.0):
				cams_pm25_this[:,:] = np.nan

			## Extract nearest-neighbor values at AirNow stations
			for nn in range(n_use_pm25):
				cams_bm3_pm25_con_use_ts[nn,vv] = cams_pm25_this[j_pm25_use[nn], i_pm25_use[nn]]
			for nn in range(n_val_pm25):
				cams_bm3_pm25_con_val_ts[nn,vv] = cams_pm25_this[j_pm25_val[nn], i_pm25_val[nn]]
		else:
			print('WARNING: File '+str(cams_fname_in_this_pm25)+' does not exist!')

	## Create 2D arrays to hold all forecast-observation pairs
	obs_pm25_con_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	obs_pm25_con_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	obs_o3_con_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	obs_o3_con_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_bm3_pm25_con_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_bm3_pm25_con_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_bm3_o3_con_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_bm3_o3_con_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cams_bm3_pm25_con_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cams_bm3_pm25_con_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cams_bm3_o3_con_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cams_bm3_o3_con_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_raw_pm25_con_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_pm25_con_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_o3_con_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_o3_con_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_raw_t2m_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_t2m_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_t2m_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_t2m_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_raw_ws10m_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_ws10m_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_ws10m_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_ws10m_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_raw_wd10m_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_wd10m_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_wd10m_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_wd10m_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_raw_swdown_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_swdown_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_swdown_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_swdown_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_raw_pbl_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_pbl_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_pbl_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_pbl_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_raw_nox_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_nox_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_nox_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_nox_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	cmaq_raw_noy_pm25_use_stn = np.full([n_obs_pm25_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_noy_pm25_val_stn = np.full([n_obs_pm25_val, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_noy_o3_use_stn = np.full([n_obs_o3_use, n_cycle, lead_hrs+1], np.nan)
	cmaq_raw_noy_o3_val_stn = np.full([n_obs_o3_val, n_cycle, lead_hrs+1], np.nan)

	## Put the time series obs, CMAQ BM, and CAMS BM into cycle-lead time arrays
	print('Putting the time series obs, CMAQ BM, and CAMS BM into cycle-lead time arrays')
	for cc in range(n_cycle):
		for ll in range(lead_hrs+1):
			ind_valid = np.where(dt_arr_valid == dt_arr_pairs[cc,ll])[0][0]
			print('cc = '+str(cc)+', ll = '+str(ll)+', ind_valid = '+str(ind_valid))
			obs_pm25_con_use_stn[:,cc,ll] = obs_pm25_con_use_ts[:,ind_valid]
			obs_pm25_con_val_stn[:,cc,ll] = obs_pm25_con_val_ts[:,ind_valid]
			obs_o3_con_use_stn[:,cc,ll] = obs_o3_con_use_ts[:,ind_valid]
			obs_o3_con_val_stn[:,cc,ll] = obs_o3_con_val_ts[:,ind_valid]

			cmaq_bm3_pm25_con_use_stn[:,cc,ll] = cmaq_bm3_pm25_con_use_ts[:,ind_valid]
			cmaq_bm3_pm25_con_val_stn[:,cc,ll] = cmaq_bm3_pm25_con_val_ts[:,ind_valid]
			cmaq_bm3_o3_con_use_stn[:,cc,ll] = cmaq_bm3_o3_con_use_ts[:,ind_valid]
			cmaq_bm3_o3_con_val_stn[:,cc,ll] = cmaq_bm3_o3_con_val_ts[:,ind_valid]

			cams_bm3_pm25_con_use_stn[:,cc,ll] = cams_bm3_pm25_con_use_ts[:,ind_valid]
			cams_bm3_pm25_con_val_stn[:,cc,ll] = cams_bm3_pm25_con_val_ts[:,ind_valid]
			cams_bm3_o3_con_use_stn[:,cc,ll] = cams_bm3_o3_con_use_ts[:,ind_valid]
			cams_bm3_o3_con_val_stn[:,cc,ll] = cams_bm3_o3_con_val_ts[:,ind_valid]

	## Loop over CMAQ raw forecast cycles
	for cc in range(n_cycle):
		this_cycle_dt = dt_arr_cycle[cc]
		this_cycle_infile = this_cycle_dt.strftime(fmt_yyyymmdd+'.'+fmt_hh+'z')
		this_cycle_outfile = this_cycle_dt.strftime(fmt_yyyymmdd_hhmm)
		this_cycle_yr = this_cycle_dt.strftime(fmt_yyyy)
		this_cycle_mo = this_cycle_dt.strftime(fmt_mm)

		## 5. Read in CMAQ files (Ju-Hye has all CMAQ hourly files in one directory)
		cmaq_fname_in_this = cmaq_raw_dir.joinpath('cmaq_9variables_'+this_cycle_infile+'.nc')

		if cmaq_fname_in_this.is_file():
			print('Reading '+str(cmaq_fname_in_this))
			ds_cmaq = xr.open_dataset(cmaq_fname_in_this)
			## Dimensions (lead_hours 73, latitude 265, longitude 442) 

#			if get_cmaq_lat_lon:
#				cmaq_lat = ds_cmaq.latitude
#				cmaq_lon = ds_cmaq.longitude
#				get_cmaq_lat_lon = False

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

		## Extract nearest-neighbor values at AirNow stations
		print('Interpolating nearest-neighbor CMAQ raw values to AirNow stations')
		for nn in range(n_obs_o3_use):
			cmaq_raw_o3_con_use_stn[nn,cc,:] = cmaq_o3_this[:,j_o3_use[nn],i_o3_use[nn]]
			cmaq_raw_t2m_o3_use_stn[nn,cc,:] = cmaq_t2m_this[:,j_o3_use[nn],i_o3_use[nn]]
			cmaq_raw_ws10m_o3_use_stn[nn,cc,:] = cmaq_ws10m_this[:,j_o3_use[nn],i_o3_use[nn]]
			cmaq_raw_wd10m_o3_use_stn[nn,cc,:] = cmaq_wd10m_this[:,j_o3_use[nn],i_o3_use[nn]]
			cmaq_raw_swdown_o3_use_stn[nn,cc,:] = cmaq_swdown_this[:,j_o3_use[nn],i_o3_use[nn]]
			cmaq_raw_pbl_o3_use_stn[nn,cc,:] = cmaq_pbl_this[:,j_o3_use[nn],i_o3_use[nn]]
			cmaq_raw_nox_o3_use_stn[nn,cc,:] = cmaq_nox_this[:,j_o3_use[nn],i_o3_use[nn]]
			cmaq_raw_noy_o3_use_stn[nn,cc,:] = cmaq_noy_this[:,j_o3_use[nn],i_o3_use[nn]]

		for nn in range(n_obs_o3_val):
			cmaq_raw_o3_con_val_stn[nn,cc,:] = cmaq_o3_this[:,j_o3_val[nn],i_o3_val[nn]]
			cmaq_raw_t2m_o3_val_stn[nn,cc,:] = cmaq_t2m_this[:,j_o3_val[nn],i_o3_val[nn]]
			cmaq_raw_ws10m_o3_val_stn[nn,cc,:] = cmaq_ws10m_this[:,j_o3_val[nn],i_o3_val[nn]]
			cmaq_raw_wd10m_o3_val_stn[nn,cc,:] = cmaq_wd10m_this[:,j_o3_val[nn],i_o3_val[nn]]
			cmaq_raw_swdown_o3_val_stn[nn,cc,:] = cmaq_swdown_this[:,j_o3_val[nn],i_o3_val[nn]]
			cmaq_raw_pbl_o3_val_stn[nn,cc,:] = cmaq_pbl_this[:,j_o3_val[nn],i_o3_val[nn]]
			cmaq_raw_nox_o3_val_stn[nn,cc,:] = cmaq_nox_this[:,j_o3_val[nn],i_o3_val[nn]]
			cmaq_raw_noy_o3_val_stn[nn,cc,:] = cmaq_noy_this[:,j_o3_val[nn],i_o3_val[nn]]

		for nn in range(n_obs_pm25_use):
			cmaq_raw_pm25_con_use_stn[nn,cc,:] = cmaq_pm25_this[:,j_pm25_use[nn],i_pm25_use[nn]]
			cmaq_raw_t2m_pm25_use_stn[nn,cc,:] = cmaq_t2m_this[:,j_pm25_use[nn],i_pm25_use[nn]]
			cmaq_raw_ws10m_pm25_use_stn[nn,cc,:] = cmaq_ws10m_this[:,j_pm25_use[nn],i_pm25_use[nn]]
			cmaq_raw_wd10m_pm25_use_stn[nn,cc,:] = cmaq_wd10m_this[:,j_pm25_use[nn],i_pm25_use[nn]]
			cmaq_raw_swdown_pm25_use_stn[nn,cc,:] = cmaq_swdown_this[:,j_pm25_use[nn],i_pm25_use[nn]]
			cmaq_raw_pbl_pm25_use_stn[nn,cc,:] = cmaq_pbl_this[:,j_pm25_use[nn],i_pm25_use[nn]]
			cmaq_raw_nox_pm25_use_stn[nn,cc,:] = cmaq_nox_this[:,j_pm25_use[nn],i_pm25_use[nn]]
			cmaq_raw_noy_pm25_use_stn[nn,cc,:] = cmaq_noy_this[:,j_pm25_use[nn],i_pm25_use[nn]]

		for nn in range(n_obs_pm25_val):
			cmaq_raw_pm25_con_val_stn[nn,cc,:] = cmaq_pm25_this[:,j_pm25_val[nn],i_pm25_val[nn]]
			cmaq_raw_t2m_pm25_val_stn[nn,cc,:] = cmaq_t2m_this[:,j_pm25_val[nn],i_pm25_val[nn]]
			cmaq_raw_ws10m_pm25_val_stn[nn,cc,:] = cmaq_ws10m_this[:,j_pm25_val[nn],i_pm25_val[nn]]
			cmaq_raw_wd10m_pm25_val_stn[nn,cc,:] = cmaq_wd10m_this[:,j_pm25_val[nn],i_pm25_val[nn]]
			cmaq_raw_swdown_pm25_val_stn[nn,cc,:] = cmaq_swdown_this[:,j_pm25_val[nn],i_pm25_val[nn]]
			cmaq_raw_pbl_pm25_val_stn[nn,cc,:] = cmaq_pbl_this[:,j_pm25_val[nn],i_pm25_val[nn]]
			cmaq_raw_nox_pm25_val_stn[nn,cc,:] = cmaq_nox_this[:,j_pm25_val[nn],i_pm25_val[nn]]
			cmaq_raw_noy_pm25_val_stn[nn,cc,:] = cmaq_noy_this[:,j_pm25_val[nn],i_pm25_val[nn]]

		'''
		## Using RectBivariateSpline for interpolation doesn't work on curvilinear grids like CAMS, WRF, etc.
		## Instead, do nearest-neighbor "interpolation" to AirNow stations
		print('Interpolating nearest-neighbor CMAQ values to AirNow stations')
		for nn in range(n_obs_o3_use):
			## Find closest j,i index
			dist = abs(cmaq_lat-obs_o3_lat_use[nn])+abs(cmaq_lon-obs_o3_lon_use[nn])
			j, i = np.unravel_index(dist.argmin(), dist.shape)
			cmaq_o3_con_use_stn[nn,cc,:] = cmaq_o3_this[:,j,i]
			cmaq_t2m_o3_use_stn[nn,cc,:] = cmaq_t2m_this[:,j,i]
			cmaq_ws10m_o3_use_stn[nn,cc,:] = cmaq_ws10m_this[:,j,i]
			cmaq_wd10m_o3_use_stn[nn,cc,:] = cmaq_wd10m_this[:,j,i]
			cmaq_swdown_o3_use_stn[nn,cc,:] = cmaq_swdown_this[:,j,i]
			cmaq_pbl_o3_use_stn[nn,cc,:] = cmaq_pbl_this[:,j,i]
			cmaq_nox_o3_use_stn[nn,cc,:] = cmaq_nox_this[:,j,i]
			cmaq_noy_o3_use_stn[nn,cc,:] = cmaq_noy_this[:,j,i]

		for nn in range(n_obs_o3_val):
			## Find closest j,i index
			dist = abs(cmaq_lat-obs_o3_lat_val[nn])+abs(cmaq_lon-obs_o3_lon_val[nn])
			j, i = np.unravel_index(dist.argmin(), dist.shape)
			cmaq_o3_con_val_stn[nn,cc,:] = cmaq_o3_this[:,j,i]
			cmaq_t2m_o3_val_stn[nn,cc,:] = cmaq_t2m_this[:,j,i]
			cmaq_ws10m_o3_val_stn[nn,cc,:] = cmaq_ws10m_this[:,j,i]
			cmaq_wd10m_o3_val_stn[nn,cc,:] = cmaq_wd10m_this[:,j,i]
			cmaq_swdown_o3_val_stn[nn,cc,:] = cmaq_swdown_this[:,j,i]
			cmaq_pbl_o3_val_stn[nn,cc,:] = cmaq_pbl_this[:,j,i]
			cmaq_nox_o3_val_stn[nn,cc,:] = cmaq_nox_this[:,j,i]
			cmaq_noy_o3_val_stn[nn,cc,:] = cmaq_noy_this[:,j,i]

		for nn in range(n_obs_pm25_use):
			## Find closest j,i index
			dist = abs(cmaq_lat-obs_pm25_lat_use[nn])+abs(cmaq_lon-obs_pm25_lon_use[nn])
			j, i = np.unravel_index(dist.argmin(), dist.shape)
			cmaq_pm25_con_use_stn[nn,cc,:] = cmaq_pm25_this[:,j,i]
			cmaq_t2m_pm25_use_stn[nn,cc,:] = cmaq_t2m_this[:,j,i]
			cmaq_ws10m_pm25_use_stn[nn,cc,:] = cmaq_ws10m_this[:,j,i]
			cmaq_wd10m_pm25_use_stn[nn,cc,:] = cmaq_wd10m_this[:,j,i]
			cmaq_swdown_pm25_use_stn[nn,cc,:] = cmaq_swdown_this[:,j,i]
			cmaq_pbl_pm25_use_stn[nn,cc,:] = cmaq_pbl_this[:,j,i]
			cmaq_nox_pm25_use_stn[nn,cc,:] = cmaq_nox_this[:,j,i]
			cmaq_noy_pm25_use_stn[nn,cc,:] = cmaq_noy_this[:,j,i]

		for nn in range(n_obs_pm25_val):
			## Find closest j,i index
			dist = abs(cmaq_lat-obs_pm25_lat_val[nn])+abs(cmaq_lon-obs_pm25_lon_val[nn])
			j, i = np.unravel_index(dist.argmin(), dist.shape)
			cmaq_pm25_con_val_stn[nn,cc,:] = cmaq_pm25_this[:,j,i]
			cmaq_t2m_pm25_val_stn[nn,cc,:] = cmaq_t2m_this[:,j,i]
			cmaq_ws10m_pm25_val_stn[nn,cc,:] = cmaq_ws10m_this[:,j,i]
			cmaq_wd10m_pm25_val_stn[nn,cc,:] = cmaq_wd10m_this[:,j,i]
			cmaq_swdown_pm25_val_stn[nn,cc,:] = cmaq_swdown_this[:,j,i]
			cmaq_pbl_pm25_val_stn[nn,cc,:] = cmaq_pbl_this[:,j,i]
			cmaq_nox_pm25_val_stn[nn,cc,:] = cmaq_nox_this[:,j,i]
			cmaq_noy_pm25_val_stn[nn,cc,:] = cmaq_noy_this[:,j,i]
		'''
		
		## Create xarray data arrays using the CMAQ and obs data
		cmaq_raw_str = 'CMAQ raw gridded forecast'
		cmaq_bm3_str = 'CMAQ BM3 bias-corrected gridded "truth"'
		cams_bm3_str = 'CAMS BM3 bias-corrected gridded "truth"'
		use_str = 'stations used for merging and training'
		val_str = 'stations withheld for validation'

		obs_o3_con_use = xr.DataArray(obs_o3_con_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':'AirNow O3 concentration observations at '+use_str, 'units':'ppbv'})
		obs_o3_con_val = xr.DataArray(obs_o3_con_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':'AirNow O3 concentration observations at '+val_str, 'units':'ppbv'})
		obs_pm25_con_use = xr.DataArray(obs_pm25_con_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':'AirNow PM2.5 concentration observations at '+use_str, 'units':'ug m-3'})
		obs_pm25_con_val = xr.DataArray(obs_pm25_con_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':'AirNow PM2.5 concentration observations at '+val_str, 'units':'ug m-3'})

		cams_bm3_o3_con_use = xr.DataArray(cams_bm3_o3_con_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cams_bm3_str+' O3 concentration interpolated to AirNow O3 '+use_str, 'units':'ppbv'})
		cams_bm3_o3_con_val = xr.DataArray(cams_bm3_o3_con_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cams_bm3_str+' O3 concentration interpolated to AirNow O3 '+val_str, 'units':'ppbv'})

		cams_bm3_pm25_con_use = xr.DataArray(cams_bm3_pm25_con_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cams_bm3_str+' PM2.5 concentration interpolated to AirNow PM2.5 '+use_str, 'units':'ug m-3'})
		cams_bm3_pm25_con_val = xr.DataArray(cams_bm3_pm25_con_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cams_bm3_str+' PM2.5 concentration interpolated to AirNow PM2.5 '+val_str, 'units':'ug m-3'})

		cmaq_bm3_o3_con_use = xr.DataArray(cmaq_bm3_o3_con_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cmaq_bm3_str+' O3 concentration interpolated to AirNow O3 '+use_str, 'units':'ppbv'})
		cmaq_bm3_o3_con_val = xr.DataArray(cmaq_bm3_o3_con_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cmaq_bm3_str+' O3 concentration interpolated to AirNow O3 '+val_str, 'units':'ppbv'})

		cmaq_bm3_pm25_con_use = xr.DataArray(cmaq_bm3_pm25_con_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cmaq_bm3_str+' PM2.5 concentration interpolated to AirNow PM2.5 '+use_str, 'units':'ug m-3'})
		cmaq_bm3_pm25_con_val = xr.DataArray(cmaq_bm3_pm25_con_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cmaq_bm3_str+' PM2.5 concentration interpolated to AirNow PM2.5 '+val_str, 'units':'ug m-3'})

		cmaq_raw_o3_con_use = xr.DataArray(cmaq_raw_o3_con_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' O3 concentration interpolated to AirNow O3 '+use_str, 'units':'ppbv'})
		cmaq_raw_o3_con_val = xr.DataArray(cmaq_raw_o3_con_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' O3 concentration interpolated to AirNow O3 '+val_str, 'units':'ppbv'})

		cmaq_raw_pm25_con_use = xr.DataArray(cmaq_raw_pm25_con_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' PM2.5 concentration interpolated to AirNow PM2.5 '+use_str, 'units':'ug m-3'})
		cmaq_raw_pm25_con_val = xr.DataArray(cmaq_raw_pm25_con_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' PM2.5 concentration interpolated to AirNow PM2.5 '+val_str, 'units':'ug m-3'})

		cmaq_raw_t2m_o3_use = xr.DataArray(cmaq_raw_t2m_o3_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' 2-m temperature interpolated to AirNow O3 '+use_str, 'units':'K'})
		cmaq_raw_t2m_o3_val = xr.DataArray(cmaq_raw_t2m_o3_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' 2-m temperature interpolated to AirNow O3 '+val_str, 'units':'K'})
		cmaq_raw_t2m_pm25_use = xr.DataArray(cmaq_raw_t2m_pm25_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' 2-m temperature interpolated to AirNow PM2.5 '+use_str, 'units':'K'})
		cmaq_raw_t2m_pm25_val = xr.DataArray(cmaq_raw_t2m_pm25_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' 2-m temperature interpolated to AirNow PM2.5 '+val_str, 'units':'K'})


		cmaq_raw_ws10m_o3_use = xr.DataArray(cmaq_raw_ws10m_o3_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' 10-m wind speed interpolated to AirNow O3 '+use_str, 'units':'m s-1'})
		cmaq_raw_ws10m_o3_val = xr.DataArray(cmaq_raw_ws10m_o3_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' 10-m wind speed interpolated to AirNow O3 '+val_str, 'units':'m s-1'})
		cmaq_raw_ws10m_pm25_use = xr.DataArray(cmaq_raw_ws10m_pm25_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' 10-m wind speed interpolated to AirNow PM2.5 '+use_str, 'units':'m s-1'})
		cmaq_raw_ws10m_pm25_val = xr.DataArray(cmaq_raw_ws10m_pm25_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' 10-m wind speed interpolated to AirNow PM2.5 '+val_str, 'units':'m s-1'})

		cmaq_raw_wd10m_o3_use = xr.DataArray(cmaq_raw_wd10m_o3_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' 10-m wind direction interpolated to AirNow O3 '+use_str, 'units':'deg'})
		cmaq_raw_wd10m_o3_val = xr.DataArray(cmaq_raw_wd10m_o3_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' 10-m wind direction interpolated to AirNow O3 '+val_str, 'units':'deg'})
		cmaq_raw_wd10m_pm25_use = xr.DataArray(cmaq_raw_wd10m_pm25_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' 10-m wind direction interpolated to AirNow PM2.5 '+use_str, 'units':'deg'})
		cmaq_raw_wd10m_pm25_val = xr.DataArray(cmaq_raw_wd10m_pm25_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' 10-m wind direction interpolated to AirNow PM2.5 '+val_str, 'units':'deg'})

		cmaq_raw_swdown_o3_use = xr.DataArray(cmaq_raw_swdown_o3_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' GHI interpolated to AirNow O3 '+use_str, 'units':'W m-2'})
		cmaq_raw_swdown_o3_val = xr.DataArray(cmaq_raw_swdown_o3_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' GHI interpolated to AirNow O3 '+val_str, 'units':'W m-2'})
		cmaq_raw_swdown_pm25_use = xr.DataArray(cmaq_raw_swdown_pm25_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' GHI interpolated to AirNow PM2.5 '+use_str, 'units':'W m-2'})
		cmaq_raw_swdown_pm25_val = xr.DataArray(cmaq_raw_swdown_pm25_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' GHI interpolated to AirNow PM2.5 '+val_str, 'units':'W m-2'})

		cmaq_raw_pbl_o3_use = xr.DataArray(cmaq_raw_pbl_o3_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' PBL depth interpolated to AirNow O3 '+use_str, 'units':'m'})
		cmaq_raw_pbl_o3_val = xr.DataArray(cmaq_raw_pbl_o3_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' PBL depth interpolated to AirNow O3 '+val_str, 'units':'m'})
		cmaq_raw_pbl_pm25_use = xr.DataArray(cmaq_raw_pbl_pm25_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' PBL depth interpolated to AirNow PM2.5 '+use_str, 'units':'m'})
		cmaq_raw_pbl_pm25_val = xr.DataArray(cmaq_raw_pbl_pm25_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' PBL depth interpolated to AirNow PM2.5 '+val_str, 'units':'m'})

		cmaq_raw_nox_o3_use = xr.DataArray(cmaq_raw_nox_o3_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' NOx interpolated to AirNow O3 '+use_str, 'units':'ppmv'})
		cmaq_raw_nox_o3_val = xr.DataArray(cmaq_raw_nox_o3_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' NOx interpolated to AirNow O3 '+val_str, 'units':'ppmv'})
		cmaq_raw_nox_pm25_use = xr.DataArray(cmaq_raw_nox_pm25_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' NOx interpolated to AirNow PM2.5 '+use_str, 'units':'ppmv'})
		cmaq_raw_nox_pm25_val = xr.DataArray(cmaq_raw_nox_pm25_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' NOx interpolated to AirNow PM2.5 '+val_str, 'units':'ppmv'})

		cmaq_raw_noy_o3_use = xr.DataArray(cmaq_raw_noy_o3_use_stn[:,cc,:],
			coords={'n_obs_use_o3':n_use_o3,'lead_hours':n_leads}, dims=['n_use_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' NOy interpolated to AirNow O3 '+use_str, 'units':'ppmv'})
		cmaq_raw_noy_o3_val = xr.DataArray(cmaq_raw_noy_o3_val_stn[:,cc,:],
			coords={'n_obs_val_o3':n_val_o3,'lead_hours':n_leads}, dims=['n_val_o3','n_leads'],
			attrs={'description':cmaq_raw_str+' NOy interpolated to AirNow O3 '+val_str, 'units':'ppmv'})
		cmaq_raw_noy_pm25_use = xr.DataArray(cmaq_raw_noy_pm25_use_stn[:,cc,:],
			coords={'n_obs_use_pm25':n_use_pm25,'lead_hours':n_leads}, dims=['n_use_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' NOy interpolated to AirNow PM2.5 '+use_str, 'units':'ppmv'})
		cmaq_raw_noy_pm25_val = xr.DataArray(cmaq_raw_noy_pm25_val_stn[:,cc,:],
			coords={'n_obs_val_pm25':n_val_pm25,'lead_hours':n_leads}, dims=['n_val_pm25','n_leads'],
			attrs={'description':cmaq_raw_str+' NOy interpolated to AirNow PM2.5 '+val_str, 'units':'ppmv'})


		## Create new xarray datasets
		ds_use = xr.Dataset(
				data_vars={
					'cams_bm3_pm25_con_use':cams_bm3_pm25_con_use, 'cams_bm3_o3_con_use':cams_bm3_o3_con_use,
					'cmaq_bm3_pm25_con_use':cmaq_bm3_pm25_con_use, 'cmaq_bm3_o3_con_use':cmaq_bm3_o3_con_use,
					'cmaq_raw_pm25_con_use':cmaq_raw_pm25_con_use, 'cmaq_raw_o3_con_use':cmaq_raw_o3_con_use,
					'cmaq_raw_t2m_pm25_use':cmaq_raw_t2m_pm25_use, 'cmaq_raw_t2m_o3_use':cmaq_raw_t2m_o3_use,
					'cmaq_raw_ws10m_pm25_use':cmaq_raw_ws10m_pm25_use, 'cmaq_raw_ws10m_o3_use':cmaq_raw_ws10m_o3_use,
					'cmaq_raw_wd10m_pm25_use':cmaq_raw_wd10m_pm25_use, 'cmaq_raw_wd10m_o3_use':cmaq_raw_wd10m_o3_use,
					'cmaq_raw_swdown_pm25_use':cmaq_raw_swdown_pm25_use, 'cmaq_raw_swdown_o3_use':cmaq_raw_swdown_o3_use,
					'cmaq_raw_pbl_pm25_use':cmaq_raw_pbl_pm25_use, 'cmaq_raw_pbl_o3_use':cmaq_raw_pbl_o3_use,
					'cmaq_raw_nox_pm25_use':cmaq_raw_nox_pm25_use, 'cmaq_raw_nox_o3_use':cmaq_raw_nox_o3_use,
					'cmaq_raw_noy_pm25_use':cmaq_raw_noy_pm25_use, 'cmaq_raw_noy_o3_use':cmaq_raw_noy_o3_use,
					'obs_pm25_con_use':obs_pm25_con_use, 'obs_o3_con_use':obs_o3_con_use,
					'obs_pm25_sid_use':obs_pm25_sid_use, 'obs_o3_sid_use':obs_o3_sid_use,
					'obs_pm25_lon_use':obs_pm25_lon_use, 'obs_o3_lon_use':obs_o3_lon_use,
					'obs_pm25_lat_use':obs_pm25_lat_use, 'obs_o3_lat_use':obs_o3_lat_use,
					'obs_pm25_cmaq_i_use':obs_pm25_cmaq_i_use, 'obs_o3_cmaq_i_use':obs_o3_cmaq_i_use,
					'obs_pm25_cmaq_j_use':obs_pm25_cmaq_j_use, 'obs_o3_cmaq_j_use':obs_o3_cmaq_j_use,
					},
				coords={'n_obs_use_pm25':n_use_pm25, 'n_obs_use_o3':n_use_o3, 'lead_hours':n_leads},
				attrs={'description':'Observations and interpolated gridded variables from CMAQ raw forecasts, bias-corrected CMAQ "truth", and bias-corrected CAMS "truth" at AirNow PM2.5 & O3 stations retained for use in merging and training'},
				)

		ds_val = xr.Dataset(
				data_vars={
					'cams_bm3_pm25_con_val':cams_bm3_pm25_con_val, 'cams_bm3_o3_con_val':cams_bm3_o3_con_val,
					'cmaq_bm3_pm25_con_val':cmaq_bm3_pm25_con_val, 'cmaq_bm3_o3_con_val':cmaq_bm3_o3_con_val,
					'cmaq_raw_pm25_con_val':cmaq_raw_pm25_con_val, 'cmaq_raw_o3_con_val':cmaq_raw_o3_con_val,
					'cmaq_raw_t2m_pm25_val':cmaq_raw_t2m_pm25_val, 'cmaq_raw_t2m_o3_val':cmaq_raw_t2m_o3_val,
					'cmaq_raw_ws10m_pm25_val':cmaq_raw_ws10m_pm25_val, 'cmaq_raw_ws10m_o3_val':cmaq_raw_ws10m_o3_val,
					'cmaq_raw_wd10m_pm25_val':cmaq_raw_wd10m_pm25_val, 'cmaq_raw_wd10m_o3_val':cmaq_raw_wd10m_o3_val,
					'cmaq_raw_swdown_pm25_val':cmaq_raw_swdown_pm25_val, 'cmaq_raw_swdown_o3_val':cmaq_raw_swdown_o3_val,
					'cmaq_raw_pbl_pm25_val':cmaq_raw_pbl_pm25_val, 'cmaq_raw_pbl_o3_val':cmaq_raw_pbl_o3_val,
					'cmaq_raw_nox_pm25_val':cmaq_raw_nox_pm25_val, 'cmaq_raw_nox_o3_val':cmaq_raw_nox_o3_val,
					'cmaq_raw_noy_pm25_val':cmaq_raw_noy_pm25_val, 'cmaq_raw_noy_o3_val':cmaq_raw_noy_o3_val,
					'obs_pm25_con_val':obs_pm25_con_val, 'obs_o3_con_val':obs_o3_con_val,
					'obs_pm25_sid_val':obs_pm25_sid_val, 'obs_o3_sid_val':obs_o3_sid_val,
					'obs_pm25_lon_val':obs_pm25_lon_val, 'obs_o3_lon_val':obs_o3_lon_val,
					'obs_pm25_lat_val':obs_pm25_lat_val, 'obs_o3_lat_val':obs_o3_lat_val,
					'obs_pm25_cmaq_i_val':obs_pm25_cmaq_i_val, 'obs_o3_cmaq_i_val':obs_o3_cmaq_i_val,
					'obs_pm25_cmaq_j_val':obs_pm25_cmaq_j_val, 'obs_o3_cmaq_j_val':obs_o3_cmaq_j_val,
					},
				coords={'n_obs_val_pm25':n_val_pm25, 'n_obs_val_o3':n_val_o3, 'lead_hours':n_leads},
				attrs={'description':'Observations and interpolated gridded variables from CMAQ raw forecasts, bias-corrected CMAQ "truth", and bias-corrected CAMS "truth" at AirNow PM2.5 & O3 stations withheld for validation'},
				)

		## Set the output paths & filenames
		out_dir = out_dir_parent.joinpath(this_cycle_yr,this_cycle_mo)
		out_dir.mkdir(parents=True, exist_ok=True)
		out_dir.chmod(0o755)
		fname_use = out_dir.joinpath('cmaq_raw_cmaq_bm_cams_bm_airnow_pm2.5_o3_'+this_cycle_outfile+'_use.nc')
		fname_val = out_dir.joinpath('cmaq_raw_cmaq_bm_cams_bm_airnow_pm2.5_o3_'+this_cycle_outfile+'_val.nc')

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
