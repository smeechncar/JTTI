'''
extract_cmaq_at_airnow_use_val.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 8 Sep 2022

This script first reads in AirNow obs NetCDF files that have already been separated into use & val
sets, for use in merging & validation, respectively. Second, it reads in a CMAQ forecast file to
get the grid coordinates of those stations and spatially interpolate the gridded values to the
AirNow sites via bilinear interpolation. CMAQ PM2.5 & O3 files are produced hourly. The output
files with CMAQ-interpolated values at AirNow sites are in NetCDF.

27 Dec 2022: Updated to use nearest neighbor interpolation instead of RectBivariateSpline.
   RectBivariateSpline gives strange results for curvilinear grids, especially far from the middle.
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
	parser.add_argument('date_beg', help='beginning date/time to be processed [YYYYMMDD_HH]')
	parser.add_argument('-l', '--date_end', default=None, help='ending date/time to be processed [YYYYMMDD_HH]')

	args = parser.parse_args()
	date_beg = args.date_beg
	date_end = args.date_end

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

	return date_beg, date_end


def main(date_beg, date_end):
	sites_vary = False
	sites_static = True

	## Set directories
	if sites_vary:
		obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_vary')
	if sites_static:
		obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static',date_beg+'-'+date_end)
	cmaq_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI', 'bcdata','cmaq_hourly')

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

	get_cmaq_lat_lon = True
	get_cmaq_proj = True

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

		## Second, read in CMAQ files (Ju-Hye has all CMAQ hourly files in one directory)
		cmaq_fname_in_this   = cmaq_dir.joinpath('cmaq_airnow_pm2.5_o3_'+this_yyyymmddhh+'.nc')

		if cmaq_fname_in_this.is_file():
#			print('Reading '+str(cmaq_fname_in_this))
			ds_cmaq = xr.open_dataset(cmaq_fname_in_this)

			if get_cmaq_lat_lon:
				cmaq_lat = ds_cmaq.latitude
				cmaq_lon = ds_cmaq.longitude
				get_cmaq_lat_lon = False

			cmaq_o3_this   = ds_cmaq.O3
			cmaq_pm25_this = ds_cmaq.PM25

		else:
			print('WARNING: File '+str(cmaq_fname_in_this)+' does not exist!')

		'''
		## Now bilinearly interpolate the CMAQ O3 & PM2.5 fields to the AirNow station locations
		## CMAQ is on a Lambert conformal grid, so extra steps are required before calculating the spline

		## First, read attributes from a template CMAQ file to get projection parameters
		if get_cmaq_proj:
			cmaq_coord_file = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
			print('Reading CMAQ coordinate data from '+str(cmaq_coord_file))
			cmaq_coord_ds = xr.open_dataset(cmaq_coord_file)
			cmaq_truelat1 = cmaq_coord_ds.attrs['P_ALP']
			cmaq_truelat2 = cmaq_coord_ds.attrs['P_BET']
			cmaq_stand_lon = cmaq_coord_ds.attrs['P_GAM']
			cmaq_cen_lat = cmaq_coord_ds.attrs['YCENT']
			cmaq_cen_lon = cmaq_coord_ds.attrs['XCENT']
			cmaq_dx = cmaq_coord_ds.attrs['XCELL']
			cmaq_dy = cmaq_coord_ds.attrs['YCELL']
			cmaq_we_dim = cmaq_coord_ds.attrs['NCOLS']
			cmaq_ns_dim = cmaq_coord_ds.attrs['NROWS']
			CMAQ_EARTH_RAD = 6370000
			get_cmaq_proj = False

		## Second, interpolate AirNow locations into CMAQ (x,y) space
		## Get the projection
		p = pyproj.Proj(proj='lcc', lon_0=cmaq_cen_lon, lat_0=cmaq_cen_lat, lat_1=cmaq_truelat1, lat_2=cmaq_truelat2,
				R=CMAQ_EARTH_RAD)
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
		cmaq_o3_ppb_spline = sp.interpolate.RectBivariateSpline(cmaq_j_arr, cmaq_i_arr, cmaq_o3_this, kx=1, ky=1)
		cmaq_o3_ppb_airnow_use = cmaq_o3_ppb_spline.ev(cmaq_j_o3_use, cmaq_i_o3_use)
		cmaq_o3_ppb_airnow_val = cmaq_o3_ppb_spline.ev(cmaq_j_o3_val, cmaq_i_o3_val)

		cmaq_pm25_ug_spline = sp.interpolate.RectBivariateSpline(cmaq_j_arr, cmaq_i_arr, cmaq_pm25_this, kx=1, ky=1)
		cmaq_pm25_ug_airnow_use = cmaq_pm25_ug_spline.ev(cmaq_j_pm25_use, cmaq_i_pm25_use)
		cmaq_pm25_ug_airnow_val = cmaq_pm25_ug_spline.ev(cmaq_j_pm25_val, cmaq_j_pm25_val)
		'''

		## Create xarray data arrays using this data
		n_use_o3 = len(obs_o3_lat_use)
		n_val_o3 = len(obs_o3_lat_val)
		n_use_pm25 = len(obs_pm25_lat_use)
		n_val_pm25 = len(obs_pm25_lat_val)

		## Using RectBivariateSpline for interpolation doesn't work on curvilinear grids like CAMS, WRF, etc.
		## Instead, do nearest-neighbor interpolation to AirNow stations
		print('Interpolating nearest-neighbor CMAQ values to AirNow stations')
		cmaq_o3_ppb_airnow_use = np.full(n_use_o3, np.nan)
		cmaq_o3_ppb_airnow_val = np.full(n_val_o3, np.nan)
		cmaq_pm25_ug_airnow_use = np.full(n_use_pm25, np.nan)
		cmaq_pm25_ug_airnow_val = np.full(n_val_pm25, np.nan)

		## Get and store the nearest-neighbor grid indices for the AirNow use & val stations for PM2.5 & O3
		## Only need to do this once
		if tt == 0:
			print('   Finding grid inds for o3_use')
			j_o3_use, i_o3_use = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_o3_lat_use, obs_o3_lon_use)
			print('   Finding grid inds for o3_val')
			j_o3_val, i_o3_val = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_o3_lat_val, obs_o3_lon_val)
			print('   Finding grid inds for pm25_use')
			j_pm25_use, i_pm25_use = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_pm25_lat_use, obs_pm25_lon_use)
			print('   Finding grid inds for pm25_val')
			j_pm25_val, i_pm25_val = gen_funcs.find_grid_inds(cmaq_lat, cmaq_lon, obs_pm25_lat_val, obs_pm25_lon_val)
			

		## Extract nearest-neighbor values at AirNow stations
		for nn in range(n_use_o3):
#			## Find closest j,i index
#			dist = abs(cmaq_lat-obs_o3_lat_use[nn])+abs(cmaq_lon-obs_o3_lon_use[nn])
#			j, i = np.unravel_index(dist.argmin(), dist.shape)
#			cmaq_o3_ppb_airnow_use[nn] = cmaq_o3_this[j,i]
			cmaq_o3_ppb_airnow_use[nn] = cmaq_o3_this[j_o3_use[nn], i_o3_use[nn]]

		for nn in range(n_val_o3):
#			dist = abs(cmaq_lat-obs_o3_lat_val[nn])+abs(cmaq_lon-obs_o3_lon_val[nn])
#			j, i = np.unravel_index(dist.argmin(), dist.shape)
#			cmaq_o3_ppb_airnow_val[nn] = cmaq_o3_this[j,i]
			cmaq_o3_ppb_airnow_val[nn] = cmaq_o3_this[j_o3_val[nn], i_o3_val[nn]]

		for nn in range(n_use_pm25):
#			dist = abs(cmaq_lat-obs_pm25_lat_use[nn])+abs(cmaq_lon-obs_pm25_lon_use[nn])
#			j, i = np.unravel_index(dist.argmin(), dist.shape)
#			cmaq_pm25_ug_airnow_use[nn] = cmaq_pm25_this[j,i]
			cmaq_pm25_ug_airnow_use[nn] = cmaq_pm25_this[j_pm25_use[nn], i_pm25_use[nn]]

		for nn in range(n_val_pm25):
#			dist = abs(cmaq_lat-obs_pm25_lat_val[nn])+abs(cmaq_lon-obs_pm25_lon_val[nn])
#			j, i = np.unravel_index(dist.argmin(), dist.shape)
#			cmaq_pm25_ug_airnow_val[nn] = cmaq_pm25_this[j,i]
			cmaq_pm25_ug_airnow_val[nn] = cmaq_pm25_this[j_pm25_val[nn], i_pm25_val[nn]]

		cmaq_str = 'CMAQ raw gridded forecast'

		cmaq_o3_con_use = xr.DataArray(cmaq_o3_ppb_airnow_use,
									coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
									attrs={'description':cmaq_str+' O3 concentration interpolated to AirNow O3 stations used for merging', 'units':'ppbv'})
		cmaq_o3_con_val = xr.DataArray(cmaq_o3_ppb_airnow_val,
									coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
									attrs={'description':cmaq_str+' O3 concentration interpolated to AirNow O3 stations withheld for validation', 'units':'ppbv'})

		cmaq_pm25_con_use = xr.DataArray(cmaq_pm25_ug_airnow_use,
									coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
									attrs={'description':cmaq_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations used for merging', 'units':'ug/m3'})
		cmaq_pm25_con_val = xr.DataArray(cmaq_pm25_ug_airnow_val,
									coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
									attrs={'description':cmaq_str+' PM2.5 concentration interpolated to AirNow PM2.5 stations withheld for validation', 'units':'ug/m3'})

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
						coords={'n_obs_use_pm25':n_use_pm25, 'n_obs_use_o3':n_use_o3},
						attrs={'description':'Subset of PM2.5 & O3 AirNow observations and interpolated '+cmaq_str+' values retained for use in merging'},
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
						coords={'n_obs_val_pm25':n_val_pm25, 'n_obs_val_o3':n_val_o3},
						attrs={'description':'Subset of PM2.5 & O3 AirNow observations and interpolated '+cmaq_str+' values withheld for validation'},
						)

		## Set the output paths & filenames
		if sites_vary:
			out_dir = cmaq_dir.joinpath('sites_vary',this_yr,this_mo)
			out_dir.mkdir(parents=True, exist_ok=True)
			out_dir.chmod(0o755)
			fname_use = out_dir.joinpath('cmaq_airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_use.nc')
			fname_val = out_dir.joinpath('cmaq_airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_val.nc')
		if sites_static:
			out_dir = cmaq_dir.joinpath('sites_static',this_yr,this_mo)
			out_dir.mkdir(parents=True, exist_ok=True)
			out_dir.chmod(0o755)
			fname_use = out_dir.joinpath('cmaq_airnow_pm2.5_o3_static_'+this_yyyymmdd_hh+'00_use.nc')
			fname_val = out_dir.joinpath('cmaq_airnow_pm2.5_o3_static_'+this_yyyymmdd_hh+'00_val.nc')

		## Write the datasets to NetCDF
		print('Writing '+str(fname_use))
		ds_use.to_netcdf(fname_use)
		fname_use.chmod(0o644)
		print('Writing '+str(fname_val))
		ds_val.to_netcdf(fname_val)
		fname_val.chmod(0o644)


if __name__ == '__main__':
	now_time_beg = dt.datetime.utcnow()
	date_beg, date_end, = parse_args()
	main(date_beg, date_end)
	now_time_end = dt.datetime.utcnow()
	run_time_tot = now_time_end - now_time_beg
	now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
	now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
	print('\nScript completed successfully.')
	print('   Beg time: '+now_time_beg_str)
	print('   End time: '+now_time_end_str)
	print('   Run time: '+str(run_time_tot)+'\n')
