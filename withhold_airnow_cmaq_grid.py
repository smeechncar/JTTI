'''
withhold_airnow_cmaq_grid.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 23 Jun 2022

This script reads in hourly NetCDF files that have valid AirNow PM2.5 and O3 observations
interpolated to the CMAQ grid. A user-specified percentage (say, 30%) are withheld for validation,
with the rest being use by downstream processes to blend into the background gridded analysis.
This script writes out two NetCDF files, one with the obs to be use, one with obs to be withheld.
'''

import sys
import argparse
import pathlib
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr

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
	## What fraction of valid observations at each time should be withheld for validation?
	val_fraction = 0.30

	## Set directories
	data_dir_in  = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','nc_pm25_o3_cmaq_grid')
	data_dir_out = data_dir_in

	fmt_yyyymmdd_hh = '%Y%m%d_%H'
	fmt_yyyymmddhh  = '%Y%m%d%H'

	## Create datetime array
	dt_beg = dt.datetime.strptime(date_beg, fmt_yyyymmdd_hh)
	dt_end = dt.datetime.strptime(date_end, fmt_yyyymmdd_hh)
	dt_array = pd.date_range(start=dt_beg, end=dt_end, freq='1H')
	n_times = len(dt_array)

	## Loop over times
	for tt in range(n_times):
		this_dt = dt_array[tt]
		this_yyyymmdd_hh = this_dt.strftime(fmt_yyyymmdd_hh)
		this_yyyymmddhh  = this_dt.strftime(fmt_yyyymmddhh)

		## Read the data from the input file
		file_in = data_dir_in.joinpath('AirNOW_regrid_cmaq_'+this_yyyymmddhh+'.nc')
		print('Reading '+str(file_in))
		ds_in = xr.open_dataset(file_in)
		latitude  = ds_in.latitude
		longitude = ds_in.longitude
		pm25_all = ds_in.PM25_O
		o3_all   = ds_in.O3_O
		shape_all = np.shape(pm25_all)

		## Make these 3D PM2.5 and O3 xarray data arrays into 1D numpy arrays for now, so it's easier to work with
		pm25_all_1d = np.ravel(pm25_all)
		o3_all_1d   = np.ravel(o3_all)

		## Which values in the gridded arrays are valid? Get those indices.
		## (Grid points not at AirNow stations are assigned nan.)
#		inds_valid_pm25_all = np.argwhere(~np.isnan(pm25_all.values))
#		inds_valid_o3_all   = np.argwhere(~np.isnan(o3_all.values))
		inds_valid_pm25_all = np.argwhere(~np.isnan(pm25_all_1d))
		inds_valid_o3_all   = np.argwhere(~np.isnan(o3_all_1d))
		n_valid_pm25_all = len(inds_valid_pm25_all)
		n_valid_o3_all   = len(inds_valid_o3_all)

		## How many obs should be withheld for validation?
		n_val_pm25 = int(n_valid_pm25_all * val_fraction)
		n_val_o3   = int(n_valid_o3_all * val_fraction)
		n_use_pm25 = n_valid_pm25_all - n_val_pm25
		n_use_o3   = n_valid_o3_all - n_val_o3

		## There may be times when there are no valid observations, or too few for any to be withheld.
		## Throw a warning in these cases. Code blocks below need to check for these conditions to avoid errors.
		if n_valid_pm25_all == 0:
			print('WARNING: There are no valid PM2.5 observations in this file!')
			print('   '+str(file_in))
		elif n_val_pm25 == 0:
			print('WARNING: No PM2.5 observations will be withheld for validation at this time.')
			print('         The file has only '+str(n_use_pm25)+' valid PM2.5 observation(s).')

		if n_valid_o3_all == 0:
			print('WARNING: There are no valid O3 observations in this file!')
			print('   '+str(file_in))
		elif n_val_o3 == 0:
			print('WARNING: No O3 observations will be withheld for validation at this time.')
			print('         The file has only '+str(n_use_o3)+' valid O3 observation(s).')

		'''
		## Helpful for debugging
		print('')
		print('n_valid_pm25_all = '+str(n_valid_pm25_all))
		print('n_use_pm25 = '+str(n_use_pm25))
		print('n_val_pm25 = '+str(n_val_pm25))
		print('')
		print('n_valid_o3_all = '+str(n_valid_o3_all))
		print('n_use_o3 = '+str(n_use_o3))
		print('n_val_o3 = '+str(n_val_o3))
		print('')
		'''

		## Create a 1D dummy array the length of the number of valid observations.
		## Randomly select without replacement the indices to be withheld from the dummy array.
		## Use all other indices from this dummy array.
		if n_val_pm25 > 0:
			inds_val_pm25 = np.sort(np.random.choice(range(n_valid_pm25_all), n_val_pm25, replace=False))

		if n_val_o3 > 0:
			inds_val_o3   = np.sort(np.random.choice(range(n_valid_o3_all), n_val_o3, replace=False))

		if n_use_pm25 > 0 and n_val_pm25 > 0:
			inds_use_pm25 = np.delete(range(n_valid_pm25_all), inds_val_pm25)
		elif n_use_pm25 > 0 and n_val_pm25 == 0:
			inds_use_pm25 = np.arange(n_valid_pm25_all)

		if n_use_o3 > 0 and n_val_o3 > 0:
			inds_use_o3 = np.delete(range(n_valid_o3_all), inds_val_o3)
		elif n_use_o3 > 0 and n_val_o3 == 0:
			inds_use_o3 = np.arange(n_valid_o3_all)

		## Now match these dummy array indices to the original (flattened) input data arrays
		if n_val_pm25 > 0:
			inds_valid_pm25_val = inds_valid_pm25_all[inds_val_pm25]
		if n_use_pm25 > 0:
			inds_valid_pm25_use = inds_valid_pm25_all[inds_use_pm25]
		if n_val_o3 > 0:
			inds_valid_o3_val = inds_valid_o3_all[inds_val_o3]
		if n_use_o3 > 0:
			inds_valid_o3_use = inds_valid_o3_all[inds_use_o3]

		## Make deep copies of the PM2.5 and O3 input data arrays for both the use and val outputs
		pm25_use = pm25_all.copy()
		pm25_val = pm25_all.copy()
		o3_use = o3_all.copy()
		o3_val = o3_all.copy()

		## Flatten these copies into 1D arrays as before
		pm25_use_1d = np.ravel(pm25_use)
		pm25_val_1d = np.ravel(pm25_val)
		o3_use_1d = np.ravel(o3_use)
		o3_val_1d = np.ravel(o3_val)

		## For the (flattened) output use data array, make the val indices nan, and vice-versa.
		## However, if there are obs to use but none to withhold (or vice-versa), do nothing--keep the array as-is.
		## This will retain only the use data in the use array, and val data in the val array.
		if n_use_pm25 > 0:
			if n_val_pm25 > 0:
				pm25_use_1d[inds_valid_pm25_val] = np.nan
		else:
			pm25_use_1d[:] = np.nan

		if n_val_pm25 > 0:
			if n_use_pm25 > 0:
				pm25_val_1d[inds_valid_pm25_use] = np.nan
		else:
			pm25_val_1d[:] = np.nan

		if n_use_o3 > 0:
			if n_val_o3 > 0:
				o3_use_1d[inds_valid_o3_val] = np.nan
		else:
			o3_val_1d[:] = np.nan

		if n_val_o3 > 0:
			if n_use_o3 > 0:
				o3_val_1d[inds_valid_o3_use] = np.nan
		else:
			o3_val_1d[:] = np.nan

#		inds_pm25_val_1d = np.argwhere(~np.isnan(pm25_val_1d))
#		inds_pm25_use_1d = np.argwhere(~np.isnan(pm25_use_1d))
#		n_pm25_val = len(inds_pm25_val_1d)
#		n_pm25_use = len(inds_pm25_use_1d)

		## Now reshape the flattened use and val arrays back to the shape of the original array
		pm25_use.values = np.reshape(pm25_use_1d, shape_all)
		pm25_val.values = np.reshape(pm25_val_1d, shape_all)
		o3_use.values = np.reshape(o3_use_1d, shape_all)
		o3_val.values = np.reshape(o3_val_1d, shape_all)
#		inds_valid_pm25_all_3d = np.argwhere(~np.isnan(pm25_all.values))
#		inds_valid_pm25_val_3d = np.argwhere(~np.isnan(pm25_val.values))
#		inds_valid_pm25_use_3d = np.argwhere(~np.isnan(pm25_use.values))

		## Create new xarray datasets for both use and val
		## Start by copying the original input dataset, then modify variables & attributes as needed
		ds_use = ds_in.copy()
		ds_val = ds_in.copy()
		ds_use['PM25'] = pm25_use
		ds_use['O3'] = o3_use
		ds_val['PM25'] = pm25_val
		ds_val['O3'] = o3_val
		del(ds_use['PM25_O'])
		del(ds_val['PM25_O'])
		del(ds_use['O3_O'])
		del(ds_val['O3_O'])
		del(ds_use.attrs['AirNOW_data_on_CMAQ_Grid'])
		del(ds_val.attrs['AirNOW_data_on_CMAQ_Grid'])
		ds_use.attrs['description'] = 'Subset of AirNow observations retained for use in merging (interpolated to the CMAQ grid)'
		ds_val.attrs['description'] = 'Subset of AirNow observations withheld for validation (interpolated to the CMAQ grid)'
#		ds_out_use = ds_use[['latitude','longitude','PM25','O3']]
#		ds_out_val = ds_val[['latitude','longitude','PM25','O3']]
#		del(ds_out_use.attrs['AirNOW_data_on_CMAQ_Grid'])
#		del(ds_out_val.attrs['AirNOW_data_on_CMAQ_Grid'])
#		ds_out_use.attrs['description'] = 'Subset of AirNow observations retained for use in merging (interpolated to the CMAQ grid)'
#		ds_out_val.attrs['description'] = 'Subset of AirNow observations withheld for validation (interpolated to the CMAQ grid)'

		## Write these output datasets to new NetCDF files
		file_out_use = data_dir_out.joinpath('use','airnow_regrid_cmaq_'+this_yyyymmdd_hh+'00_use.nc')
		file_out_val = data_dir_out.joinpath('val','airnow_regrid_cmaq_'+this_yyyymmdd_hh+'00_val.nc')
		print('Writing '+str(file_out_use))
		ds_use.to_netcdf(path=file_out_use, format='NETCDF4')
#		ds_out_use.to_netcdf(path=file_out_use, format='NETCDF4')
		print('Writing '+str(file_out_val))
		ds_val.to_netcdf(path=file_out_val, format='NETCDF4')
		

if __name__ == '__main__':
	now_time_beg = dt.datetime.utcnow()
	date_beg, date_end = parse_args()
	main(date_beg, date_end)
	now_time_end = dt.datetime.utcnow()
	run_time_tot = now_time_end - now_time_beg
	now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
	now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
	print('\nScript completed successfully.')
	print('   Beg time: '+now_time_beg_str)
	print('   End time: '+now_time_end_str)
	print('   Run time: '+str(run_time_tot)+'\n')
