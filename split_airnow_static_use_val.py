'''
split_airnow_static_use_val.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 6 Oct 2022

This script reads in the static-site use/val netCDF files created by tally_airnow_valid_obs.py
for a given multi-month time period (the sites chosen for use/val will differ over different time
periods), and splits them up into individual files for each valid time.
'''

import sys
import argparse
import pathlib
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt

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
	obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static')

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

	## Open single files that contain all times
	all_file = obs_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_all.nc')
	use_file = obs_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_use.nc')
	val_file = obs_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_val.nc')

	print('Reading '+str(all_file))
	ds_all = xr.open_dataset(all_file)
	pm25_con_all = ds_all.pm25_con_all
	pm25_sid_all = ds_all.pm25_sid_all
	pm25_lon_all = ds_all.pm25_lon_all
	pm25_lat_all = ds_all.pm25_lat_all
	pm25_cmaq_i_all = ds_all.pm25_cmaq_i_all
	pm25_cmaq_j_all = ds_all.pm25_cmaq_j_all
	o3_con_all = ds_all.o3_con_all
	o3_sid_all = ds_all.o3_sid_all
	o3_lon_all = ds_all.o3_lon_all
	o3_lat_all = ds_all.o3_lat_all
	o3_cmaq_i_all = ds_all.o3_cmaq_i_all
	o3_cmaq_j_all = ds_all.o3_cmaq_j_all
	n_all_pm25 = len(pm25_sid_all)
	n_all_o3   = len(o3_sid_all)

	print('Reading '+str(use_file))
	ds_use = xr.open_dataset(use_file)
	pm25_con_use = ds_use.pm25_con_use
	pm25_sid_use = ds_use.pm25_sid_use
	pm25_lon_use = ds_use.pm25_lon_use
	pm25_lat_use = ds_use.pm25_lat_use
	pm25_cmaq_i_use = ds_use.pm25_cmaq_i_use
	pm25_cmaq_j_use = ds_use.pm25_cmaq_j_use
	o3_con_use = ds_use.o3_con_use
	o3_sid_use = ds_use.o3_sid_use
	o3_lon_use = ds_use.o3_lon_use
	o3_lat_use = ds_use.o3_lat_use
	o3_cmaq_i_use = ds_use.o3_cmaq_i_use
	o3_cmaq_j_use = ds_use.o3_cmaq_j_use
	n_use_pm25 = len(pm25_sid_use)
	n_use_o3   = len(o3_sid_use)

	print('Reading '+str(val_file))
	ds_val = xr.open_dataset(val_file)
	pm25_con_val = ds_val.pm25_con_val
	pm25_sid_val = ds_val.pm25_sid_val
	pm25_lon_val = ds_val.pm25_lon_val
	pm25_lat_val = ds_val.pm25_lat_val
	pm25_cmaq_i_val = ds_val.pm25_cmaq_i_val
	pm25_cmaq_j_val = ds_val.pm25_cmaq_j_val
	o3_con_val = ds_val.o3_con_val
	o3_sid_val = ds_val.o3_sid_val
	o3_lon_val = ds_val.o3_lon_val
	o3_lat_val = ds_val.o3_lat_val
	o3_cmaq_i_val = ds_val.o3_cmaq_i_val
	o3_cmaq_j_val = ds_val.o3_cmaq_j_val
	n_val_pm25 = len(pm25_sid_val)
	n_val_o3   = len(o3_sid_val)

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

		## Build the DataArrays to be written out
		pm25_all_sid = pm25_sid_all
		pm25_all_lon = pm25_lon_all
		pm25_all_lat = pm25_lat_all
		pm25_all_con = xr.DataArray(pm25_con_all[tt,:].values,
								coords={'n_obs_all_pm25':n_all_pm25}, dims=['n_all_pm25'],
								attrs={'description':'AirNow PM2.5 station concentration', 'units':'ug/m3'})
		pm25_all_cmaq_i = pm25_cmaq_i_all
		pm25_all_cmaq_j = pm25_cmaq_j_all

		o3_all_sid = o3_sid_all
		o3_all_lon = o3_lon_all
		o3_all_lat = o3_lat_all
		o3_all_con = xr.DataArray(o3_con_all[tt,:].values,
							coords={'n_obs_all_o3':n_all_o3}, dims=['n_all_o3'],
							attrs={'description':'AirNow O3 station concentration', 'units':'ppbv'})
		o3_all_cmaq_i = o3_cmaq_i_all
		o3_all_cmaq_j = o3_cmaq_j_all

		pm25_use_sid = pm25_sid_use
		pm25_use_lon = pm25_lon_use
		pm25_use_lat = pm25_lat_use
		pm25_use_con = xr.DataArray(pm25_con_use[tt,:].values,
								coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
								attrs={'description':'AirNow PM2.5 station concentration', 'units':'ug/m3'})
		pm25_use_cmaq_i = pm25_cmaq_i_use
		pm25_use_cmaq_j = pm25_cmaq_j_use

		o3_use_sid = o3_sid_use
		o3_use_lon = o3_lon_use
		o3_use_lat = o3_lat_use
		o3_use_con = xr.DataArray(o3_con_use[tt,:].values,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station concentration', 'units':'ppbv'})
		o3_use_cmaq_i = o3_cmaq_i_use
		o3_use_cmaq_j = o3_cmaq_j_use

		pm25_val_sid = pm25_sid_val
		pm25_val_lon = pm25_lon_val
		pm25_val_lat = pm25_lat_val
		pm25_val_con = xr.DataArray(pm25_con_val[tt,:].values,
								coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
								attrs={'description':'AirNow PM2.5 station concentration', 'units':'ug/m3'})
		pm25_val_cmaq_i = pm25_cmaq_i_val
		pm25_val_cmaq_j = pm25_cmaq_j_val

		o3_val_sid = o3_sid_val
		o3_val_lon = o3_lon_val
		o3_val_lat = o3_lat_val
		o3_val_con = xr.DataArray(o3_con_val[tt,:].values,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station concentration', 'units':'ppbv'})
		o3_val_cmaq_i = o3_cmaq_i_val
		o3_val_cmaq_j = o3_cmaq_j_val

		all_ds = xr.Dataset(
						data_vars={ 'pm25_sid_all':pm25_all_sid, 'o3_sid_all':o3_all_sid,
										'pm25_lon_all':pm25_all_lon, 'o3_lon_all':o3_all_lon,
										'pm25_lat_all':pm25_all_lat, 'o3_lat_all':o3_all_lat,
										'pm25_con_all':pm25_all_con, 'o3_con_all':o3_all_con,
										'pm25_cmaq_i_all':pm25_all_cmaq_i, 'o3_cmaq_i_all':o3_all_cmaq_i,
										'pm25_cmaq_j_all':pm25_all_cmaq_j, 'o3_cmaq_j_all':o3_all_cmaq_j,
							},
						coords={'n_obs_all_pm25':n_all_pm25, 'n_obs_all_o3':n_all_o3},
						attrs={'description':'All PM2.5 & O3 AirNow observations for all sites with any valid reports during '+date_beg+' to '+date_end},
						)

		use_ds = xr.Dataset(
						data_vars={ 'pm25_sid_use':pm25_use_sid, 'o3_sid_use':o3_use_sid,
										'pm25_lon_use':pm25_use_lon, 'o3_lon_use':o3_use_lon,
										'pm25_lat_use':pm25_use_lat, 'o3_lat_use':o3_use_lat,
										'pm25_con_use':pm25_use_con, 'o3_con_use':o3_use_con,
										'pm25_cmaq_i_use':pm25_use_cmaq_i, 'o3_cmaq_i_use':o3_use_cmaq_i,
										'pm25_cmaq_j_use':pm25_use_cmaq_j, 'o3_cmaq_j_use':o3_use_cmaq_j,
							},
						coords={'n_obs_use_pm25':n_use_pm25, 'n_obs_use_o3':n_use_o3},
						attrs={'description':'Subset of PM2.5 & O3 AirNow observations retained for use in merging & AnEn training. All stations have >= 80% uptime during '+date_beg+' to '+date_end+'.'},
						)

		val_ds = xr.Dataset(
						data_vars={ 'pm25_sid_val':pm25_val_sid, 'o3_sid_val':o3_val_sid,
										'pm25_lon_val':pm25_val_lon, 'o3_lon_val':o3_val_lon,
										'pm25_lat_val':pm25_val_lat, 'o3_lat_val':o3_val_lat,
										'pm25_con_val':pm25_val_con, 'o3_con_val':o3_val_con,
										'pm25_cmaq_i_val':pm25_val_cmaq_i, 'o3_cmaq_i_val':o3_val_cmaq_i,
										'pm25_cmaq_j_val':pm25_val_cmaq_j, 'o3_cmaq_j_val':o3_val_cmaq_j,
							},
						coords={'n_obs_val_pm25':n_val_pm25, 'n_obs_val_o3':n_val_o3},
						attrs={'description':'Subset of PM2.5 & O3 AirNow observations withheld for validation. These are all the stations that have any valid reports during '+date_beg+' to '+date_end+' that were not used for AnEn training.'},
						)

		## Now write out the files for this valid time
		out_dir = obs_dir.joinpath(date_beg+'-'+date_end, this_yr, this_mo)
		out_dir.mkdir(parents=True, exist_ok=True)

		fname_all = out_dir.joinpath('airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_all.nc')
		fname_use = out_dir.joinpath('airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_use.nc')
		fname_val = out_dir.joinpath('airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_val.nc')
		print('Writing '+str(fname_all))
		all_ds.to_netcdf(fname_all)
		print('Writing '+str(fname_use))
		use_ds.to_netcdf(fname_use)
		print('Writing '+str(fname_val))
		val_ds.to_netcdf(fname_val)


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
