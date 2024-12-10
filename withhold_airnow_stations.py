'''
withhold_airnow_stations.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 18 Aug 2022

This script reads ASCII files that contain AirNow observations. PM2.5 & O3 are in different files.
In each file, there is one row per station, and columns with metadata + 24 hourly observations.
For each hour/column, the script first identifies the non-missing observations. A user-specified
percentage (say, 30%) of the valid obs are randomly withheld for validation, with the rest being
used by downstream processes to blend into the background gridded analysis. The random seed is set
to be the valid time YYYYMMDDHH integer for repeatability. NetCDF output files are written hourly.
Each output file combines the PM2.5 & O3 obs, with separate files for used, withheld, & missing obs
at that valid hour.
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
	data_dir_in_master = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow')
#	data_dir_in_pm25 = data_dir_in_master.joinpath('csv_pm25')
	data_dir_in_pm25 = data_dir_in_master.joinpath('csv_pm25_noQC')
	data_dir_in_o3   = data_dir_in_master.joinpath('csv_o3')
	data_dir_out = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split')

	fmt_yyyymmdd_hh = '%Y%m%d_%H'
	fmt_yyyymmddhh  = '%Y%m%d%H'
	fmt_yyyymmdd = '%Y%m%d'
	fmt_yyyy = '%Y'
	fmt_mm = '%m'
	fmt_hh = '%H'

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
		this_yyyymmdd = this_dt.strftime(fmt_yyyymmdd)
		this_yyyy = this_dt.strftime(fmt_yyyy)
		this_mm = this_dt.strftime(fmt_mm)
		this_hh = this_dt.strftime(fmt_hh)

		file_in_pm25 = data_dir_in_pm25.joinpath('airnow_'+this_yyyymmdd+'_pm25.csv')
		file_in_o3   = data_dir_in_o3.joinpath(  'airnow_'+this_yyyymmdd+'_o3.csv')

		## Read in the file to a pandas dataframe
		df_in_pm25 = pd.read_csv(file_in_pm25, header=None)
		df_in_o3   = pd.read_csv(file_in_o3, header=None)

		## Col 0: Station ID
		## Col 1: CMAQ i-grid point
		## Col 2: CMAQ j-grid point
		## Col 3: Station longitude (deg E)
		## Col 4: Station latitude (deg N)
		## Cols 5-28: Hourly data values
		sid_all_pm25    = df_in_pm25[0].to_numpy()
		cmaq_i_all_pm25 = df_in_pm25[1].to_numpy()
		cmaq_j_all_pm25 = df_in_pm25[2].to_numpy()
		lon_all_pm25    = df_in_pm25[3].to_numpy()
		lat_all_pm25    = df_in_pm25[4].to_numpy()
		sid_all_o3    = df_in_o3[0].to_numpy()
		cmaq_i_all_o3 = df_in_o3[1].to_numpy()
		cmaq_j_all_o3 = df_in_o3[2].to_numpy()
		lon_all_o3    = df_in_o3[3].to_numpy()
		lat_all_o3    = df_in_o3[4].to_numpy()

		## Extract data values for this hour
		pm25_all = df_in_pm25[int(this_hh)+5].to_numpy()
		o3_all   = df_in_o3[  int(this_hh)+5].to_numpy()
		n_pm25_all = len(pm25_all)
		n_o3_all   = len(o3_all)

		## Get indices of missing and valid values (note that Ju-Hye used -999 for pm2.5 and -999000 for o3)
#		inds_mis_pm25 = np.where(vals_all_pm25 == -999.0)[0]
#		inds_mis_o3   = np.where(vals_all_o3   == -999000.0)[0]
#		n_mis_pm25 = len(inds_mis_pm25)
#		n_mis_o3   = len(inds_mis_o3)
		inds_valid_pm25_all = np.where(pm25_all != -999.0)[0]
		inds_valid_o3_all   = np.where(o3_all   != -999000.0)[0]
		n_valid_pm25_all = len(inds_valid_pm25_all)
		n_valid_o3_all   = len(inds_valid_o3_all)

		## How many of the valid obs should be withheld for validation?
		n_val_pm25 = int(n_valid_pm25_all * val_fraction)
		n_val_o3   = int(n_valid_o3_all   * val_fraction)
		n_use_pm25 = n_valid_pm25_all - n_val_pm25
		n_use_o3   = n_valid_o3_all   - n_val_o3

		## There may be times when there are no valid observations, or too few for any to be withheld.
		## Throw a warning in these cases. Code blocks below need to check for these conditions to avoid errors.
		if n_valid_pm25_all == 0:
			print('WARNING: There are no valid PM2.5 observations in this file at time '+this_yyyymmdd_hh+'00:')
			print('   '+str(file_in_pm25))
		elif n_val_pm25 == 0:
			print('WARNING: No PM2.5 observations will be withheld for validation at time '+this_yyyymmdd_hh+'00.')
			print('         The file has only '+str(n_use_pm25)+' valid PM2.5 observation(s).')

		if n_valid_o3_all == 0:
			print('WARNING: There are no valid O3 observations in this file at time '+this_yyyymmdd_hh+'00:')
			print('   '+str(file_in_o3))
		elif n_val_o3 == 0:
			print('WARNING: No O3 observations will be withheld for validation at time '+this_yyyymmdd_hh+'00.')
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
		## For reproducibility, use a random seed that is the integer value of the valid date.
		## Randomly select without replacement the indices to be withheld from the dummy array.
		## Use all other indices from this dummy array.
		if n_val_pm25 > 0:
			np.random.seed(int(this_yyyymmddhh)) ;
			inds_val_pm25 = np.sort(np.random.choice(range(n_valid_pm25_all), n_val_pm25, replace=False))

		if n_val_o3 > 0:
			np.random.seed(int(this_yyyymmddhh)) ;
			inds_val_o3   = np.sort(np.random.choice(range(n_valid_o3_all), n_val_o3, replace=False))

		if n_use_pm25 > 0 and n_val_pm25 > 0:
			inds_use_pm25 = np.delete(range(n_valid_pm25_all), inds_val_pm25)
		elif n_use_pm25 > 0 and n_val_pm25 == 0:
			inds_use_pm25 = np.arange(n_valid_pm25_all)

		if n_use_o3 > 0 and n_val_o3 > 0:
			inds_use_o3 = np.delete(range(n_valid_o3_all), inds_val_o3)
		elif n_use_o3 > 0 and n_val_o3 == 0:
			inds_use_o3 = np.arange(n_valid_o3_all)

#		print(inds_val_pm25)
#		print(inds_val_o3)

		## Now match these dummy array indices to the original input data arrays
		## Then create the arrays from that indexed data
		if n_val_pm25 > 0:
			inds_valid_pm25_val = inds_valid_pm25_all[inds_val_pm25]
			pm25_val_con_data = pm25_all[inds_valid_pm25_val]
			pm25_val_sid_data = sid_all_pm25[inds_valid_pm25_val]
			pm25_val_lon_data = lon_all_pm25[inds_valid_pm25_val]
			pm25_val_lat_data = lat_all_pm25[inds_valid_pm25_val]
			pm25_val_cmaq_i_data = cmaq_i_all_pm25[inds_valid_pm25_val]
			pm25_val_cmaq_j_data = cmaq_j_all_pm25[inds_valid_pm25_val]
		else:
			n_val_pm25 = 1		# need to set to 1 so that the variable can be written to NetCDF
			pm25_val_con_data = np.nan
			pm25_val_sid_data = -999
			pm25_val_lon_data = np.nan
			pm25_val_lat_data = np.nan
			pm25_val_cmaq_i_data = -999
			pm25_val_cmaq_j_data = -999

		pm25_val_con = xr.DataArray(pm25_val_con_data,
								coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
								attrs={'description':'AirNow PM2.5 station concentration', 'units':'ug/m3'})
		pm25_val_sid = xr.DataArray(pm25_val_sid_data,
								coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
								attrs={'description':'AirNow PM2.5 station ID', '_FillValue':-999})
		pm25_val_lon = xr.DataArray(pm25_val_lon_data,
								coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
								attrs={'description':'AirNow PM2.5 station longitude', 'units':'deg E'})
		pm25_val_lat = xr.DataArray(pm25_val_lat_data,
								coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
								attrs={'description':'AirNow PM2.5 station latitude', 'units':'deg N'})
		pm25_val_cmaq_i = xr.DataArray(pm25_val_cmaq_i_data,
									coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
									attrs={'description':'AirNow PM2.5 station nearest CMAQ grid i-index','_FillValue':-999})
		pm25_val_cmaq_j = xr.DataArray(pm25_val_cmaq_j_data,
									coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
									attrs={'description':'AirNow PM2.5 station nearest CMAQ grid j-index','_FillValue':-999})

		if n_use_pm25 > 0:
			inds_valid_pm25_use = inds_valid_pm25_all[inds_use_pm25]
			pm25_use_con_data = pm25_all[inds_valid_pm25_use]
			pm25_use_sid_data = sid_all_pm25[inds_valid_pm25_use]
			pm25_use_lon_data = lon_all_pm25[inds_valid_pm25_use]
			pm25_use_lat_data = lat_all_pm25[inds_valid_pm25_use]
			pm25_use_cmaq_i_data = cmaq_i_all_pm25[inds_valid_pm25_use]
			pm25_use_cmaq_j_data = cmaq_j_all_pm25[inds_valid_pm25_use]
		else:
			n_use_pm25 = 1		# need to set to 1 so that the variable can be written to NetCDF
			pm25_use_con_data = np.nan
			pm25_use_sid_data = -999
			pm25_use_lon_data = np.nan
			pm25_use_lat_data = np.nan
			pm25_use_cmaq_i_data = -999
			pm25_use_cmaq_j_data = -999

		pm25_use_con = xr.DataArray(pm25_use_con_data,
								coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
								attrs={'description':'AirNow PM2.5 station concentration', 'units':'ug/m3'})
		pm25_use_sid = xr.DataArray(pm25_use_sid_data,
								coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
								attrs={'description':'AirNow PM2.5 station ID', '_FillValue':-999})
		pm25_use_lon = xr.DataArray(pm25_use_lon_data,
								coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
								attrs={'description':'AirNow PM2.5 station longitude', 'units':'deg E'})
		pm25_use_lat = xr.DataArray(pm25_use_lat_data,
								coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
								attrs={'description':'AirNow PM2.5 station latitude', 'units':'deg N'})
		pm25_use_cmaq_i = xr.DataArray(pm25_use_cmaq_i_data,
									coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
									attrs={'description':'AirNow PM2.5 station nearest CMAQ grid i-index','_FillValue':-999})
		pm25_use_cmaq_j = xr.DataArray(pm25_use_cmaq_j_data,
									coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
									attrs={'description':'AirNow PM2.5 station nearest CMAQ grid j-index','_FillValue':-999})

		if n_val_o3 > 0:
			inds_valid_o3_val = inds_valid_o3_all[inds_val_o3]
			o3_val_con_data = o3_all[inds_valid_o3_val]
			o3_val_sid_data = sid_all_o3[inds_valid_o3_val]
			o3_val_lon_data = lon_all_o3[inds_valid_o3_val]
			o3_val_lat_data = lat_all_o3[inds_valid_o3_val]
			o3_val_cmaq_i_data = cmaq_i_all_o3[inds_valid_o3_val]
			o3_val_cmaq_j_data = cmaq_j_all_o3[inds_valid_o3_val]
		else:
			n_val_o3 = 1	# need to set to 1 so that the variable can be written to NetCDF
			o3_val_con_data = np.nan
			o3_val_sid_data = -999
			o3_val_lon_data = np.nan
			o3_val_lat_data = np.nan
			o3_val_cmaq_i_data = -999
			o3_val_cmaq_j_data = -999

		o3_val_con = xr.DataArray(o3_val_con_data,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station concentration', 'units':'ppbv'})
		o3_val_sid = xr.DataArray(o3_val_sid_data,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station ID', '_FillValue':-999})
		o3_val_lon = xr.DataArray(o3_val_lon_data,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station longitude', 'units':'deg E'})
		o3_val_lat = xr.DataArray(o3_val_lat_data,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station latitude', 'units':'deg N'})
		o3_val_cmaq_i = xr.DataArray(o3_val_cmaq_i_data,
								coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
								attrs={'description':'AirNow O3 station nearest CMAQ grid i-index', '_FillValue':-999})
		o3_val_cmaq_j = xr.DataArray(o3_val_cmaq_j_data,
								coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
								attrs={'description':'AirNow O3 station nearest CMAQ grid j-index', '_FillValue':-999})

		if n_use_o3 > 0:
			inds_valid_o3_use = inds_valid_o3_all[inds_use_o3]
			o3_use_con_data = o3_all[inds_valid_o3_use]
			o3_use_sid_data = sid_all_o3[inds_valid_o3_use]
			o3_use_lon_data = lon_all_o3[inds_valid_o3_use]
			o3_use_lat_data = lat_all_o3[inds_valid_o3_use]
			o3_use_cmaq_i_data = cmaq_i_all_o3[inds_valid_o3_use]
			o3_use_cmaq_j_data = cmaq_j_all_o3[inds_valid_o3_use]
		else:
			n_use_o3 = 1	# need to set to 1 so that the variable can be written to NetCDF
			o3_use_con_data = np.nan
			o3_use_sid_data = -999
			o3_use_lon_data = np.nan
			o3_use_lat_data = np.nan
			o3_use_cmaq_i_data = -999
			o3_use_cmaq_j_data = -999

		o3_use_con = xr.DataArray(o3_use_con_data,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station concentration', 'units':'ppbv'})
		o3_use_sid = xr.DataArray(o3_use_sid_data,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station ID', '_FillValue':-999})
		o3_use_lon = xr.DataArray(o3_use_lon_data,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station longitude', 'units':'deg E'})
		o3_use_lat = xr.DataArray(o3_use_lat_data,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station latitude', 'units':'deg N'})
		o3_use_cmaq_i = xr.DataArray(o3_use_cmaq_i_data,
								coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
								attrs={'description':'AirNow O3 station nearest CMAQ grid i-index', '_FillValue':-999})
		o3_use_cmaq_j = xr.DataArray(o3_use_cmaq_j_data,
								coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
								attrs={'description':'AirNow O3 station nearest CMAQ grid j-index', '_FillValue':-999})

#		print(pm25_val_con)
#		print(pm25_use_con)
#		print(o3_val_con)
#		print(o3_use_con)

#		print(o3_val_sid)
#		print(o3_val_lon)
#		print(o3_val_lat)
#		print(o3_val_cmaq_i)
#		print(o3_val_cmaq_j)
#		print(o3_val_con.size)
#		print(o3_val_sid.size)
#		print(o3_val_lon.size)
#		print(o3_val_lat.size)
#		print(o3_val_cmaq_i.size)
#		print(o3_val_cmaq_j.size)

		## Create new xarray datasets for both use and val
		ds_val = xr.Dataset(
						data_vars={ 'pm25_con_val':pm25_val_con, 'o3_con_val':o3_val_con,
										'pm25_sid_val':pm25_val_sid, 'o3_sid_val':o3_val_sid,
										'pm25_lon_val':pm25_val_lon, 'o3_lon_val':o3_val_lon,
										'pm25_lat_val':pm25_val_lat, 'o3_lat_val':o3_val_lat,
										'pm25_cmaq_i_val':pm25_val_cmaq_i, 'o3_cmaq_i_val':o3_val_cmaq_i,
										'pm25_cmaq_j_val':pm25_val_cmaq_j, 'o3_cmaq_j_val':o3_val_cmaq_j,
							},
						coords={'n_obs_val_pm25':n_val_pm25, 'n_obs_val_o3':n_val_o3},
						attrs={'description':'Subset of PM2.5 & O3 AirNow observations withheld for validation'},
						)

		ds_use = xr.Dataset(
						data_vars={ 'pm25_con_use':pm25_use_con, 'o3_con_use':o3_use_con,
										'pm25_sid_use':pm25_use_sid, 'o3_sid_use':o3_use_sid,
										'pm25_lon_use':pm25_use_lon, 'o3_lon_use':o3_use_lon,
										'pm25_lat_use':pm25_use_lat, 'o3_lat_use':o3_use_lat,
										'pm25_cmaq_i_use':pm25_use_cmaq_i, 'o3_cmaq_i_use':o3_use_cmaq_i,
										'pm25_cmaq_j_use':pm25_use_cmaq_j, 'o3_cmaq_j_use':o3_use_cmaq_j,
							},
						coords={'n_obs_use_pm25':n_use_pm25, 'n_obs_use_o3':n_use_o3},
						attrs={'description':'Subset of PM2.5 & O3 AirNow observations retained for use in merging'},
						)

		## Set the output paths & filenames
		out_dir = data_dir_out.joinpath('sites_vary',this_yyyy,this_mm)
		out_dir.mkdir(parents=True, exist_ok=True)
		out_dir.chmod(0o777)
		fname_val = out_dir.joinpath('airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_val.nc')
		fname_use = out_dir.joinpath('airnow_pm2.5_o3_'+this_yyyymmdd_hh+'00_use.nc')

		## Write the datasets to NetCDF
		print('Writing file '+str(fname_val))
		ds_val.to_netcdf(fname_val)
		fname_val.chmod(0o644)
		print('Writing file '+str(fname_use))
		ds_use.to_netcdf(fname_use)
		fname_use.chmod(0o644)
		
		

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
