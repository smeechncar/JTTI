'''
tally_airnow_valid_obs.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 30 Sep 2022

This script reads in the AirNow use & val hourly netCDF files for O3 & PM2.5, and tracks the count
of valid hourly observations for each AirNow station. When a new station is encountered, it is
appended to the station list. Statistics are then presented as maps or histograms. A subset of
stations that have good reporting up-time are then set aside for use in AnEn training or
validation. Lists of the two sets of stations are then written out.
'''

import os
import sys
import argparse
import pathlib
import copy
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import wrf
import cartopy
import cartopy.crs as ccrs
from functions import map_functions

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


def make_hist_plot(fname, var, suptitle, title, xlabel, ylabel, suptitle_y=0.95):
	print('-- Creating plot: '+str(fname))
	counts, edges, bars = plt.hist(var)
	plt.bar_label(bars)
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.suptitle(suptitle, y=suptitle_y)
	plt.title(title)
	plt.grid(True, axis='both')
	plt.savefig(fname)
	plt.close()


def main(date_beg, date_end):
	## Set directories
	obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_vary')
	out_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static')

	plot_type = 'png'
	plot_hist = True
	plot_maps = True

	## Uptime thresholds for selection of stations for use in AnEn training
	up_thresh_pm25 = 80.0
	up_thresh_o3 = 80.0

	## For the stations exceeding the uptime threshold, what percentage should be randomly withheld for training?
	use_fraction = 0.70
	val_fraction = 1.0 - use_fraction

	pm25_str = '$\mathregular{PM_{2.5}}$'
	o3_str = '$\mathregular{O_{3}}$'

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

		## First, read in the use & val AirNow PM2.5 & O3 observations for this time
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

		## On the first time through the loop, create new arrays
		if tt == 0:
			obs_pm25_sid_all = copy.deepcopy(obs_pm25_sid_use)
			obs_pm25_lon_all = copy.deepcopy(obs_pm25_lon_use)
			obs_pm25_lat_all = copy.deepcopy(obs_pm25_lat_use)
			obs_pm25_con_all = copy.deepcopy(obs_pm25_con_use)
			obs_pm25_cmaq_i_all = copy.deepcopy(obs_pm25_cmaq_i_use)
			obs_pm25_cmaq_j_all = copy.deepcopy(obs_pm25_cmaq_j_use)

			obs_o3_sid_all = copy.deepcopy(obs_o3_sid_use)
			obs_o3_lon_all = copy.deepcopy(obs_o3_lon_use)
			obs_o3_lat_all = copy.deepcopy(obs_o3_lat_use)
			obs_o3_con_all = copy.deepcopy(obs_o3_con_use)
			obs_o3_cmaq_i_all = copy.deepcopy(obs_o3_cmaq_i_use)
			obs_o3_cmaq_j_all = copy.deepcopy(obs_o3_cmaq_j_use)

			obs_pm25_sid_all = np.append(obs_pm25_sid_all, obs_pm25_sid_val)
			obs_pm25_lon_all = np.append(obs_pm25_lon_all, obs_pm25_lon_val)
			obs_pm25_lat_all = np.append(obs_pm25_lat_all, obs_pm25_lat_val)
			obs_pm25_con_all = np.append(obs_pm25_con_all, obs_pm25_con_val)
			obs_pm25_con_all = np.expand_dims(obs_pm25_con_all, axis=0)
			obs_pm25_cmaq_i_all = np.append(obs_pm25_cmaq_i_all, obs_pm25_cmaq_i_val)
			obs_pm25_cmaq_j_all = np.append(obs_pm25_cmaq_j_all, obs_pm25_cmaq_j_val)

			obs_o3_sid_all = np.append(obs_o3_sid_all, obs_o3_sid_val)
			obs_o3_lon_all = np.append(obs_o3_lon_all, obs_o3_lon_val)
			obs_o3_lat_all = np.append(obs_o3_lat_all, obs_o3_lat_val)
			obs_o3_con_all = np.append(obs_o3_con_all, obs_o3_con_val)
			obs_o3_con_all = np.expand_dims(obs_o3_con_all, axis=0)
			obs_o3_cmaq_i_all = np.append(obs_o3_cmaq_i_all, obs_o3_cmaq_i_val)
			obs_o3_cmaq_j_all = np.append(obs_o3_cmaq_j_all, obs_o3_cmaq_j_val)

			n_pm25_this = len(obs_pm25_sid_all)
			n_o3_this = len(obs_o3_sid_all)

			obs_pm25_cnt_all = np.full(n_pm25_this, 1)
			obs_o3_cnt_all = np.full(n_o3_this, 1)

		## On all subsequent times through the loop, search element by element and increment or append as needed
		else:
			## First, append a new row to the time axis for concentrations
			n_pm25_prev = len(obs_pm25_sid_all)
			n_o3_prev = len(obs_o3_sid_all)
			obs_pm25_con_all = np.append(obs_pm25_con_all, np.full([1,n_pm25_prev], np.nan), axis=0)
			obs_o3_con_all = np.append(obs_o3_con_all, np.full([1,n_o3_prev], np.nan), axis=0)

			## Second, build temporary arrays as in the first time case
			obs_pm25_sid_this = copy.deepcopy(obs_pm25_sid_use)
			obs_pm25_lon_this = copy.deepcopy(obs_pm25_lon_use)
			obs_pm25_lat_this = copy.deepcopy(obs_pm25_lat_use)
			obs_pm25_con_this = copy.deepcopy(obs_pm25_con_use)
			obs_pm25_cmaq_i_this = copy.deepcopy(obs_pm25_cmaq_i_use)
			obs_pm25_cmaq_j_this = copy.deepcopy(obs_pm25_cmaq_j_use)

			obs_o3_sid_this = copy.deepcopy(obs_o3_sid_use)
			obs_o3_lon_this = copy.deepcopy(obs_o3_lon_use)
			obs_o3_lat_this = copy.deepcopy(obs_o3_lat_use)
			obs_o3_con_this = copy.deepcopy(obs_o3_con_use)
			obs_o3_cmaq_i_this = copy.deepcopy(obs_o3_cmaq_i_use)
			obs_o3_cmaq_j_this = copy.deepcopy(obs_o3_cmaq_j_use)

			obs_pm25_sid_this = np.append(obs_pm25_sid_this, obs_pm25_sid_val)
			obs_pm25_lon_this = np.append(obs_pm25_lon_this, obs_pm25_lon_val)
			obs_pm25_lat_this = np.append(obs_pm25_lat_this, obs_pm25_lat_val)
			obs_pm25_con_this = np.append(obs_pm25_con_this, obs_pm25_con_val)
			obs_pm25_cmaq_i_this = np.append(obs_pm25_cmaq_i_this, obs_pm25_cmaq_i_val)
			obs_pm25_cmaq_j_this = np.append(obs_pm25_cmaq_j_this, obs_pm25_cmaq_j_val)

			obs_o3_sid_this = np.append(obs_o3_sid_this, obs_o3_sid_val)
			obs_o3_lon_this = np.append(obs_o3_lon_this, obs_o3_lon_val)
			obs_o3_lat_this = np.append(obs_o3_lat_this, obs_o3_lat_val)
			obs_o3_con_this = np.append(obs_o3_con_this, obs_o3_con_val)
			obs_o3_cmaq_i_this = np.append(obs_o3_cmaq_i_this, obs_o3_cmaq_i_val)
			obs_o3_cmaq_j_this = np.append(obs_o3_cmaq_j_this, obs_o3_cmaq_j_val)

			n_pm25_this = len(obs_pm25_sid_this)
			n_o3_this = len(obs_o3_sid_this)

			## Third, step through the obs from this time
			for ii in range(n_pm25_this):
				if np.isnan(obs_pm25_sid_this[ii]):
					## If the element is a nan, then skip to the next iteration
					continue
				else:
					## Is this station ID already in the array?
					if (obs_pm25_sid_this[ii] in obs_pm25_sid_all):
						## Increment the counter for this site
						ind = np.where(obs_pm25_sid_all == obs_pm25_sid_this[ii])
						obs_pm25_cnt_all[ind] += 1

						## Add the concentration for this time for this site
						obs_pm25_con_all[tt,ind] = obs_pm25_con_this[ii]
					else:
						## Append a new element to the arrays for this new site
						obs_pm25_sid_all = np.append(obs_pm25_sid_all, obs_pm25_sid_this[ii])
						obs_pm25_lon_all = np.append(obs_pm25_lon_all, obs_pm25_lon_this[ii])
						obs_pm25_lat_all = np.append(obs_pm25_lat_all, obs_pm25_lat_this[ii])
						obs_pm25_cnt_all = np.append(obs_pm25_cnt_all, 1)
						obs_pm25_cmaq_i_all = np.append(obs_pm25_cmaq_i_all, obs_pm25_cmaq_i_this[ii])
						obs_pm25_cmaq_j_all = np.append(obs_pm25_cmaq_j_all, obs_pm25_cmaq_j_this[ii])

						## Add a new site to the concentration array, then add the concentration for this time
						obs_pm25_con_all = np.append(obs_pm25_con_all, np.full([tt+1,1], np.nan), axis=1)
						obs_pm25_con_all[tt,-1] = obs_pm25_con_this[ii]

			for ii in range(n_o3_this):
				if np.isnan(obs_o3_sid_this[ii]):
					## If the element is a nan, then skip to the next iteration
					continue
				else:
					## Is this station ID already in the array?
					if (obs_o3_sid_this[ii] in obs_o3_sid_all):
						## Increment the counter for this site
						ind = np.where(obs_o3_sid_all == obs_o3_sid_this[ii])
						obs_o3_cnt_all[ind] += 1

						## Add the concentration for this time for this site
						obs_o3_con_all[tt,ind] = obs_o3_con_this[ii]
					else:
						## Append a new element to the arrays for this new site
						obs_o3_sid_all = np.append(obs_o3_sid_all, obs_o3_sid_this[ii])
						obs_o3_lon_all = np.append(obs_o3_lon_all, obs_o3_lon_this[ii])
						obs_o3_lat_all = np.append(obs_o3_lat_all, obs_o3_lat_this[ii])
						obs_o3_cnt_all = np.append(obs_o3_cnt_all, 1)
						obs_o3_cmaq_i_all = np.append(obs_o3_cmaq_i_all, obs_o3_cmaq_i_this[ii])
						obs_o3_cmaq_j_all = np.append(obs_o3_cmaq_j_all, obs_o3_cmaq_j_this[ii])

						# Add a new site to the concentration array, then add the concentration for this time
						obs_o3_con_all = np.append(obs_o3_con_all, np.full([tt+1,1], np.nan), axis=1)
						obs_o3_con_all[tt,-1] = obs_o3_con_this[ii]

	## Convert raw counts into percentages
	obs_pm25_pct_all = obs_pm25_cnt_all / n_times * 100.0
	obs_o3_pct_all = obs_o3_cnt_all / n_times * 100.0

	## Get total number of stations with any valid reports
	n_all_pm25 = len(obs_pm25_pct_all)
	n_all_o3 = len(obs_o3_pct_all)

	n_stn_pm25_gt80 = (obs_pm25_pct_all >= 80.0).sum()
	n_stn_o3_gt80 = (obs_o3_pct_all >= 80.0).sum()

	print('Number (percentage) of stations with any valid reports reporting at least 80% of the time:')
	print('PM2.5: '+str(n_stn_pm25_gt80)+' ('+str(100.0*n_stn_pm25_gt80/n_all_pm25)+'%)')
	print('   O3: '+str(n_stn_o3_gt80)+' ('+str(100.0*n_stn_o3_gt80/n_all_o3)+'%)')

	## Write file with observation counts & percentages
	obs_times = xr.DataArray(dt_array,
						coords={'n_obs_times':n_times}, dims=['n_times'],
						attrs={'description':'valid dates/times of AirNow observations'})

	pm25_all_sid = xr.DataArray(obs_pm25_sid_all,
							coords={'n_obs_all_pm25':n_all_pm25}, dims=['n_all_pm25'],
							attrs={'description':'AirNow PM2.5 station ID', '_FillValue':-999})
	pm25_all_lon = xr.DataArray(obs_pm25_lon_all,
							coords={'n_obs_all_pm25':n_all_pm25}, dims=['n_all_pm25'],
							attrs={'description':'AirNow PM2.5 station longitude', 'units':'deg E'})
	pm25_all_lat = xr.DataArray(obs_pm25_lat_all,
							coords={'n_obs_all_pm25':n_all_pm25}, dims=['n_all_pm25'],
							attrs={'description':'AirNow PM2.5 station latitude', 'units':'deg N'})
	pm25_all_cnt = xr.DataArray(obs_pm25_cnt_all,
							coords={'n_obs_all_pm25':n_all_pm25}, dims=['n_all_pm25'],
							attrs={'description':'AirNow PM2.5 station count of valid reports'})
	pm25_all_pct = xr.DataArray(obs_pm25_pct_all,
							coords={'n_obs_all_pm25':n_all_pm25}, dims=['n_all_pm25'],
							attrs={'description':'AirNow PM2.5 station percent of hours with valid reports','units':'%'})
	pm25_all_con = xr.DataArray(obs_pm25_con_all,
							coords={'n_obs_times':n_times, 'n_obs_all_pm25':n_all_pm25},
							dims=['n_times', 'n_all_pm25'],
							attrs={'description':'AirNow PM2.5 concentration observations', 'units':'ug/m3'})
	pm25_all_cmaq_i = xr.DataArray(obs_pm25_cmaq_i_all,
							coords={'n_obs_all_pm25':n_all_pm25}, dims=['n_all_pm25'],
							attrs={'description':'AirNow PM2.5 station nearest CMAQ grid i-index'})
	pm25_all_cmaq_j = xr.DataArray(obs_pm25_cmaq_j_all,
							coords={'n_obs_all_pm25':n_all_pm25}, dims=['n_all_pm25'],
							attrs={'description':'AirNow PM2.5 station nearest CMAQ grid j-index'})

	o3_all_sid = xr.DataArray(obs_o3_sid_all,
							coords={'n_obs_all_o3':n_all_o3}, dims=['n_all_o3'],
							attrs={'description':'AirNow O3 station ID', '_FillValue':-999})
	o3_all_lon = xr.DataArray(obs_o3_lon_all,
							coords={'n_obs_all_o3':n_all_o3}, dims=['n_all_o3'],
							attrs={'description':'AirNow O3 station longitude', 'units':'deg E'})
	o3_all_lat = xr.DataArray(obs_o3_lat_all,
							coords={'n_obs_all_o3':n_all_o3}, dims=['n_all_o3'],
							attrs={'description':'AirNow O3 station latitude', 'units':'deg N'})
	o3_all_cnt = xr.DataArray(obs_o3_cnt_all,
							coords={'n_obs_all_o3':n_all_o3}, dims=['n_all_o3'],
							attrs={'description':'AirNow O3 station count of valid reports'})
	o3_all_pct = xr.DataArray(obs_o3_pct_all,
							coords={'n_obs_all_o3':n_all_o3}, dims=['n_all_o3'],
							attrs={'description':'AirNow O3 station percent of hours with valid reports','units':'%'})
	o3_all_con = xr.DataArray(obs_o3_con_all,
							coords={'n_obs_times':n_times, 'n_obs_all_o3':n_all_o3},
							dims=['n_times', 'n_all_o3'],
							attrs={'description':'AirNow O3 concentration observations', 'units':'ppbv'})
	o3_all_cmaq_i = xr.DataArray(obs_o3_cmaq_i_all,
							coords={'n_obs_all_o3':n_all_o3}, dims=['n_all_o3'],
							attrs={'description':'AirNow O3 station nearest CMAQ grid i-index'})
	o3_all_cmaq_j = xr.DataArray(obs_o3_cmaq_j_all,
							coords={'n_obs_all_o3':n_all_o3}, dims=['n_all_o3'],
							attrs={'description':'AirNow O3 station nearest CMAQ grid j-index'})

	ds_all = xr.Dataset(
					data_vars={ 'pm25_sid_all':pm25_all_sid, 'o3_sid_all':o3_all_sid,
									'pm25_lon_all':pm25_all_lon, 'o3_lon_all':o3_all_lon,
									'pm25_lat_all':pm25_all_lat, 'o3_lat_all':o3_all_lat,
									'pm25_cnt_all':pm25_all_cnt, 'o3_cnt_all':o3_all_cnt,
									'pm25_pct_all':pm25_all_pct, 'o3_pct_all':o3_all_pct,
									'pm25_con_all':pm25_all_con, 'o3_con_all':o3_all_con,
									'pm25_cmaq_i_all':pm25_all_cmaq_i, 'o3_cmaq_i_all':o3_all_cmaq_i,
									'pm25_cmaq_j_all':pm25_all_cmaq_j, 'o3_cmaq_j_all':o3_all_cmaq_j,
									'obs_times':obs_times,
						},
					coords={'n_obs_all_pm25':n_all_pm25, 'n_obs_all_o3':n_all_o3, 'n_obs_times':n_times},
					attrs={'description':'Raw counts and percent uptime for all PM2.5 & O3 AirNow observations from '+date_beg+' to '+date_end},
					)

	fname_out = out_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_all.nc')
	print('Writing file '+str(fname_out))
	ds_all.to_netcdf(fname_out)
	
	## Everything from here and below could also be done in a separate script that reads in the file just written.
	## But processing 1-2 years of data still only takes a few minutes, so it's not a big deal or time-saver.

	## Randomly select a fraction of stations to withhold for AnEn training that are above uptime threshold
	## First, find the indices of stations that exceed the uptime threshold
	inds_uptime_pm25 = np.where(obs_pm25_pct_all >= up_thresh_pm25)[0]
	inds_uptime_o3 = np.where(obs_o3_pct_all >= up_thresh_o3)[0]
	n_uptime_pm25 = len(inds_uptime_pm25)
	n_uptime_o3 = len(inds_uptime_o3)

	## Second, randomly select indices from that group to be used for AnEn training
	n_use_pm25 = int(n_uptime_pm25 * use_fraction)
	n_use_o3 = int(n_uptime_o3 * use_fraction)

	np.random.seed(int(dt_beg.strftime(fmt_yyyymmddhh))) ;
	inds_uptime_use_pm25 = np.sort(np.random.choice(range(n_uptime_pm25), n_use_pm25, replace=False))
	inds_use_pm25 = np.take(inds_uptime_pm25, inds_uptime_use_pm25)

	np.random.seed(int(dt_beg.strftime(fmt_yyyymmddhh))) ;
	inds_uptime_use_o3 = np.sort(np.random.choice(range(n_uptime_o3), n_use_o3, replace=False))
	inds_use_o3 = np.take(inds_uptime_o3, inds_uptime_use_o3)

	inds_val_pm25 = np.delete(range(n_all_pm25), inds_use_pm25)
	inds_val_o3 = np.delete(range(n_all_o3), inds_use_o3)
	n_val_pm25 = len(inds_val_pm25)
	n_val_o3 = len(inds_val_o3)

	print('\npm2.5 stations:')
	print('all = '+str(n_all_pm25)+', use = '+str(n_use_pm25)+', val = '+str(n_val_pm25))
	print('\no3 stations:')
	print('all = '+str(n_all_o3)+', use = '+str(n_use_o3)+', val = '+str(n_val_o3))
	print('\n')

	obs_pm25_sid_use = obs_pm25_sid_all[inds_use_pm25]
	obs_pm25_lon_use = obs_pm25_lon_all[inds_use_pm25]
	obs_pm25_lat_use = obs_pm25_lat_all[inds_use_pm25]
	obs_pm25_cnt_use = obs_pm25_cnt_all[inds_use_pm25]
	obs_pm25_pct_use = obs_pm25_pct_all[inds_use_pm25]
	obs_pm25_con_use = obs_pm25_con_all[:,inds_use_pm25]
	obs_pm25_cmaq_i_use = obs_pm25_cmaq_i_all[inds_use_pm25]
	obs_pm25_cmaq_j_use = obs_pm25_cmaq_j_all[inds_use_pm25]

	obs_pm25_sid_val = obs_pm25_sid_all[inds_val_pm25]
	obs_pm25_lon_val = obs_pm25_lon_all[inds_val_pm25]
	obs_pm25_lat_val = obs_pm25_lat_all[inds_val_pm25]
	obs_pm25_cnt_val = obs_pm25_cnt_all[inds_val_pm25]
	obs_pm25_pct_val = obs_pm25_pct_all[inds_val_pm25]
	obs_pm25_con_val = obs_pm25_con_all[:,inds_val_pm25]
	obs_pm25_cmaq_i_val = obs_pm25_cmaq_i_all[inds_val_pm25]
	obs_pm25_cmaq_j_val = obs_pm25_cmaq_j_all[inds_val_pm25]

	obs_o3_sid_use = obs_o3_sid_all[inds_use_o3]
	obs_o3_lon_use = obs_o3_lon_all[inds_use_o3]
	obs_o3_lat_use = obs_o3_lat_all[inds_use_o3]
	obs_o3_cnt_use = obs_o3_cnt_all[inds_use_o3]
	obs_o3_pct_use = obs_o3_pct_all[inds_use_o3]
	obs_o3_con_use = obs_o3_con_all[:,inds_use_o3]
	obs_o3_cmaq_i_use = obs_o3_cmaq_i_all[inds_use_o3]
	obs_o3_cmaq_j_use = obs_o3_cmaq_j_all[inds_use_o3]

	obs_o3_sid_val = obs_o3_sid_all[inds_val_o3]
	obs_o3_lon_val = obs_o3_lon_all[inds_val_o3]
	obs_o3_lat_val = obs_o3_lat_all[inds_val_o3]
	obs_o3_cnt_val = obs_o3_cnt_all[inds_val_o3]
	obs_o3_pct_val = obs_o3_pct_all[inds_val_o3]
	obs_o3_con_val = obs_o3_con_all[:,inds_val_o3]
	obs_o3_cmaq_i_val = obs_o3_cmaq_i_all[inds_val_o3]
	obs_o3_cmaq_j_val = obs_o3_cmaq_j_all[inds_val_o3]

	## Third, write files with observation counts & percentages to give a permanent reference for use/val sites
	pm25_use_sid = xr.DataArray(obs_pm25_sid_use,
							coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
							attrs={'description':'AirNow PM2.5 station ID', '_FillValue':-999})
	pm25_use_lon = xr.DataArray(obs_pm25_lon_use,
							coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
							attrs={'description':'AirNow PM2.5 station longitude', 'units':'deg E'})
	pm25_use_lat = xr.DataArray(obs_pm25_lat_use,
							coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
							attrs={'description':'AirNow PM2.5 station latitude', 'units':'deg N'})
	pm25_use_cnt = xr.DataArray(obs_pm25_cnt_use,
							coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
							attrs={'description':'AirNow PM2.5 station count of valid reports'})
	pm25_use_pct = xr.DataArray(obs_pm25_pct_use,
							coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
							attrs={'description':'AirNow PM2.5 station percent of hours with valid reports','units':'%'})
	pm25_use_con = xr.DataArray(obs_pm25_con_use,
							coords={'n_obs_times':n_times, 'n_obs_use_pm25':n_use_pm25},
							dims=['n_times', 'n_use_pm25'],
							attrs={'description':'AirNow PM2.5 concentration observations', 'units':'ug/m3'})
	pm25_use_cmaq_i = xr.DataArray(obs_pm25_cmaq_i_use,
							coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
							attrs={'description':'AirNow PM2.5 station nearest CMAQ grid i-index'})
	pm25_use_cmaq_j = xr.DataArray(obs_pm25_cmaq_j_use,
							coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
							attrs={'description':'AirNow PM2.5 station nearest CMAQ grid j-index'})

	o3_use_sid = xr.DataArray(obs_o3_sid_use,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station ID', '_FillValue':-999})
	o3_use_lon = xr.DataArray(obs_o3_lon_use,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station longitude', 'units':'deg E'})
	o3_use_lat = xr.DataArray(obs_o3_lat_use,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station latitude', 'units':'deg N'})
	o3_use_cnt = xr.DataArray(obs_o3_cnt_use,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station count of valid reports'})
	o3_use_pct = xr.DataArray(obs_o3_pct_use,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station percent of hours with valid reports','units':'%'})
	o3_use_con = xr.DataArray(obs_o3_con_use,
							coords={'n_obs_times':n_times, 'n_obs_use_o3':n_use_o3},
							dims=['n_times', 'n_use_o3'],
							attrs={'description':'AirNow O3 concentration observations', 'units':'ppbv'})
	o3_use_cmaq_i = xr.DataArray(obs_o3_cmaq_i_use,
							coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station nearest CMAQ grid i-index'})
	o3_use_cmaq_j = xr.DataArray(obs_o3_cmaq_j_use,
							coords={'n_obs_ise_o3':n_use_o3}, dims=['n_use_o3'],
							attrs={'description':'AirNow O3 station nearest CMAQ grid j-index'})

	pm25_val_sid = xr.DataArray(obs_pm25_sid_val,
							coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
							attrs={'description':'AirNow PM2.5 station ID', '_FillValue':-999})
	pm25_val_lon = xr.DataArray(obs_pm25_lon_val,
							coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
							attrs={'description':'AirNow PM2.5 station longitude', 'units':'deg E'})
	pm25_val_lat = xr.DataArray(obs_pm25_lat_val,
							coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
							attrs={'description':'AirNow PM2.5 station latitude', 'units':'deg N'})
	pm25_val_cnt = xr.DataArray(obs_pm25_cnt_val,
							coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
							attrs={'description':'AirNow PM2.5 station count of valid reports'})
	pm25_val_pct = xr.DataArray(obs_pm25_pct_val,
							coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
							attrs={'description':'AirNow PM2.5 station percent of hours with valid reports','units':'%'})
	pm25_val_con = xr.DataArray(obs_pm25_con_val,
							coords={'n_obs_times':n_times, 'n_obs_val_pm25':n_val_pm25},
							dims=['n_times', 'n_val_pm25'],
							attrs={'description':'AirNow PM2.5 concentration observations', 'units':'ug/m3'})
	pm25_val_cmaq_i = xr.DataArray(obs_pm25_cmaq_i_val,
							coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
							attrs={'description':'AirNow PM2.5 station nearest CMAQ grid i-index'})
	pm25_val_cmaq_j = xr.DataArray(obs_pm25_cmaq_j_val,
							coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
							attrs={'description':'AirNow PM2.5 station nearest CMAQ grid j-index'})

	o3_val_sid = xr.DataArray(obs_o3_sid_val,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station ID', '_FillValue':-999})
	o3_val_lon = xr.DataArray(obs_o3_lon_val,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station longitude', 'units':'deg E'})
	o3_val_lat = xr.DataArray(obs_o3_lat_val,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station latitude', 'units':'deg N'})
	o3_val_cnt = xr.DataArray(obs_o3_cnt_val,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station count of valid reports'})
	o3_val_pct = xr.DataArray(obs_o3_pct_val,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station percent of hours with valid reports','units':'%'})
	o3_val_con = xr.DataArray(obs_o3_con_val,
							coords={'n_obs_times':n_times, 'n_obs_val_o3':n_val_o3},
							dims=['n_times', 'n_val_o3'],
							attrs={'description':'AirNow O3 concentration observations', 'units':'ppbv'})
	o3_val_cmaq_i = xr.DataArray(obs_o3_cmaq_i_val,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station nearest CMAQ grid i-index'})
	o3_val_cmaq_j = xr.DataArray(obs_o3_cmaq_j_val,
							coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
							attrs={'description':'AirNow O3 station nearest CMAQ grid j-index'})

	ds_use = xr.Dataset(
					data_vars={ 'pm25_sid_use':pm25_use_sid, 'o3_sid_use':o3_use_sid,
									'pm25_lon_use':pm25_use_lon, 'o3_lon_use':o3_use_lon,
									'pm25_lat_use':pm25_use_lat, 'o3_lat_use':o3_use_lat,
									'pm25_cnt_use':pm25_use_cnt, 'o3_cnt_use':o3_use_cnt,
									'pm25_pct_use':pm25_use_pct, 'o3_pct_use':o3_use_pct,
									'pm25_con_use':pm25_use_con, 'o3_con_use':o3_use_con,
									'pm25_cmaq_i_use':pm25_use_cmaq_i, 'o3_cmaq_i_use':o3_use_cmaq_i,
									'pm25_cmaq_j_use':pm25_use_cmaq_j, 'o3_cmaq_j_use':o3_use_cmaq_j,
									'obs_times':obs_times,
						},
					coords={'n_obs_use_pm25':n_use_pm25, 'n_obs_use_o3':n_use_o3, 'n_obs_times':n_times},
					attrs={'description':'Raw counts and percent uptime for the PM2.5 & O3 AirNow observations from '+date_beg+' to '+date_end+' that are used for AnEn training'},
					)

	ds_val = xr.Dataset(
					data_vars={ 'pm25_sid_val':pm25_val_sid, 'o3_sid_val':o3_val_sid,
									'pm25_lon_val':pm25_val_lon, 'o3_lon_val':o3_val_lon,
									'pm25_lat_val':pm25_val_lat, 'o3_lat_val':o3_val_lat,
									'pm25_cnt_val':pm25_val_cnt, 'o3_cnt_val':o3_val_cnt,
									'pm25_pct_val':pm25_val_pct, 'o3_pct_val':o3_val_pct,
									'pm25_con_val':pm25_val_con, 'o3_con_val':o3_val_con,
									'pm25_cmaq_i_val':pm25_val_cmaq_i, 'o3_cmaq_i_val':o3_val_cmaq_i,
									'pm25_cmaq_j_val':pm25_val_cmaq_j, 'o3_cmaq_j_val':o3_val_cmaq_j,
									'obs_times':obs_times
						},
					coords={'n_obs_val_pm25':n_val_pm25, 'n_obs_val_o3':n_val_o3, 'n_obs_times':n_times},
					attrs={'description':'Raw counts and percent uptime for the PM2.5 & O3 AirNow observations from '+date_beg+' to '+date_end+' that are withheld for AnEn validation'},
					)

	fname_use = out_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_use.nc')
	fname_val = out_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_val.nc')
	print('Writing file '+str(fname_use))
	ds_use.to_netcdf(fname_use)
	print('Writing file '+str(fname_val))
	ds_val.to_netcdf(fname_val)

	## ================
	## PLOTTING SECTION
	## ================

	suptitle_pm25_all = 'AirNow '+pm25_str+' Obs, All Sites ('+str(n_all_pm25)+' Stations)'
	suptitle_o3_all   = 'AirNow '+o3_str+' Obs, All Sites ('+str(n_all_o3)+' Stations)'
	suptitle_pm25_use = 'AirNow '+pm25_str+' Obs, Sites Used for Training ('+str(n_use_pm25)+' Stations)'
	suptitle_o3_use   = 'AirNow '+o3_str+' Obs, Sites Used for Training ('+str(n_use_o3)+' Stations)'
	suptitle_pm25_val = 'AirNow '+pm25_str+' Obs, Sites Held for Validation ('+str(n_val_pm25)+' Stations)'
	suptitle_o3_val   = 'AirNow '+o3_str+' Obs, Sites Held for Validation ('+str(n_val_o3)+' Stations)'

	if plot_hist:
		fontsize = 12
		suptitle_y = 0.95
		mpl.rcParams['figure.figsize'] = (10,8)
		mpl.rcParams['figure.titlesize'] = fontsize+4
		mpl.rcParams['grid.color'] = 'gray'
		mpl.rcParams['grid.linestyle'] = ':'
		mpl.rcParams['font.size'] = fontsize
		mpl.rcParams['savefig.bbox'] = 'tight'
		mpl.rcParams['hist.bins'] = np.arange(0.0, 100.1, 5.0)

		xlabel = 'Percentage of Hours with Valid Reports'
		ylabel = 'Count'
		title = date_beg+' to '+date_end

		## Make histograms with percentage bins for all stations within the period
		fname = out_dir.joinpath('hist_pm2.5_valid_obs_'+date_beg+'-'+date_end+'_all.'+plot_type)
		var = obs_pm25_pct_all
		make_hist_plot(fname, var, suptitle_pm25_all, title, xlabel, ylabel)

		fname = out_dir.joinpath('hist_o3_valid_obs_'+date_beg+'-'+date_end+'_all.'+plot_type)
		var = obs_o3_pct_all
		make_hist_plot(fname, var, suptitle_o3_all, title, xlabel, ylabel)

		## Make histograms with percentage bins for use stations within the period
		fname = out_dir.joinpath('hist_pm2.5_valid_obs_'+date_beg+'-'+date_end+'_use.'+plot_type)
		var = obs_pm25_pct_use
		make_hist_plot(fname, var, suptitle_pm25_use, title, xlabel, ylabel)

		fname = out_dir.joinpath('hist_o3_valid_obs_'+date_beg+'-'+date_end+'_use.'+plot_type)
		var = obs_o3_pct_use
		make_hist_plot(fname, var, suptitle_o3_use, title, xlabel, ylabel)

		## Make histograms with percentage bins for val stations within the period
		fname = out_dir.joinpath('hist_pm2.5_valid_obs_'+date_beg+'-'+date_end+'_val.'+plot_type)
		var = obs_pm25_pct_val
		make_hist_plot(fname, var, suptitle_pm25_val, title, xlabel, ylabel)

		fname = out_dir.joinpath('hist_o3_valid_obs_'+date_beg+'-'+date_end+'_val.'+plot_type)
		var = obs_o3_pct_val
		make_hist_plot(fname, var, suptitle_o3_val, title, xlabel, ylabel)

	if plot_maps:
		fontsize = 12
		suptitle_y = 0.88

		## Make maps with percentage valid reports for all stations within the period
		## First, open the CMAQ	sample/coordinate file to get map projection parameters
		cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
		print('Reading CMAQ coordinate data from '+str(cmaq_fname))
		cmaq_ds = xr.open_dataset(cmaq_fname)
		cmaq_lon = cmaq_ds.LON[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'}) # range (-180, 180]
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

		'''
		## Or use a CAMS file regridded to CMAQ grid
		cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst','2020','08','20200801_00','cams_o3_pm25_regrid_cmaq_20200801_0000.nc')
		print('Reading CMAQ coordinate data from '+str(cmaq_fname))
		cmaq_ds = xr.open_dataset(cmaq_fname)
		cmaq_lat = cmaq_ds.latitude
		cmaq_lon = cmaq_ds.longitude
		cmaq_lat = cmaq_lat.rename({'latitude':'XLAT', 'longitude':'XLONG'})
		cmaq_lon = cmaq_lon.rename({'latitude':'XLAT', 'longitude':'XLONG'})
		map_proj = ds_xr.attrs['MAP_PROJ']
		truelat1 = ds_xr.attrs['TRUELAT1']
		truelat2 = ds_xr.attrs['TRUELAT2']
		stand_lon = ds_xr.attrs['STAND_LON']
		cen_lat = ds_xr.attrs['CEN_LAT']
		cen_lon = ds_xr.attrs['CEN_LON']
		pole_lat = ds_xr.attrs['POLE_LAT']
		pole_lon = ds_xr.attrs['POLE_LON']
		moad_cen_lat = ds_xr.attrs['MOAD_CEN_LAT']
		dx = ds_xr.attrs['DX']
		dy = ds_xr.attrs['DY']
		'''

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

		title = date_beg+' to '+date_end
		cbar_lab = 'Percentage of Hours with Valid Reports'
		extend = 'max'
		cmap = mpl.cm.viridis
		bounds = np.arange(0.0, 95.1, 5.0)
		norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)

		## Make maps showing all stations in the period
		suptitle = suptitle_pm25_all
		fname = out_dir.joinpath('map_pm2.5_valid_obs_'+date_beg+'-'+date_end+'_all.'+plot_type)
		var = obs_pm25_pct_all
		map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
			cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
			borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
			marker_lat=obs_pm25_lat_all, marker_lon=obs_pm25_lon_all, marker='o', marker_val_fill=True, marker_size=64,
			title_c=title, suptitle_y=suptitle_y)

		suptitle = suptitle_o3_all
		fname = out_dir.joinpath('map_o3_valid_obs_'+date_beg+'-'+date_end+'_all.'+plot_type)
		var = obs_o3_pct_all
		map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
			cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
			borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
			marker_lat=obs_o3_lat_all, marker_lon=obs_o3_lon_all, marker='o', marker_val_fill=True, marker_size=64,
			title_c=title, suptitle_y=suptitle_y)

		## Make maps showing use stations in the period
		suptitle = suptitle_pm25_use
		fname = out_dir.joinpath('map_pm2.5_valid_obs_'+date_beg+'-'+date_end+'_use.'+plot_type)
		var = obs_pm25_pct_use
		map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
			cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
			borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
			marker_lat=obs_pm25_lat_use, marker_lon=obs_pm25_lon_use, marker='o', marker_val_fill=True, marker_size=64,
			title_c=title, suptitle_y=suptitle_y)

		suptitle = suptitle_o3_use
		fname = out_dir.joinpath('map_o3_valid_obs_'+date_beg+'-'+date_end+'_use.'+plot_type)
		var = obs_o3_pct_use
		map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
			cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
			borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
			marker_lat=obs_o3_lat_use, marker_lon=obs_o3_lon_use, marker='o', marker_val_fill=True, marker_size=64,
			title_c=title, suptitle_y=suptitle_y)

		## Make maps showing val stations in the period
		suptitle = suptitle_pm25_val
		fname = out_dir.joinpath('map_pm2.5_valid_obs_'+date_beg+'-'+date_end+'_val.'+plot_type)
		var = obs_pm25_pct_val
		map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
			cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
			borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
			marker_lat=obs_pm25_lat_val, marker_lon=obs_pm25_lon_val, marker='o', marker_val_fill=True, marker_size=64,
			title_c=title, suptitle_y=suptitle_y)

		suptitle = suptitle_o3_val
		fname = out_dir.joinpath('map_o3_valid_obs_'+date_beg+'-'+date_end+'_val.'+plot_type)
		var = obs_o3_pct_val
		map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
			cmaq_lon, cmaq_lat, cmap, bounds, norm, extend,
			borders=borders, states=states, lakes=lakes, oceans=oceans, water_color='lightblue',
			marker_lat=obs_o3_lat_val, marker_lon=obs_o3_lon_val, marker='o', marker_val_fill=True, marker_size=64,
			title_c=title, suptitle_y=suptitle_y)



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
