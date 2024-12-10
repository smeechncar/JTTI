'''
build_cams_at_airnow_nc_files.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 10 Dec 2021

This script reads in CSV files of CAMS O3 and PM2.5 data at AirNOW sites (from extract_cams_at_airnow_sites.py).
It is assumed that these CSV files are written out at the same master list of sites.
CAMS forecast data has 1-h resolution for PM2.5 and 3-h resolution for O3.
CAMS reanalysis data has 3-h resolution for both PM2.5 and O3.
This script interpolates CAMS data to 1-h resolution (if needed) to assemble a time series.
This time series for the entire time range and a master site list is written, one file for O3 and one file for PM2.5.
'''

import os
import sys
import pathlib
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4
import scipy as sp
import scipy.interpolate

def nan_helper(y):
	'''
	Helper to handle indices and logical indices of Nans.
	See: https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array

	Input:
	- y: 1d numpy array with possible NaNs
	Output:
	- nans: logical indices of NaNs
	- index: a function, with signature indices= index(logical_indices)
	  to convert logical indices of NaNs to 'equivalent' indices
	Example:
	>> # linear interpolation of NaNs
	>> nans, x = nan_helper(y)
	>> y[nans] = np.interp(x(nans), x(~nans), y[~nans])
	'''
	return np.isnan(y), lambda z: z.nonzero()[0]

## Main program:

beg_date = '2020-08-01_00'
end_date = '2022-01-01_00'

cams_fcst = True
cams_eac4 = False
if cams_fcst and cams_eac4:
	print('ERROR: Both cams_fcst and cams_eac4 set to True. Please only set one of them to True.')
	sys.exit()
if not cams_fcst and not cams_eac4:
	print('ERROR: Both cams_fcst and cams_eac4 set to False. Please set one of them to True.')
	sys.exit()

fmt_date = '%Y-%m-%d'
fmt_date_hh = '%Y-%m-%d_%H'
fmt_yyyymmdd = '%Y%m%d'
fmt_yyyymmdd_hh = '%Y%m%d_%H'

## Build array of dates
beg_dt = pd.to_datetime(beg_date, format=fmt_date_hh)
end_dt = pd.to_datetime(end_date, format=fmt_date_hh)
all_dt = pd.date_range(start=beg_dt, end=end_dt, freq='1H')
n_times = len(all_dt)
all_hrs = np.arange(0, n_times, 1)

## Get master list of stations (id, lat, lon)
sites_pm25	= pathlib.Path('/','glade','campaign','ral','nsap','JTTI','AIRNOW_MASTER_COORDI_PM25.csv')
sites_o3		= pathlib.Path('/','glade','campaign','ral','nsap','JTTI','AIRNOW_MASTER_COORDI_O3.csv')
if cams_fcst:
	out_dir	= pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')
elif cams_eac4:
	out_dir	= pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4')

df_pm25 = pd.read_csv(sites_pm25, header=0)
pm25_site_ids = np.asarray(df_pm25['Site_ID'])
pm25_site_lon = np.asarray(df_pm25[' Site_Longitude'])
pm25_site_lat = np.asarray(df_pm25[' Site_latitude'])
df_o3 = pd.read_csv(sites_o3, header=0)
o3_site_ids = np.asarray(df_o3['Site_ID'])
o3_site_lon = np.asarray(df_o3[' Site_Longitude'])
o3_site_lat = np.asarray(df_o3[' Site_latitude'])

n_sites_pm25 = len(pm25_site_ids)
n_sites_o3   = len(o3_site_ids)

## Initialize arrays
ts_cams_pm25 = np.full([n_sites_pm25, n_times], np.nan)
ts_cams_o3   = np.full([n_sites_o3, n_times], np.nan)

## Fill the time series arrays by reading through the CSV files
for tt in range(n_times):
	this_dt = all_dt[tt]
	this_yr = this_dt.strftime('%Y')
	this_mo = this_dt.strftime('%m')
	this_dy = this_dt.strftime('%d')
	this_hr = this_dt.strftime('%H')
	this_yyyymmdd = this_dt.strftime(fmt_yyyymmdd)
	this_date_str = this_dt.strftime(fmt_date)
	this_date_hh_str = this_dt.strftime(fmt_date_hh)
	if int(this_hr) < 12:
		this_cycle = this_yyyymmdd+'_00'
	else:
		this_cycle = this_yyyymmdd+'_12'

	if cams_fcst:
		this_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst',this_yr,this_mo,this_cycle)
	elif cams_eac4:
		this_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4',this_yr,this_mo)

	## If file exists, open it
	pm25_fname = this_dir.joinpath('cams_pm2.5_at_airnow_'+this_date_hh_str+'00.csv')
	o3_fname   = this_dir.joinpath('cams_o3sfc_at_airnow_'+this_date_hh_str+'00.csv')

	if o3_fname.is_file():
		print('Reading file '+str(o3_fname))
		df_o3 = pd.read_csv(o3_fname)
		o3_vals = df_o3['cams_o3_ppb'].values
		ts_cams_o3[:,tt] = o3_vals

	if pm25_fname.is_file():
		print('Reading file '+str(pm25_fname))
		df_pm25 = pd.read_csv(pm25_fname)
		pm25_vals = df_pm25['cams_pm25_ug'].values
		ts_cams_pm25[:,tt] = pm25_vals

## Interpolate the ozone data from 3-hourly to 1-hourly
for ss in range(n_sites_o3):
	nans, x = nan_helper(ts_cams_o3[ss,:])
	ts_cams_o3[ss,nans] = np.interp(x(nans), x(~nans), ts_cams_o3[ss,~nans])

if cams_eac4:
	for ss in range(n_sites_pm25):
		nans, x = nan_helper(ts_cams_pm25[ss,:])
		ts_cams_pm25[ss,nans] = np.interp(x(nans), x(~nans), ts_cams_pm25[ss,~nans])

## Build an xarray dataset
if cams_fcst:
	attrs_desc = 'Time series of CAMS forecast O3 and PM2.5 values'
elif cams_eac4:
	attrs_desc = 'Time series of CAMS EAC4 reanalysis O3 and PM2.5 values'

ds = xr.Dataset(
	data_vars = dict(
		o3_cams = (['o3_sites', 'time'], ts_cams_o3),
		pm25_cams = (['pm25_sites', 'time'], ts_cams_pm25),
		o3_site_ids = (['o3_sites'], o3_site_ids),
		pm25_site_ids = (['pm25_sites'], pm25_site_ids),
	),
	coords = dict(
		o3_lon = (['o3_sites'], o3_site_lon),
		o3_lat = (['o3_sites'], o3_site_lat),
		pm25_lon = (['pm25_sites'], pm25_site_lon),
		pm25_lat = (['pm25_sites'], pm25_site_lat),
		time = all_dt
	),
	attrs = dict(description=attrs_desc),
)

## Write the data to a netcdf file
beg_date_str = beg_dt.strftime(fmt_yyyymmdd)
end_date_str = end_dt.strftime(fmt_yyyymmdd)
if cams_fcst:
	fname = out_dir.joinpath('cams_fcst_o3_pm25_'+beg_date_str+'_to_'+end_date_str+'.nc')
elif cams_eac4:
	fname = out_dir.joinpath('cams_eac4_o3_pm25_'+beg_date_str+'_to_'+end_date_str+'.nc')
print('Saving data to '+str(fname))
ds.to_netcdf(fname)

## Verify that the data as written is identical
#ds_new = xr.open_dataset(fname)
#print(ds_new)
#time_new = ds_new.time
#print(all_dt)
#print(time_new)
	
