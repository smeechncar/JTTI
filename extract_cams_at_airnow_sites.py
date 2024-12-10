'''
extract_cams_at_airnow_sites.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 9 Dec 2021

This script first reads in a coordinates file to obtain the lat/lon locations of AirNOW stations.
Second, it reads in either a CAMS reanalysis (EAC4) or forecast file to get the grid coordinates of those stations.
Finally, for each forecast/reanalysis file (one per valid time), a CSV file of CAMS values at AirNOW stations
is written out.
A subsequent script will then read in these CSV files and then generate a NetCDF file over a desired time span.

Edits:
7 Feb 2022 by JAL
 - Add optional plotting routine for the CAMS fields.
 - Fixed mismatch of CAMS longitude [0.0, 360.0) and AIRNow longitude (-180.0, 180.0] that caused squirrelly
   behavior in the spline interpolation. Also changed the spline from a degree-3 to degree-1.
'''

import os
import sys
import pathlib
import datetime as dt
import pandas as pd
import numpy as np
import xarray as xr
import scipy as sp
import scipy.interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
from functions import map_functions

beg_date = '2020-08-13_00'
end_date = '2020-08-15_23'

cams_fcst = False
cams_eac4 = True

plot_o3   = True
plot_pm25 = True
plot_type = 'png'

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
fmt_ddmmyyyyhhmm = '%d %b %Y/%H%M'

mpl_pm25 = '$\mathregular{PM_{2.5}}$'
mpl_ugm3 = '\u03bcg $\mathregular{m^{-3}}$'

min_o3 = 0.0
max_o3 = 100.0
int_o3 = 5.0

min_pm25 = 0.0
max_pm25 = 150.0
int_pm25 = 5.0

## Ju-Hye built a master site list of all the stations that report during our period of interest
sites_pm25	= pathlib.Path('/','glade','campaign','ral','nsap','JTTI','AIRNOW_MASTER_COORDI_PM25.csv')
sites_o3		= pathlib.Path('/','glade','campaign','ral','nsap','JTTI','AIRNOW_MASTER_COORDI_O3.csv')

df_pm25 = pd.read_csv(sites_pm25, header=0)
pm25_site_ids = np.asarray(df_pm25['Site_ID'])
pm25_site_lon = np.asarray(df_pm25[' Site_Longitude'])
pm25_site_lat = np.asarray(df_pm25[' Site_latitude'])
df_o3 = pd.read_csv(sites_o3, header=0)
o3_site_ids = np.asarray(df_o3['Site_ID'])
o3_site_lon = np.asarray(df_o3[' Site_Longitude'])
o3_site_lat = np.asarray(df_o3[' Site_latitude'])

## The station longitudes are in range (-180.0, 180.0]
## The CAMS longitudes are in range [0.0, 360.0)
## This mismatch will mess up RectBivariateSpline interpolation & extracted values
## So one of these needs to be converted. For now go with converting station longitudes to match CAMS range.
#pm25_site_lon = np.where(pm25_site_lon < 0.0, pm25_site_lon + 360.0, pm25_site_lon)      
#o3_site_lon   = np.where(o3_site_lon < 0.0, o3_site_lon + 360.0, o3_site_lon)

n_sites_pm25 = len(pm25_site_ids)
n_sites_o3 = len(o3_site_ids)

## Build array of dates
beg_dt = pd.to_datetime(beg_date, format=fmt_date_hh)
end_dt = pd.to_datetime(end_date, format=fmt_date_hh)

## Get Cartopy plotting features if any plots are desired
if plot_o3 or plot_pm25:
	borders, states, oceans, lakes, rivers, land = map_functions.get_cartopy_features()
	suptitle_y = 0.77

## Loop through dates looking for files
this_dt = beg_dt
while this_dt <= end_dt:
	this_yr = this_dt.strftime('%Y')
	this_mo = this_dt.strftime('%m')
	this_dy = this_dt.strftime('%d')
	this_hr = this_dt.strftime('%H')
	this_yyyymmdd = this_dt.strftime(fmt_yyyymmdd)
	this_date_str = this_dt.strftime(fmt_date)
	this_date_hh_str = this_dt.strftime(fmt_date_hh)
	this_date_plot = this_dt.strftime(fmt_ddmmyyyyhhmm)+' UTC'

	if cams_fcst:
		if int(this_hr) < 12:
			this_cycle = this_yyyymmdd+'_00'
		else:
			this_cycle = this_yyyymmdd+'_12'

		this_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst',this_yr,this_mo,this_cycle)
		this_cycle_dt = pd.to_datetime(this_cycle, format=fmt_yyyymmdd_hh)
		this_cycle_plot = this_cycle_dt.strftime(fmt_ddmmyyyyhhmm)+' UTC'
	elif cams_eac4:
		this_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4',this_yr,this_mo)

	## If file exists, open it
	pm25_fname = this_dir.joinpath('pm2.5_'+this_date_hh_str+'00.nc')
	o3_fname   = this_dir.joinpath('o3sfc_'+this_date_hh_str+'00.nc')
	get_cams_lat_lon = True

	if o3_fname.is_file():
		print('Reading file '+str(o3_fname))
		ds_o3 = xr.open_dataset(o3_fname)
		if get_cams_lat_lon:
			## Need to reverse cams latitude array to be strictly increasing, rather than strictly decreasing
			cams_lat = ds_o3.latitude.values[::-1]
			cams_lon = ds_o3.longitude.values
			## The station longitudes are in range (-180.0, 180.0]
			## The CAMS longitudes are in range [0.0, 360.0)
			## This mismatch will mess up RectBivariateSpline interpolation & extracted values
			## So one of these needs to be converted. For now convert CAMS longitudes to station longitudes range.
			cams_lon = np.where(cams_lon > 180.0, cams_lon - 360.0, cams_lon)
			get_cams_lat_lon = False
		## Need to reverse the latitude dimension just as cams_lat was reversed
		cams_o3 = ds_o3.go3[:,::-1,:]
		## Convert O3 mass mixing ratio (kg/kg) to volume mixing ratio (ppb)
		## https://confluence.ecmwf.int/pages/viewpage.action?pageId=153391710
		cams_o3_ppb = cams_o3.values * 1e9 * (28.9644 / 47.9982)

		if plot_o3:
			fname = this_dir.joinpath('map_o3sfc_'+this_date_hh_str+'00.'+plot_type)
			var = cams_o3_ppb[0,:,:]
			var_name = 'Surface-level Ozone Concentration'
			var_unit = 'ppbv'
			var_min = np.min(var)
			var_max = np.max(var)
			if cams_fcst:
				suptitle = 'CAMS Forecast'
			elif cams_eac4:
				suptitle = 'CAMS Reanalysis (EAC4)'
			cbar_lab = var_name+' ['+var_unit+']'
			cart_proj = ccrs.PlateCarree()
			cart_xlim = ([np.min(cams_lon), np.max(cams_lon)])
			cart_ylim = ([np.min(cams_lat), np.max(cams_lat)])
			extend = 'max'
			cmap = mpl.cm.rainbow
			bounds = np.arange(min_o3, max_o3+0.01, int_o3)
			norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)
			title_l = var_name+f'\nMin: {var_min:.2f} '+var_unit+f', Max: {var_max:.2f} '+var_unit
			if cams_fcst:
				title_r = 'Start: '+this_cycle_plot+'\nValid: '+this_date_plot
			elif cams_eac4:
				title_r = 'Valid: '+this_date_plot
			map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
				cams_lon, cams_lat, cmap, bounds, norm, extend,
				borders=borders, states=states, lakes=lakes, water_color='none',
				title_l=title_l, title_r=title_r, suptitle_y=suptitle_y)

		## Interpolate O3 data using bivariate spline. Already on rectangular lat-lon grid, so no transformation required.
		## There cannot be any nans in the input 2D array of values, so use array where nan was replaced with _FillValue.
		## This will necessarily impact values near the boundary of real data, so be careful at the edges.
		## To keep extracted values of the spline positive, set the degree of the spline to 1 (true bilinear interp??)
		cams_o3_ppb_spline = sp.interpolate.RectBivariateSpline(cams_lat, cams_lon, cams_o3_ppb[0,:,:], kx=1, ky=1)
		cams_o3_ppb_airnow = cams_o3_ppb_spline.ev(o3_site_lat, o3_site_lon)

		## Create pandas dataframe
		df_o3_cols = ['site_id', 'site_lat', 'site_lon', 'cams_o3_ppb']
		df_cams_o3_ppb_airnow = pd.DataFrame(data=[o3_site_ids, o3_site_lat, o3_site_lon, cams_o3_ppb_airnow])
		df_cams_o3_ppb_airnow = df_cams_o3_ppb_airnow.T
		df_cams_o3_ppb_airnow.columns = df_o3_cols
		df_cams_o3_ppb_airnow.site_id = df_cams_o3_ppb_airnow.site_id.astype(int)

		## Write out CSV file for this valid time
		out_file = this_dir.joinpath('cams_o3sfc_at_airnow_'+this_date_hh_str+'00.csv')
		print('Writing interpolated data to '+str(out_file))
		df_cams_o3_ppb_airnow.to_csv(out_file, index=False)
		
	if pm25_fname.is_file():
		print('Reading file '+str(pm25_fname))
		ds_pm25 = xr.open_dataset(pm25_fname)
		if get_cams_lat_lon:
			## Need to reverse cams latitude array to be strictly increasing, rather than strictly decreasing
			cams_lat = ds_pm25.latitude.values[::-1]
			cams_lon = ds_pm25.longitude.values
			## The station longitudes are in range (-180.0, 180.0]
			## The CAMS longitudes are in range [0.0, 360.0)
			## This mismatch will mess up RectBivariateSpline interpolation & extracted values
			## So one of these needs to be converted. For now convert CAMS longitudes to station longitudes range.
			cams_lon = np.where(cams_lon > 180.0, cams_lon - 360.0, cams_lon)
			get_cams_lat_lon = False
		cams_pm25 = ds_pm25.pm2p5
		## Convert PM2.5 concentration from kg/m3 to ug/m3
		cams_pm25_ug = cams_pm25.values[:,::-1,:] * 1e9

		if plot_pm25:
			fname = this_dir.joinpath('map_pm2.5_'+this_date_hh_str+'00.'+plot_type)
			var = cams_pm25_ug[0,:,:]
			var_name = 'Surface-level '+mpl_pm25+' Concentration'
			var_unit = mpl_ugm3
			var_min = np.min(var)
			var_max = np.max(var)
			if cams_fcst:
				suptitle = 'CAMS Forecast'
			elif cams_eac4:
				suptitle = 'CAMS Reanalysis (EAC4)'
			cbar_lab = var_name+' ['+var_unit+']'
			cart_proj = ccrs.PlateCarree()
			cart_xlim = ([np.min(cams_lon), np.max(cams_lon)])
			cart_ylim = ([np.min(cams_lat), np.max(cams_lat)])
			extend = 'max'
			cmap = mpl.cm.rainbow
			bounds = np.arange(min_pm25, max_pm25+0.01, int_pm25)
			norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)
			title_l = var_name+f'\nMin: {var_min:.2f} '+var_unit+f', Max: {var_max:.2f} '+var_unit
			if cams_fcst:
				title_r = 'Start: '+this_cycle_plot+'\nValid: '+this_date_plot
			elif cams_eac4:
				title_r = 'Valid: '+this_date_plot
			map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
				cams_lon, cams_lat, cmap, bounds, norm, extend,
				borders=borders, states=states, lakes=lakes, water_color='none',
				title_l=title_l, title_r=title_r, suptitle_y=suptitle_y)

		## Interpolate O3 data using bivariate spline. Already on rectangular lat-lon grid, so no transformation required.
		## There cannot be any nans in the input 2D array of values, so use array where nan was replaced with _FillValue.
		## This will necessarily impact values near the boundary of real data, so be careful at the edges.
		## To keep extracted values of the spline positive, set the degree of the spline to 1 (true bilinear interp??)
		cams_pm25_ug_spline = sp.interpolate.RectBivariateSpline(cams_lat, cams_lon, cams_pm25_ug[0,:,:], kx=1, ky=1)
		cams_pm25_ug_airnow = cams_pm25_ug_spline.ev(pm25_site_lat, pm25_site_lon)

		## Create pandas dataframe
		df_pm25_cols = ['site_id', 'site_lat', 'site_lon', 'cams_pm25_ug']
		df_cams_pm25_ug_airnow = pd.DataFrame(data=[pm25_site_ids, pm25_site_lat, pm25_site_lon, cams_pm25_ug_airnow])
		df_cams_pm25_ug_airnow = df_cams_pm25_ug_airnow.T
		df_cams_pm25_ug_airnow.columns = df_pm25_cols
		df_cams_pm25_ug_airnow.site_id = df_cams_pm25_ug_airnow.site_id.astype(int)

		## Write out CSV file for this valid time
		out_file = this_dir.joinpath('cams_pm2.5_at_airnow_'+this_date_hh_str+'00.csv')
		print('Writing interpolated data to '+str(out_file))
		df_cams_pm25_ug_airnow.to_csv(out_file, index=False)

	## Increment loop
	this_dt = this_dt + dt.timedelta(hours=1)
