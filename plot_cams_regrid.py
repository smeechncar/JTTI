'''
plot_cams_regrid.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 7 Apr 2022

This script plots CAMS or other model data that has been regridded to the CMAQ grid.
'''

import os
import sys
import pathlib
import datetime as dt
import pandas as pd
import numpy as np
import xarray as xr
import wrf
import netCDF4
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
from functions import map_functions

def main():
	beg_date = '2020-08-13_00'
	end_date = '2020-08-15_23'

	cams_fc = True
	cams_ra = False

	plot_o3   = True
	plot_pm25 = True
	plot_type = 'png'

	if plot_o3 or plot_pm25:
		set_cart_params = True

	fmt_date = '%Y-%m-%d'
	fmt_date_hh = '%Y-%m-%d_%H'
	fmt_date_hhmm = '%Y-%m-%d_%H%M'
	fmt_yyyymmdd = '%Y%m%d'
	fmt_yyyymmdd_hh = '%Y%m%d_%H'
	fmt_yyyymmdd_hhmm = '%Y%m%d_%H%M'
	fmt_ddmmyyyyhhmm = '%d %b %Y/%H%M'

	mpl_pm25 = '$\mathregular{PM_{2.5}}$'
	mpl_ugm3 = '\u03bcg $\mathregular{m^{-3}}$'

	min_o3 = 0.0
	max_o3 = 90.1
	int_o3 = 5.0

	min_pm25 = 0.0
	max_pm25 = 50.1
	int_pm25 = 5.0
    cont_pm25 = np.arange(min_pm25, max_pm25, int_pm25)
    cont_pm25 = np.append(cont_pm25, np.arange(60.0, 100.1, 10.0))
    cont_pm25 = np.append(cont_pm25, np.arange(125.0, 200.1, 25.0))

	data_dir_parent_cams_fc = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')
	data_dir_parent_cams_ra = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4')

	## Build array of dates
	beg_dt = pd.to_datetime(beg_date, format=fmt_date_hh)
	end_dt = pd.to_datetime(end_date, format=fmt_date_hh)
	all_dt = pd.date_range(start=beg_dt, end=end_dt, freq='1H')
	n_dates = len(all_dt)

	## Loop through dates looking for files
	for tt in range(n_dates):
		this_dt = all_dt[tt]
		this_yr = this_dt.strftime('%Y')
		this_mo = this_dt.strftime('%m')
		this_dy = this_dt.strftime('%d')
		this_hr = this_dt.strftime('%H')
		this_yyyymmdd = this_dt.strftime(fmt_yyyymmdd)
		this_yyyymmdd_hhmm = this_dt.strftime(fmt_yyyymmdd_hhmm)
		this_date_str = this_dt.strftime(fmt_date)
		this_date_hh_str = this_dt.strftime(fmt_date_hh)
		this_date_hhmm_str = this_dt.strftime(fmt_date_hhmm)
		this_date_plot = this_dt.strftime(fmt_ddmmyyyyhhmm)+' UTC'

		if cams_fc:
			if int(this_hr) < 12:
				this_cycle = this_yyyymmdd+'_00'
			else:
				this_cycle = this_yyyymmdd+'_12'

			this_dir = data_dir_parent_cams_fc.joinpath(this_yr, this_mo, this_cycle)
			this_cycle_dt = pd.to_datetime(this_cycle, format=fmt_yyyymmdd_hh)
			this_cycle_plot = this_cycle_dt.strftime(fmt_ddmmyyyyhhmm)+' UTC'
			fname = this_dir.joinpath('cams_o3_pm25_regrid_cmaq_'+this_yyyymmdd_hhmm+'.nc')
		elif cams_ra:
			this_dir = data_dir_parent_cams_ra.joinpath(this_yr, this_mo)
			fname = this_dir.joinpath('cams_o3_pm25_regrid_cmaq_'+this_yyyymmdd_hhmm+'.nc')

		## If file exists, open it
		if fname.is_file():
			print('Reading file '+str(fname))
			ds_xr = xr.open_dataset(fname)
#			ds_nc = netCDF4.Dataset(fname, mode='r')
			if cams_fc or cams_ra:
				latitude = ds_xr.latitude
				longitude = ds_xr.longitude
				latitude = latitude.rename({'latitude':'XLAT', 'longitude':'XLONG'})
				longitude = longitude.rename({'latitude':'XLAT', 'longitude':'XLONG'})
				if set_cart_params:
					borders, states, oceans, lakes, rivers, land = map_functions.get_cartopy_features()
					data_crs = ccrs.PlateCarree()
					suptitle_y = 0.92

					## Get the WRF-style projection attributes
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

					## Now get the cartopy projection x and y limits
					latitude = latitude.assign_attrs(projection=wrf_proj)
					longitude = longitude.assign_attrs(projection=wrf_proj)
					cart_xlim = wrf.cartopy_xlim(var=latitude)
					cart_ylim = wrf.cartopy_ylim(var=latitude)

					set_cart_params = False    # don't go through this loop again

				o3 = ds_xr.cams_o3[0,:,:]
				o3 = o3.rename({'latitude':'XLAT', 'longitude':'XLONG'})
				o3 = o3.assign_attrs(projection=wrf_proj)
				pm25 = ds_xr.cams_pm25[0,:,:]
				pm25 = pm25.rename({'latitude':'XLAT', 'longitude':'XLONG'})
				pm25 = pm25.assign_attrs(projection=wrf_proj)
#				latname, lonname, _ = wrf.util._find_coord_names(o3.coords)

				if cams_fc:
					suptitle = 'CAMS Forecast, Regridded to CMAQ Grid'
					title_r  = 'Start: '+this_cycle_plot+'\nValid: '+this_date_plot
				elif cams_ra:
					suptitle = 'CAMS Reanalysis (EAC4), Regridded to CMAQ Grid'
					title_r  = 'Valid: '+this_date_plot

				if plot_o3:
					fname = this_dir.joinpath('map_regrid_o3sfc_'+this_date_hhmm_str+'.'+plot_type)
					var = o3
					var_name = 'Surface-level Ozone Concentration'
					var_unit = 'ppbv'
					var_min = np.nanmin(var)
					var_max = np.nanmax(var)
					cbar_lab = var_name+' ['+var_unit+']'
					extend = 'max'
					cmap = mpl.cm.rainbow
					bounds = np.arange(min_o3, max_o3, int_o3)
					norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)
					title_l = var_name+f'\nMin: {var_min:.2f} '+var_unit+f', Max: {var_max:.2f} '+var_unit
					map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
						longitude, latitude, cmap, bounds, norm, extend,
						borders=borders, states=states, lakes=lakes, water_color='none',
						title_l=title_l, title_r=title_r, suptitle_y=suptitle_y)

				if plot_pm25:
					fname = this_dir.joinpath('map_regrid_pm2.5_'+this_date_hhmm_str+'.'+plot_type)
					var = pm25
					var_name = 'Surface-level '+mpl_pm25+' Concentration'
					var_unit = mpl_ugm3
					var_min = np.nanmin(var)
					var_max = np.nanmax(var)
					cbar_lab = var_name+' ['+var_unit+']'
					extend = 'max'
					cmap = mpl.cm.rainbow
#					bounds = np.arange(min_pm25, max_pm25+0.0001, int_pm25)
                    bounds = cont_pm25
					norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend=extend)
					title_l = var_name+f'\nMin: {var_min:.2f} '+var_unit+f', Max: {var_max:.2f} '+var_unit
					map_functions.mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim,
						longitude, latitude, cmap, bounds, norm, extend,
						borders=borders, states=states, lakes=lakes, water_color='none',
						title_l=title_l, title_r=title_r, suptitle_y=suptitle_y)
				

if __name__ == '__main__':
	now_time_beg = dt.datetime.utcnow()
	main()
	now_time_end = dt.datetime.utcnow()
	run_time_tot = now_time_end - now_time_beg
	now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
	now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
	print('\nScript completed successfully.')
	print('   Beg time: '+now_time_beg_str)
	print('   End time: '+now_time_end_str)
	print('   Run time: '+str(run_time_tot)+'\n')
