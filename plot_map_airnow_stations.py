'''
plot_map_airnow_stations.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 5 Jan 2023

This script reads in a CMAQ file for map/projection info and AirNow obs file(s)
and plots the PM2.5 and O3 stations on a map.
'''

import sys
import pathlib
import argparse
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import wrf
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
from functions import map_functions

def main():

	plot_type = 'png'
	plot_subdomain = False

	if plot_subdomain:
		i_beg, i_end = 20, 60
		j_beg, j_end = 130, 170
		lat_labels = [36, 37, 38, 39, 40, 41]
		lon_labels = [-124, -123, -122, -121, -120, -119, -118]
		suptitle_y = 1.01
	else:
		i_beg, i_end = 0, -1
		j_beg, j_end = 0, -1
		lat_labels = [25, 30, 35, 40, 45, 50]
		lon_labels = [-130, -120, -110, -110, -100, -90, -80, -70, -60]
		suptitle_y = 0.85

	mpl_o3   = '$\mathregular{O_3}$'
	mpl_pm25 = '$\mathregular{PM_{2.5}}$'
	mpl_ugm3 = '\u03bcg $\mathregular{m^{-3}}$'

	out_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')

	obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')
	date_str = '20200801_to_20220101'
	obs_file = obs_dir.joinpath('cams_fcst_o3_pm25_'+date_str+'.nc')
	print('Reading AirNow station metadata from '+str(obs_file))
	obs_ds = xr.open_dataset(obs_file)
	obs_lat_pm25 = obs_ds.pm25_lat
	obs_lon_pm25 = obs_ds.pm25_lon
	obs_lat_o3 = obs_ds.o3_lat
	obs_lon_o3 = obs_ds.o3_lon

	## Open a CMAQ sample/coordinate file to get map projection parameters
	cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
	print('Reading CMAQ coordinate data from '+str(cmaq_fname))
	cmaq_ds = xr.open_dataset(cmaq_fname)
	if plot_subdomain:
		cmaq_lon = cmaq_ds.LON[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
		cmaq_lat = cmaq_ds.LAT[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
	else:
		cmaq_lon = cmaq_ds.LON[0,0,j_beg:,i_beg:].rename({'ROW':'latitude', 'COL':'longitude'})
		cmaq_lat = cmaq_ds.LAT[0,0,j_beg:,i_beg:].rename({'ROW':'latitude', 'COL':'longitude'})
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

	## PLOTTING SECTION
	lat_labels = [25, 30, 35, 40, 45, 50]
	lon_labels = [-130, -120, -110, -110, -100, -90, -80, -70, -60]
	suptitle = 'AirNow Stations'

	pm25_mark = 'v'
	o3_mark = 'o'
	pm25_size = 81
	o3_size = 81
	pm25_col = 'tab:orange'
	o3_col = 'tab:green'

	plot_types = ['png', 'eps']
	for plot_type in plot_types:
		## Make the PM2.5 map
		fname = 'map_all_pm2.5_stations_'+date_str+'.'+plot_type
		plot_file = out_dir.joinpath(fname)
		lg_text = [mpl_pm25]
		map_functions.mpl_map_stations(plot_file, suptitle, cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat,
			borders=borders, states=states, oceans=oceans, lakes=lakes,
			lat_labels=lat_labels, lon_labels=lon_labels, suptitle_y=suptitle_y,
			mark1_lon=obs_lon_pm25, mark1_lat=obs_lat_pm25, mark1=pm25_mark, mark1_size=pm25_size,mark1_color=pm25_col,
			lg_text=lg_text, )

		## Make the O3 map
		fname = 'map_all_o3_stations_'+date_str+'.'+plot_type
		plot_file = out_dir.joinpath(fname)
		lg_text = [mpl_o3]
		map_functions.mpl_map_stations(plot_file, suptitle, cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat,
			borders=borders, states=states, oceans=oceans, lakes=lakes,
			lat_labels=lat_labels, lon_labels=lon_labels, suptitle_y=suptitle_y,
			mark1_lon=obs_lon_o3, mark1_lat=obs_lat_o3, mark1=o3_mark, mark1_size=o3_size, mark1_color=o3_col,
			lg_text=lg_text, )

		## Make the combined map
		fname = 'map_all_pm2.5+o3_stations_'+date_str+'.'+plot_type
		plot_file = out_dir.joinpath(fname)
		lg_text = [mpl_pm25, mpl_o3]
		map_functions.mpl_map_stations(plot_file, suptitle, cart_proj, cart_xlim, cart_ylim, cmaq_lon, cmaq_lat,
			borders=borders, states=states, oceans=oceans, lakes=lakes,
			lat_labels=lat_labels, lon_labels=lon_labels, suptitle_y=suptitle_y,
			mark1_lon=obs_lon_pm25, mark1_lat=obs_lat_pm25, mark1=pm25_mark, mark1_size=pm25_size,mark1_color=pm25_col,
			mark2_lon=obs_lon_o3, mark2_lat=obs_lat_o3, mark2=o3_mark, mark2_size=o3_size, mark2_color=o3_col,
			lg_text=lg_text, )

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
