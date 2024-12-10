'''
plot_use_val_map_conus_calif.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 20 Dec 2022

This script reads in PM2.5/O3 use/val files, and creates maps indicating which sites are use/val.
On the CONUS map, a box showing the outline of the California subdomain is optionally shown.
Separate maps of just the California subdomain are also optionally created.
'''

import sys
import pathlib
import wrf
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import shapely
from functions import map_functions

def plot_map(fname, lon_use, lat_use, lon_val, lat_val, suptitle, cart_proj, cart_xlim, cart_ylim,
	borders=None, states=None, oceans=None, lakes=None, water_color='lightblue',
	lat_labels=None, lon_labels=None, polygon_verts=None, fontsize=12,
	suptitle_y=0.98, title_l=None, title_c=None, title_r=None, marker_size=81, legend_loc='lower left',
	):

	mpl.rcParams['figure.figsize'] = (10, 8)
	mpl.rcParams['grid.color'] = 'gray'
	mpl.rcParams['grid.linestyle'] = ':'
	mpl.rcParams['font.size'] = fontsize+2
	mpl.rcParams['figure.titlesize'] = fontsize+2
	mpl.rcParams['savefig.bbox'] = 'tight'
	mpl.rcParams['legend.framealpha'] = 0.95
	ll_size = fontsize-2
	data_crs = ccrs.PlateCarree()

	color_use = 'tab:blue'
	color_val = 'tab:orange'
	marker_use = 'v'
	marker_val = 'o'

	print('-- Plotting: '+str(fname))
	fig = plt.figure()
	ax = plt.subplot(projection=cart_proj)

	## If cart_xlim and cart_ylim tuples are not provided, then set plot limits from lat/lon data directly
	if cart_xlim is not None and cart_ylim is not None:
		ax.set_xlim(cart_xlim)
		ax.set_ylim(cart_ylim)
	else:
		ax.set_extent([np.min(lons), np.max(lons), np.min(lats), np.max(lats)], crs=cart_proj)

	## Sometimes longitude labels show up on y-axis, and latitude labels on x-axis in older versions of cartopy
	## Print lat/lon labels only for a specified set (determined by trial & error) to avoid this problem for now
	gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
	gl.rotate_labels = False
	## If specific lat/lon labels are not specified, then just label the default gridlines
	if lon_labels is None:
		gl.top_labels = True
		gl.bottom_labels = True
	else:
		gl.xlocator = mticker.FixedLocator(lon_labels)
	if lat_labels is None:
		gl.left_labels = True
		gl.right_labels = True
	else:
		gl.ylocator = mticker.FixedLocator(lat_labels)
	gl.xlabel_style = {'size':ll_size}
	gl.ylabel_style = {'size':ll_size}

	if borders != None:
#		ax.add_feature(borders, linestyle='-', zorder=3)
		ax.add_feature(borders, linestyle='-')
	if states != None:
#		ax.add_feature(states, linewidth=0.25, edgecolor='black', zorder=4)
		ax.add_feature(states, linewidth=0.25, edgecolor='black')
	if oceans != None:
#		ax.add_feature(oceans, facecolor=water_color, zorder=2)
		ax.add_feature(oceans, facecolor=water_color)
	if lakes != None:
#		ax.add_feature(lakes, facecolor=water_color, linewidth=0.25, edgecolor='black', zorder=5)
		ax.add_feature(lakes, facecolor=water_color, linewidth=0.25, edgecolor='black')
#	ax.coastlines(zorder=6)
	ax.coastlines()

	## Add markers to the plot
	ax.scatter(lon_val, lat_val, marker=marker_val, color=color_val, s=marker_size,
		edgecolors='black', transform=data_crs, zorder=10)
	ax.scatter(lon_use, lat_use, marker=marker_use, color=color_use, s=marker_size,
		edgecolors='black', transform=data_crs, zorder=10)

	## Add the legend
	ax.legend(['Validation', 'Training'], loc=legend_loc)

	## Optional: Draw a polygon
	if polygon_verts is not None:
		n_vertices = polygon_verts.shape[0]
		## Draw a line between each vertex of the polygon
		for vv in range(n_vertices):
			lon_beg = polygon_verts[vv][0]
			lat_beg = polygon_verts[vv][1]
			if vv < n_vertices-1:
				lon_end = polygon_verts[vv+1][0]
				lat_end = polygon_verts[vv+1][1]
			else: # loop back around to the first vertex
				lon_end = polygon_verts[0][0]
				lat_end = polygon_verts[0][1]
			## Plot each side of the polygon (could optionally add a marker at each vertex, but not doing that here)
			ax.plot([lon_beg,lon_end], [lat_beg,lat_end], color='tab:green', linewidth=3, transform=data_crs, zorder=20)

	## Add the plot overall title
	plt.suptitle(suptitle, y=suptitle_y)

	## Optional: Add titles to the subplot
	## In this case, one title will be on the left side, one on the right side
	if title_l != None:
		ax.set_title(title_l, fontsize=fontsize, loc='left')
	if title_r != None:
		ax.set_title(title_r, fontsize=fontsize, loc='right')
	if title_c != None:
		ax.set_title(title_c, fontsize=fontsize, loc='center')

	## Save and close the figure
	plt.savefig(fname)
	plt.close()


def main():
	obs_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','airnow','split','sites_static')
	out_dir = obs_dir

	date_beg = '20200801_00'
	date_end = '20211231_23'
	plot_type = 'png'

	i_beg = 20
	i_end = 60
	j_beg = 130
	j_end = 170

	fmt_yyyymmdd_hh = '%Y%m%d_%H'
	fmt_yyyymmddhh  = '%Y%m%d%H'
	fmt_yyyymmdd = '%Y%m%d'
	fmt_yyyy = '%Y'
	fmt_mm = '%m'
	fmt_hh = '%H'
	fmt_date = '%Y-%m-%d'
	fmt_date_hh = '%Y-%m-%d_%H'

	pm25_str = '$\mathregular{PM_{2.5}}$'
	o3_str = '$\mathregular{O_{3}}$'

	en_dash = u'\u2013'
	em_dash = u'\u2014'

	## Create datetime array
	dt_beg = dt.datetime.strptime(date_beg, fmt_yyyymmdd_hh)
	dt_end = dt.datetime.strptime(date_end, fmt_yyyymmdd_hh)
	dt_array = pd.date_range(start=dt_beg, end=dt_end, freq='1H')
	n_times = len(dt_array)

	## Open single files that contain all times
	use_file = obs_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_use.nc')
	val_file = obs_dir.joinpath('airnow_pm2.5_o3_valid_uptime_'+date_beg+'-'+date_end+'_val.nc')

	print('Reading '+str(use_file))
	ds_use = xr.open_dataset(use_file)
	pm25_sid_use = ds_use.pm25_sid_use
	pm25_lon_use = ds_use.pm25_lon_use
	pm25_lat_use = ds_use.pm25_lat_use
	o3_sid_use = ds_use.o3_sid_use
	o3_lon_use = ds_use.o3_lon_use
	o3_lat_use = ds_use.o3_lat_use
	n_use_pm25 = len(pm25_sid_use)
	n_use_o3   = len(o3_sid_use)

	print('Reading '+str(val_file))
	ds_val = xr.open_dataset(val_file)
	pm25_sid_val = ds_val.pm25_sid_val
	pm25_lon_val = ds_val.pm25_lon_val
	pm25_lat_val = ds_val.pm25_lat_val
	o3_sid_val = ds_val.o3_sid_val
	o3_lon_val = ds_val.o3_lon_val
	o3_lat_val = ds_val.o3_lat_val
	n_val_pm25 = len(pm25_sid_val)
	n_val_o3   = len(o3_sid_val)

	## ================
   ## PLOTTING SECTION
   ## ================

	## Now open the CMAQ sample/coordinate file to get map projection parameters
	cmaq_fname = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','aqm.t12z.grdcro2d.ncf')
	print('Reading CMAQ coordinate data from '+str(cmaq_fname))
	cmaq_ds = xr.open_dataset(cmaq_fname)
	cmaq_lon = cmaq_ds.LON[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})
	cmaq_lat = cmaq_ds.LAT[0,0,:,:].rename({'ROW':'latitude', 'COL':'longitude'})
	cmaq_lon = cmaq_lon.assign_coords(coords={'XLAT':cmaq_lat, 'XLONG':cmaq_lon})
	cmaq_lat = cmaq_lat.assign_coords(coords={'XLAT':cmaq_lat, 'XLONG':cmaq_lon})
	cmaq_lon.attrs['long_name'] = 'longitude'
	cmaq_lat.attrs['long_name'] = 'latitude'
	cmaq_lon.attrs['units'] = 'degrees_east'
	cmaq_lat.attrs['units'] = 'degrees_north'
	cmaq_lon_sub = cmaq_ds.LON[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
	cmaq_lat_sub = cmaq_ds.LAT[0,0,j_beg:j_end,i_beg:i_end].rename({'ROW':'latitude', 'COL':'longitude'})
	cmaq_lon_sub = cmaq_lon_sub.assign_coords(coords={'XLAT':cmaq_lat_sub, 'XLONG':cmaq_lon_sub})
	cmaq_lat_sub = cmaq_lat_sub.assign_coords(coords={'XLAT':cmaq_lat_sub, 'XLONG':cmaq_lon_sub})
	cmaq_lon_sub.attrs['long_name'] = 'longitude'
	cmaq_lat_sub.attrs['long_name'] = 'latitude'
	cmaq_lon_sub.attrs['units'] = 'degrees_east'
	cmaq_lat_sub.attrs['units'] = 'degrees_north'
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

	## Polygon vertices
	lats_vect = np.array([cmaq_lat[j_beg,i_beg],cmaq_lat[j_end,i_beg],cmaq_lat[j_end,i_end],cmaq_lat[j_beg,i_end]])
	lons_vect = np.array([cmaq_lon[j_beg,i_beg],cmaq_lon[j_end,i_beg],cmaq_lon[j_end,i_end],cmaq_lon[j_beg,i_end]])
	ll_verts = np.column_stack((lons_vect,lats_vect))
	polygon = shapely.geometry.polygon.Polygon(ll_verts)

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

	cmaq_lat_sub = cmaq_lat_sub.assign_attrs(projection=wrf_proj)
	cmaq_lon_sub = cmaq_lon_sub.assign_attrs(projection=wrf_proj)
	cart_xlim_sub = wrf.cartopy_xlim(var=cmaq_lat_sub)
	cart_ylim_sub = wrf.cartopy_ylim(var=cmaq_lat_sub)

	## Get Cartopy features
	borders, states, oceans, lakes, rivers, land = map_functions.get_cartopy_features()

	## Set some other plot attributes
	suptitle_pm25_use_val = 'AirNow '+pm25_str+' Observation Sites'
	suptitle_o3_use_val   = 'AirNow '+o3_str+' Observation Sites'

	fontsize = 12
	marker_size = 64
	suptitle_y = 0.85

	## First, make maps of CONUS with use stations one color and val stations another color
	lon_labels = [-130, -120, -110, -110, -100, -90, -80, -70, -60]
	lat_labels = [25, 30, 35, 40, 45, 50]

	fname = out_dir.joinpath('map_conus_pm2.5_use+val_'+date_beg+'-'+date_end+'.'+plot_type)
	plot_map(fname, pm25_lon_use, pm25_lat_use, pm25_lon_val, pm25_lat_val,
		suptitle_pm25_use_val, cart_proj, cart_xlim, cart_ylim,
		borders=borders, states=states, oceans=oceans, lakes=lakes,
		lat_labels=lat_labels, lon_labels=lon_labels, fontsize=fontsize, polygon_verts=None,
		suptitle_y=suptitle_y, marker_size=marker_size)

	fname = out_dir.joinpath('map_conus_o3_use+val_'+date_beg+'-'+date_end+'.'+plot_type)
	plot_map(fname, o3_lon_use, o3_lat_use, o3_lon_val, o3_lat_val,
		suptitle_o3_use_val, cart_proj, cart_xlim, cart_ylim,
		borders=borders, states=states, oceans=oceans, lakes=lakes,
		lat_labels=lat_labels, lon_labels=lon_labels, fontsize=fontsize, polygon_verts=None,
		suptitle_y=suptitle_y, marker_size=marker_size)

	## Second, add a polygon outlining the California subdomain
	fname = out_dir.joinpath('map_conus+calif_pm2.5_use+val_'+date_beg+'-'+date_end+'.'+plot_type)
	plot_map(fname, pm25_lon_use, pm25_lat_use, pm25_lon_val, pm25_lat_val,
		suptitle_pm25_use_val, cart_proj, cart_xlim, cart_ylim,
		borders=borders, states=states, oceans=oceans, lakes=lakes,
		lat_labels=lat_labels, lon_labels=lon_labels, fontsize=fontsize, polygon_verts=ll_verts,
		suptitle_y=suptitle_y, marker_size=marker_size)

	fname = out_dir.joinpath('map_conus+calif_o3_use+val_'+date_beg+'-'+date_end+'.'+plot_type)
	plot_map(fname, o3_lon_use, o3_lat_use, o3_lon_val, o3_lat_val,
		suptitle_o3_use_val, cart_proj, cart_xlim, cart_ylim,
		borders=borders, states=states, oceans=oceans, lakes=lakes,
		lat_labels=lat_labels, lon_labels=lon_labels, fontsize=fontsize, polygon_verts=ll_verts,
		suptitle_y=suptitle_y, marker_size=marker_size)

	## Third, make maps of just the California subdomain
	suptitle_y = 0.95
	lat_labels = [36, 37, 38, 39, 40, 41]
	lon_labels = [-124, -123, -122, -121, -120, -119, -118]
	fname = out_dir.joinpath('map_calif_pm2.5_use+val_'+date_beg+'-'+date_end+'.'+plot_type)
	plot_map(fname, pm25_lon_use, pm25_lat_use, pm25_lon_val, pm25_lat_val,
		suptitle_pm25_use_val, cart_proj, cart_xlim_sub, cart_ylim_sub,
		borders=borders, states=states, oceans=oceans, lakes=lakes,
		lat_labels=lat_labels, lon_labels=lon_labels, fontsize=fontsize, polygon_verts=None,
		suptitle_y=suptitle_y, marker_size=marker_size, legend_loc='upper right')

	fname = out_dir.joinpath('map_calif_o3_use+val_'+date_beg+'-'+date_end+'.'+plot_type)
	plot_map(fname, o3_lon_use, o3_lat_use, o3_lon_val, o3_lat_val,
		suptitle_o3_use_val, cart_proj, cart_xlim_sub, cart_ylim_sub,
		borders=borders, states=states, oceans=oceans, lakes=lakes,
		lat_labels=lat_labels, lon_labels=lon_labels, fontsize=fontsize, polygon_verts=None,
		suptitle_y=suptitle_y, marker_size=marker_size, legend_loc='upper right')

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
