'''
map_functions.py

Created by: Jared A. Lee
            jaredlee@ucar.edu
Created on: 4 Apr 2019

This script contains commonly-used functions that can be called by other Python scripts.
The primary function is to set some basic map resources for contour plots.

calc_bearing
    -- Function to calculate the bearing angle between two points (lon1,lat1) and (lon2,lat2).

get_cartopy_features
   -- Function to get some commonly used Cartopy features (borders, states, oceans, lakes, rivers, land)
        for use in Matplotlib/Cartopy plots.

truncate_cmap
    -- Function to truncate a matplotlib colormap. Particularly useful when plotting terrain height.

mpl_cross_plot
   -- Procedure to make a WRF vertical cross-section plot with filled contours using matplotlib.
        Optionally plots a small inset map in an upper corner showing the cross-section path.
        Optionally shows a contour-filled field in the inset map.
        Optionally overlays the cross-section plot with wind barbs.

mpl_map_plot
    -- Procedure to make a map plot with filled contours using matplotlib and Cartopy.
        Optionally overlays the map with wind barbs, markers, text labels, a polygon, or a cross-section path.

mpl_map_stations
    -- Procedure to make a map plot with station markers, using matplotlib and Cartopy.

make_skewt
    -- Procedure to make a skew-T plot using MetPy. Optionally overlay a second sounding for comparison.

make_ngl_map_plot_poly
   -- Procedure to make an empty map plot with polymarkers and optional text labels.

make_ngl_contour_map_plot
   -- Procedure to make a plot using Ngl.contour_map.

make_ngl_contour_map_plot_vec
   -- Procedure to make a plot using Ngl.contour_map, overlaid with wind vectors using Ngl.contour_map & Ngl.overlay

make_ngl_contour_map_plot_poly
   -- Procedure to make a plot using Ngl.contour_map, with optional polymarkers & text added.

make_ngl_contour_map_plot_poly_vec
   -- Procedure to make a plot using Ngl.contour_map, overlaid with wind vectors using
        Ngl.contour_map & Ngl.overlay, with optional polymarkers & text added.

ngl_strings
    -- Adds annotations (left, right, or center) above a PyNGL plot.
        Corresponds to NCL's gsnLeftString, gsnCenterString, gsnRightString.

set_ngl_map_res
    -- Sets numerous basic resources for a contour map plot.

set_ngl_vcres
    -- Sets several basic resources for adding wind barbs on top of a contour map plot.

set_ngl_polyres
    -- Sets several basic resources for adding polymarker(s) at location(s) of interest on a map plot.

set_ngl_txres
    -- Sets a basic resource for adding text at location(s) of interest on a map plot.

class Meteogram functions:
    -- plot_ws_wd
    -- plot_t_td
    -- plot_t_td_rh
    -- plot_rh
    -- plot_p
    -- plot_ghi_dni
    -- plot_ghi
    -- plot_rain
    -- plot_hfx_qfx
    -- plot_th_e
    -- plot_lcl_lfc_pblh
'''

#import Ngl
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import wrf
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.cm import get_cmap
import datetime as dt
import metpy
import metpy.calc as mpcalc
import metpy.plots as mpplots
from metpy.units import units as mpunits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes

def calc_bearing(lon1, lat1, lon2, lat2):
    '''
    -- Function to calculate the bearing angle between two points (lon1,lat1) and (lon2,lat2).
    -- Added by JAL on 22 Jan 2022
    -- Inputs:
        - lon1, lat1: floats, longitude and latitude of point 1
        - lon2, lat2: floats, longitude and latitude of point 2
    -- Ouputs:
        - bearing_geog: float, geographical bearing (0 deg = north, 90 deg = east, etc.)
        - bearing_math: float, mathematical bearing (0 deg = east, 90 deg = north, etc.)
    -- Reference:
        - https://www.igismap.com/formula-to-find-bearing-or-heading-angle-between-two-points-latitude-longitude/
    '''
    DEG2RAD = np.pi / 180.0
    RAD2DEG = 180.0 / np.pi
    x = np.cos(lat2 * DEG2RAD) * np.sin((lon2 - lon1) * DEG2RAD)
    y = (np.cos(lat1 * DEG2RAD) * np.sin(lat2 * DEG2RAD)) - (np.sin(lat1 * DEG2RAD) * np.cos(lat2 * DEG2RAD) * np.cos((lon2-lon1)*DEG2RAD))
    bearing_geog = np.arctan2(x, y) * RAD2DEG
    bearing_math = 90.0 - bearing_geog
    ## Probably not necessary, but keeps it in the range (-180, 180]
    if bearing_math <= -180.0:
        bearing_math += 360.0

    return bearing_geog, bearing_math


def get_cartopy_features():
    '''
   -- Function to get some commonly used Cartopy features (borders, states, oceans, lakes, rivers, land)
        for use in Matplotlib/Cartopy plots.
   -- Added by JAL on 13 Jul 2021
   -- Inputs: None.
   -- Outputs:
      - borders, states, oceans, lakes, rivers, land
    '''

    borders = cartopy.feature.BORDERS
    states  = cartopy.feature.NaturalEarthFeature(category='cultural',scale='10m',facecolor='none',name='admin_1_states_provinces_lakes')
    oceans  = cartopy.feature.OCEAN
    lakes   = cartopy.feature.LAKES
    rivers  = cartopy.feature.RIVERS
    land    = cartopy.feature.LAND

    return borders, states, oceans, lakes, rivers, land

def truncate_cmap(cmap, minval=0.0, maxval=1.0, n=100):
    '''
    -- Function to truncate a matplotlib colormap. Particularly useful when plotting terrain.
    -- Added by JAL on 31 Jan 2022
    -- Required Positional Inputs:
        - cmap: Original colormap
    -- Optional Inputs:
        - minval: Fraction into the original colormap to start the new colormap (default: 0.0)
        - maxval: Fraction into the original colormap to end the new colormap (default: 1.0)
        - n: Number of colors to create in the new colormap (default: 100)
    -- Output:
        - new_cmap: New, truncated colormap
    '''
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def mpl_cross_plot(fname, var, var_filled, levs, cmap, norm, extend, xlab, ylab, ymin, ymax, suptitle, cbar_lab, ter_line,
    title_l=None, title_r=None, title_c=None, suptitle_y=0.95, cbar_loc='right', fontsize=14,
    x=None, y=None, u=None, v=None,
    cart_proj=None, cart_xlim=None, cart_ylim=None, lat_beg=None, lon_beg=None, lat_end=None, lon_end=None,
    ins_var=None, ins_bounds=None, ins_cmap=None, ins_norm=None, wrf_lats=None, wrf_lons=None, ins_loc='upper left',
    borders=None, states=None, lakes=None, water_color='lightblue'):

    '''
    -- Procedure to make a WRF vertical cross-section plot with filled contours using matplotlib.
      Optionally plots a small inset map in an upper corner showing the cross-section path.
      Optionally shows a contour-filled field in the inset map.
      Optionally overlays the cross-section plot with wind barbs.
   -- Added by JAL on 21 Jan 2022 (adapted, improved, & overhauled from earlier functions written in July 2021)
   -- Required Positional Inputs:
      - fname: string or pathlib object
      - var: 2-D variable array to be plotted with filled contours <-- NOTE: Might be able to get rid of this
      - var_filled: 2-D variable array to be plotted with filled contours, including solid fill for below terrain
      - levs: 1-D array of contour levels
      - cmap: Matplotlib colormap
      - norm: Matplotlib colormap norm
      - extend: string for colorbar cap ('max', 'min', 'both')
      - xlab: string for x-axis label
      - ylab: string for y-axis label
      - ymin: numeric, y-axis minimum value
      - ymax: numeric, y-axis maxiumum value
      - suptitle: string for overall plot title (usually one line)
      - cbar_lab: string for colorbar label (e.g., Wind Speed [m/s])
      - ter_line: 1-D array of WRF terrain heights along the cross-section line
    -- Optional Inputs:
        - title_l: string, plot subtitle (1 or 2 lines) that gets placed above the top-left corner of the plot axes
        - title_r: string, plot subtitle (1 or 2 lines) that gets placed above the top-right corner of the plot axes
        - title_c: string, plot subtitle (1 or 2 lines) that gets placed above the center of the plot axes
        - suptitle_y: float, y-axis position of the suptitle <-- NOTE: Might be able to delete this
        - cbar_loc: string, identifier for positioning of the colorbar ('top', 'bottom', 'right', 'left')
        - fontsize: text label font size (default: 14). Other plot text sizes are anchored off this value.
      - x, y: array-like, define the barb locations
      - u, v: array-like, define the barb directions
      - cart_proj: Cartopy projection object
      - cart_xlim, cart_ylim: Cartopy objects, x- & y-axis limits for the domain from which the cross-section comes
      - lat_beg, lon_beg, lat_end, lon_end: floats defining the beginning & ending points of the cross-section path
      - ins_var: 2-D variable array to plot in the inset
      - ins_bounds: 1-D array of the contour bounds for the inset variable
      - ins_cmap: Matplotlib colormap for the inset variable
      - ins_norm: Matplotlib colormap norm for the inset variable
      - wrf_lats, wrf_lons: Latitude & longitude cooridnates from an xarray.DataArray object read in by wrf.getvar
        - ins_loc: string setting the location of the inset map ('upper left' or 'upper right')
        - borders, states, lakes: Cartopy feature objects for the inset map
        - water_color: string defining the water color for the inset map (currently not used)
    '''

    mpl.rcParams['figure.figsize'] = (10, 8)
    mpl.rcParams['figure.titlesize'] = fontsize+2
    mpl.rcParams['grid.color'] = 'gray'
    mpl.rcParams['grid.linestyle'] = ':'
    mpl.rcParams['font.size'] = fontsize
    mpl.rcParams['savefig.bbox'] = 'tight'
    mpl.rcParams['lines.markersize'] = 4

    data_crs = ccrs.PlateCarree()

    print('-- Creating plot '+str(fname))
    fig = plt.figure()
    ax_cross = plt.axes()
    xs = np.arange(0, var.shape[-1], 1)
    ys = wrf.to_np(var.coords['vertical'])
#   contours = ax_cross.contourf(xs, ys, wrf.to_np(var_filled), levels=levs, cmap=cmap, norm=norm, extend=extend)
    contours = ax_cross.contourf(xs, ys, wrf.to_np(var), levels=levs, cmap=cmap, norm=norm, extend=extend)
    cbar = fig.colorbar(contours, ax=ax_cross, label=cbar_lab)
    ht_fill = ax_cross.fill_between(xs, 0, wrf.to_np(ter_line), facecolor='saddlebrown')

    coord_pairs = wrf.to_np(var.coords['xy_loc'])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str() for pair in wrf.to_np(coord_pairs)]
    num_ticks = 5
    thin = int((len(x_ticks) / num_ticks) + 0.5)
    ax_cross.set_xticks(x_ticks[::thin])
    ax_cross.set_xticklabels(x_labels[::thin], rotation=10, fontsize=fontsize)

    ax_cross.set_xlabel(xlab, fontsize=fontsize+2)
    ax_cross.set_ylabel(ylab, fontsize=fontsize+2)

    ## Plot the overall title
    ax_cross.set_title(suptitle)                # This will center it on the plot (not counting colorbar or axes labels)
#   plt.suptitle(suptitle, y=suptitle_y)    # This will center it on the plot (including the colorbar & axes labels)

    ## Optional: Add titles to the subplot
    ## In this case, one title will be on the left side, one on the right side
    if title_l != None:
        plt.title(title_l, fontsize=fontsize-2, loc='left')
    if title_r != None:
        plt.title(title_r, fontsize=fontsize-2, loc='right')
    if title_c != None:
        plt.title(title_c, fontsize=fontsize-2, loc='center')

    plt.ylim(ymin, ymax)
    plt.grid(True, axis='both')

    ## Optional: Overlay wind barbs
    if x is not None and y is not None and u is not None and v is not None:
        ## Assume winds input to here are in m/s, so reduce the barb_increments from 5/10/50 to 2.5/5/25
        plt.barbs(x, y, u, v, length=6, linewidth=1.00, barb_increments={'half':2.5, 'full':5, 'flag':25})

    ## Optional: Add inset showing where the cross-section is located within the domain
    if cart_proj is not None and cart_xlim is not None and cart_ylim is not None:
#       borders, states, oceans, lakes, rivers, land = get_cartopy_features()
        if ins_loc == 'upper left':
            ax_inset = fig.add_axes([0.125, 0.668, 0.25, 0.25], projection=cart_proj)   # upper-left corner inset
        elif ins_loc == 'upper right':
            ax_inset = fig.add_axes([0.495, 0.668, 0.25, 0.25], projection=cart_proj)   # upper-right corner inset
        ax_inset.set_xlim(cart_xlim)
        ax_inset.set_ylim(cart_ylim)
        ax_inset.coastlines()
        if borders is not None:
            ax_inset.add_feature(borders, linestyle='-')
        if states is not None:
            ax_inset.add_feature(states, linewidth=0.25, edgecolor='black')
        if lakes is not None:
            ax_inset.add_feature(lakes, linewidth=0.25, edgecolor='black')
        ax_inset.set_title('')

        ## Optional: Plot a contour-filled variable field in the inset
        if ins_var is not None and wrf_lons is not None and wrf_lats is not None and ins_cmap is not None and ins_norm is not None:
            ax_inset.contourf(wrf.to_np(wrf_lons), wrf.to_np(wrf_lats), wrf.to_np(ins_var),
                ins_bounds, transform=data_crs, cmap=ins_cmap, norm=ins_norm)

        ## Optional: Plot the path of the cross section
        if lon_beg is not None and lon_end is not None and lat_beg is not None and lat_end is not None:
            ax_inset.plot([lon_beg, lon_end], [lat_beg, lat_end],
                color='black', linewidth=2, marker='o', transform=data_crs)

    plt.savefig(fname)
    plt.close()

def mpl_map_plot(fname, var, suptitle, cbar_lab, cart_proj, cart_xlim, cart_ylim, lons, lats,
        cmap, bounds, norm, extend, cbar_loc='bottom', cbar_fac=0.047,
        borders=None, states=None, oceans=None, lakes=None, water_color='none', border_width=2.0,
        lat_labels=None, lon_labels=None, polygon_verts=None, polygon_color='black',
        suptitle_y=0.95, title_l=None, title_r=None, title_c=None,
        marker_lat=None, marker_lon=None, marker='*', marker_size=100, marker_color='None', marker_edgecolor='black', marker_width=1.5,
        marker_lat2=None,marker_lon2=None,marker2='*',marker_size2=100,marker_color2='None',marker_edgecolor2='black', marker_width2=1.5,
        marker_val_fill=False, marker_val_fill2=False, marker_var=None, marker_var2=None, marker_zorder=10,
        text_lat=None, text_lon=None, text_lab=None, fontsize=14, text_lab_wt='bold',
        text_lat2=None, text_lon2=None, text_lab2=None, text_lab_wt2='normal',
        cont_var=None, cont_levs=None, cont_color='black', cont_width=0.75, cont_lab=False,
        cross_lon_beg=None, cross_lat_beg=None, cross_lon_end=None, cross_lat_end=None,
        map_x_thin=25, map_y_thin=25, u=None, v=None,
        lg_text=None, lg_loc='lower left', lg_fontsize=14,):
    ## NOTE: There are also ShapelyDeprecationWarnings that are thrown with Cartopy 0.20.1/Shapely 1.8.0.
    ##       These warnings involve issues with multi-part geometries.
    ##       These (and other) warnings can be ignored by executing the script with python3 -W ignore script.py.
    ##       These warnings were resolved with Cartopy 0.20.2.

    '''
    -- Procedure to make a map plot with filled contours using matplotlib and Cartopy.
        Optionally overlays the map with wind barbs, markers, text labels, a polygon, or cross-section path.
    -- Added by JAL on 4 Feb 2022
    -- Required Positional Inputs:
        - fname: string or pathlib object
        - var: 2-D variable array to be plotted with filled contours
        - suptitle: string for overall plot title (usually one line)
        - cbar_lab: string for colorbar label (e.g., Wind Speed [m/s])
        - cart_proj: Cartopy object, map projection
        - cart_xlim, cart_ylim: Cartopy objects, x- & y-axis limits
        - lons, lats: 1D or 2D arrays of longitude and latitude values
        - cmap: Matplotlib colormap
        - bounds: Matplotlib colormap bounds
        - norm: Matplotlib colormap norm
        - extend: string for colorbar cap ('max', 'min', 'both')
    -- Optional Inputs:
        - cbar_loc: string, identifier for positioning of the colorbar ('top', 'bottom', 'right', 'left')
        - cbar_fac: scale factor for colorbar size (default: 0.047)
        - borders, states, oceans, lakes: Cartopy feature objects for the map
        - border_width: line thickness for national borders & coastlines (default: 2.0)
        - water_color: string defining water color for the map (usually only for terrain plots; 'none'=transparent)
        - lat_labels, lon_labels: arrays of latitude & longitude values to label explicitly on the map
        - polygon_verts: vertexes of a polygon to overlay on the map
        - polygon_color: line color of the polygon (default: black)
        - suptitle_y: float, y-axis position of the suptitle
        - title_l: string, plot subtitle (1 or 2 lines) that gets placed above the top-left corner of the plot axes
        - title_r: string, plot subtitle (1 or 2 lines) that gets placed above the top-right corner of the plot axes
        - title_c: string, plot subtitle (1 or 2 lines) that gets placed above the center of the plot axes
        - marker_lat, marker_lon: arrays of the latitude and longitude values for markers
        - marker: marker style (if multiple marker types are desired, this will need to be modified to handle arrays)
        - marker_size: default marker size (default: 100)
        - marker_color: marker color ('none' means transparent)
        - marker_edgecolor: marker edge color (default: black)
        - marker_width: marker linewidth (default: 1.5)
        - marker_val_fill: flag to fill markers from data value, cmap, and norm (default: False)
        - marker_var: variable containing data to plot in the markers (use this to plot station data in a model field)
        - marker_zorder: integer indicator of the matplotlib draw order (default: 10)
        - cont_var: 2D variable array for unfilled contour overlays
        - cont_levs: specified levels for plotting of cont_var
        - cont_color, cont_width: contour overlay color (default: black) and linewidth (default: 0.5)
        - cont_lab: boolean flag to turn off/on contour labels (default: False)
        - text_lat, text_lon: arrays of the latitude and longitude values for text labels
        - text_lab: array of text labels
        - text_lab_wt: text label weight (normal, bold, heavy, light, ultrabold, ultralight)
        - fontsize: text label font size (default: 14). Other plot text sizes are anchored off this value.
        - cross_lon_beg, cross_lat_beg, cross_lon_end, cross_lat_end: Start & end lat/lon points of cross-section
        - map_x_thin, map_y_thin: integer for thinning wind barb location overlays (every Nth grid point)
        - u, v: array-like, define the barb directions
        - lg_text: array of legend labels
        - lg_loc: string defining legend placement (default: lower left)
        - lg_fontsize: integer fontsize for legend labels
   '''
    mpl.rcParams['figure.figsize'] = (10, 8)
    mpl.rcParams['grid.color'] = 'gray'
    mpl.rcParams['grid.linestyle'] = ':'
    mpl.rcParams['font.size'] = fontsize+2
    mpl.rcParams['figure.titlesize'] = fontsize+2
    mpl.rcParams['savefig.bbox'] = 'tight'
    ll_size = fontsize-2
    data_crs = ccrs.PlateCarree()

    print('-- Plotting: '+str(fname))
    fig = plt.figure()
    ## These definitions for ax are functionally equivalent, it seems...
    ## https://stackoverflow.com/questions/43326680/what-are-the-differences-between-add-axes-and-add-subplot
#   ax = plt.axes(projection=cart_proj)
    ax = plt.subplot(projection=cart_proj)

    ## If cart_xlim and cart_ylim tuples are not provided, then set plot limits from lat/lon data directly
    if cart_xlim is not None and cart_ylim is not None:
        ax.set_xlim(cart_xlim)
        ax.set_ylim(cart_ylim)
    else:
        ax.set_extent([np.min(lons), np.max(lons), np.min(lats), np.max(lats)], crs=cart_proj)

    ## Optional: Add various cartopy features
    if borders != None:
        ax.add_feature(borders, linewidth=border_width, linestyle='-', zorder=3)
    if states != None:
        ax.add_feature(states, linewidth=border_width/2.0, edgecolor='black', zorder=4)
    if oceans != None:
        ## Drawing the oceans is VERY slow with Cartopy 0.20+ for some domains, so may want to skip it
        ## Can set the facecolor for the axes to water_color instead (usually we want this 'none' except for terrain)
        ax.add_feature(oceans, facecolor=water_color, zorder=2)
#       ax.set_facecolor(water_color)
    if lakes != None:
        ## Unless facecolor='none', lakes w/ facecolor will appear above filled contour plot, which is undesireable
        ax.add_feature(lakes, facecolor=water_color, linewidth=0.25, edgecolor='black', zorder=5)
    ax.coastlines(zorder=6, linewidth=border_width)

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
            else:   # loop back around to the first vertex
                lon_end = polygon_verts[0][0]
                lat_end = polygon_verts[0][1]
            ## Plot each side of the polygon (could optionally add a marker at each vertex, but not doing that here)
            ax.plot([lon_beg, lon_end], [lat_beg, lat_end], color=polygon_color, linewidth=2, transform=data_crs)

    if lons.ndim == 1 and lats.ndim == 1:
        ll2d = np.meshgrid(lons, lats)
        lons = ll2d[0]
        lats = ll2d[1]

    ## Optional: Draw unfilled contours of another variable
    if cont_var is not None:
        if not cont_lab:
            plt.contour(wrf.to_np(lons), wrf.to_np(lats), wrf.to_np(cont_var), cont_levs,
                colors=cont_color, linewidths=cont_width, transform=data_crs, transform_first=(ax,True))
        if cont_lab:
            conts = plt.contour(wrf.to_np(lons), wrf.to_np(lats), wrf.to_np(cont_var), cont_levs,
                colors=cont_color, linewidths=cont_width, transform=data_crs, transform_first=(ax,True))
            plt.clabel(conts, conts.levels, inline=True, fontsize=fontsize)

    ## Draw the actual filled contour plot
    ## NOTE: Sometimes a cartopy contourf plot may fail with a Shapely TopologicalError.
    ##       It appears to happen sometimes when projecting a dataset onto a different projection.
    ##       This error occurs more often if nans are present.
    ##       Replacing nans with scipy.interpolate.griddata solves some of these errors, but not all.
    ##       Interestingly, plt.contour still seems to work in these situations as a (suboptimal) workaround.
    ##       Hopefully this bug will be resolved eventually.
    ##       It occurs with Shapely 1.8.0, Cartopy 0.20.1.
    ##       This issue was resolved with Cartopy 0.20.2 (https://github.com/SciTools/cartopy/issues/1936).
    ## Using the transform_first argument results in a noticeable speed-up:
    ##    https://scitools.org.uk/cartopy/docs/latest/gallery/scalar_data/contour_transforms.html
    ##  - This also requires the X and Y variables to be 2-D arrays.
    ##  - This also resolves TopologyException: side location conflict errors that can occur when NaNs are present.

    ## If the variable has the same shape as lats, then plot the filled contour field
    if (var.shape == lats.shape):
        contourf = True
        plt.contourf(wrf.to_np(lons), wrf.to_np(lats), wrf.to_np(var), bounds, cmap=cmap, norm=norm, extend=extend,
            transform=data_crs, transform_first=(ax,True))
    ## Otherwise, we presumably need to plot an empty map and then plot markers
    else:
        contourf = False
        if marker_lon is None or marker_lat is None:
            print('ERROR: mpl_map_plot in functions/map_functions.py:')
            print('   var does not match the shape of lons or lats.')
            print('   If plotting a contour map, fix the mismatch in shape between var and lons/lats.')
            print('   Otherwise, this could imply a need for a blank map with markers.')
            print('   However, marker_lon and/or marker_lat are not provided either.')
            print('   If a filled contour map is desired, ensure var is the same shape as lons and lats.')
            print('   If a blank map with markers is desired instead, provide marker_lon and marker_lat.')
            print('   If the markers should be filled according to a data value, then set marker_val_fill=True.')
            print('   Or if a future use case is identified that requires a blank map & no markers, then modify code.')
            print('   Exiting!')
            sys.exit()

    ## Optional: Add markers to the plot
    if marker_lon is not None and marker_lat is not None:
        if lg_text is None:
            lg_lab1 = None
        else:
            lg_lab1 = lg_text[0]
        if marker_val_fill:
            ## Fill markers according to their data value, cmap, and norm
            ## NOTE: Must be plt.scatter and not ax.scatter to avoid an error about the lack of a mappable when
            ##       drawing the colorbar below.
            if (var.shape == marker_lat.shape):
                if contourf:
                    ax.scatter(marker_lon, marker_lat, c=var, marker=marker, s=marker_size,
                        cmap=cmap, norm=norm,
                        edgecolors=marker_edgecolor, transform=data_crs, zorder=marker_zorder, label=lg_lab1,
                        linewidths=marker_width)
                else:
                    plt.scatter(marker_lon, marker_lat, c=var, marker=marker, s=marker_size,
                        cmap=cmap, norm=norm,
                        edgecolors=marker_edgecolor, transform=data_crs, zorder=marker_zorder, label=lg_lab1,
                        linewidths=marker_width)
            elif (marker_var.shape == marker_lat.shape):
                if contourf:
                    ax.scatter(marker_lon, marker_lat, c=marker_var, marker=marker, s=marker_size,
                        cmap=cmap, norm=norm,
                        edgecolors=marker_edgecolor, transform=data_crs, zorder=marker_zorder, label=lg_lab1,
                        linewidths=marker_width)
                else:
                    plt.scatter(marker_lon, marker_lat, c=var, marker=marker, s=marker_size,
                        cmap=cmap, norm=norm,
                        edgecolors=marker_edgecolor, transform=data_crs, zorder=marker_zorder, label=lg_lab1,
                        linewidths=marker_width)
        else:
            ## These both seem functionally equivalent:
#           plt.scatter(marker_lon, marker_lat, marker=marker, color=marker_color, s=marker_size,
#               edgecolors=marker_edgecolor, transform=data_crs, zorder=10)
            ax.scatter(marker_lon, marker_lat, marker=marker, color=marker_color, s=marker_size,
                edgecolors=marker_edgecolor, transform=data_crs, zorder=marker_zorder, label=lg_lab1,
                linewidths=marker_width)

    ## Optional: Add a second set of markers to the plot
    if marker_lon2 is not None and marker_lat2 is not None:
        if lg_text is None:
            lg_lab2 = None
        else:
            lg_lab2 = lg_text[1]
        if marker_val_fill2:
            ax.scatter(marker_lon2, marker_lat2, c=marker_var2, marker=marker2, s=marker_size2,
                cmap=cmap, norm=norm,
                edgecolors=marker_edgecolor2, transform=data_crs, zorder=marker_zorder, label=lg_lab2,
                linewidths=marker_width2)
        else:
            ax.scatter(marker_lon2, marker_lat2, marker=marker2, color=marker_color2, s=marker_size2,
                edgecolors=marker_edgecolor2, transform=data_crs, zorder=marker_zorder, label=lg_lab2,
                linewidths=marker_width2)

    ## Draw the colorbar
    '''
    ## Want to ensure the colorbar length matches the plot length (for horizontal cbar) or height (for vertical cbar)
    ## First, calculate the proper ratio (width/height for horizontal colorbar, height/width for vertical colorbar)
    ## See: https://www.geeksforgeeks.org/set-matplotlib-colorbar-size-to-match-graph/
    ## NOTE: This method is not robust for all domain sizes & aspect ratios!
    if cbar_loc == 'top' or cbar_loc == 'bottom':
        ## Horizontal colorbar
        im_ratio = var.shape[1] / var.shape[0]
    elif cbar_loc == 'left' or cbar_loc == 'right':
        ## Vertical colorbar
        im_ratio = var.shape[0] / var.shape[1]
    plt.colorbar(ax=ax, location=cbar_loc, fraction=cbar_fac*im_ratio, label=cbar_lab, extend=extend, pad=0.05)
    '''
    ## Credit: https://stackoverflow.com/questions/30030328/correct-placement-of-colorbar-relative-to-geo-axes-cartopy
    ## Create colorbar axes (temporarily) anywhere
    cax = fig.add_axes([0,0,0.1,0.1])
    ## Find the location of the main plot axes
    posn = ax.get_position()
    ## Adjust the positioning and orientation of the colorbar and then draw it
    ## The colorbar will inherit the extend property from the plt.contourf call above and its norm/extend attributes
    if cbar_loc == 'bottom':
        cax.set_position([posn.x0, posn.y0-0.09, posn.width, 0.05])
        plt.colorbar(cax=cax, orientation='horizontal', label=cbar_lab)
    elif cbar_loc == 'right':
        cax.set_position([posn.x0+posn.width+0.05, posn.y0, 0.04, posn.height])
        plt.colorbar(cax=cax, orientation='vertical', label=cbar_lab)
    elif cbar_loc == 'top' or cbar_loc == 'left':
        print('WARNING: cbar_loc='+cbar_loc+' requested. Unsupported option. Colorbar will not be drawn.')
        print('   Add directives in map_functions.mpl_map_plot to handle that option and draw the colorbar.')

#   ## Optional: Draw unfilled contours of another variable
#   if cont_var is not None:
#       plt.contour(wrf.to_np(lons), wrf.to_np(lats), wrf.to_np(cont_var), cont_levs,
#           transform=data_crs, transform_first=(ax,True))

    ## Add the plot overall title
    plt.suptitle(suptitle, y=suptitle_y)

    ## Optional: Add titles to the subplot
    ## In this case, one title will be on the left side, one on the right side
    if title_l != None:
#       plt.title(title_l, fontsize=10, loc='left')
        ax.set_title(title_l, fontsize=fontsize, loc='left')
    if title_r != None:
#       plt.title(title_r, fontsize=10, loc='right')
        ax.set_title(title_r, fontsize=fontsize, loc='right')
    if title_c != None:
#       plt.title(title_c, fontsize=10, loc='center')
        ax.set_title(title_c, fontsize=fontsize, loc='center')

    ## Optional: Add text labels to the plot
    if text_lab is not None and text_lat is not None and text_lon is not None:
        n_text = len(text_lab)
        for xx in range(n_text):
#           plt.text(text_lon[xx], text_lat[xx], text_lab[xx], horizontalalignment='center',
#               transform=data_crs, size=10, zorder=11)
            ax.text(text_lon[xx], text_lat[xx], text_lab[xx], horizontalalignment='center',
                transform=data_crs, size=fontsize, zorder=11, weight=text_lab_wt)

    ## Optional: Add a second set of text labels to the plot
    if text_lab2 is not None and text_lat2 is not None and text_lon2 is not None:
        n_text = len(text_lab2)
        for xx in range(n_text):
            ax.text(text_lon2[xx], text_lat2[xx], text_lab2[xx], horizontalalignment='center',
                transform=data_crs, size=fontsize, zorder=11, weight=text_lab_wt2)

    ## Optional: Plot the path of the cross section
    if cross_lon_beg is not None and cross_lon_end is not None and cross_lat_beg is not None and cross_lat_end is not None:
#       plt.plot([cross_lon_beg, cross_lon_end], [cross_lat_beg, cross_lat_end],
#           color='black', linewidth=2, marker='o', transform=data_crs)
        ax.plot([cross_lon_beg, cross_lon_end], [cross_lat_beg, cross_lat_end],
            color='black', linewidth=2, marker='o', transform=data_crs)

    ## Optional: Draw wind barbs
    if u is not None and v is not None:
        if isinstance(lons, np.ndarray):
            x_thin = lons[::map_y_thin, ::map_x_thin]
        else:
            x_thin = lons[::map_y_thin, ::map_x_thin].values
        if isinstance(lats, np.ndarray):
            y_thin = lats[::map_y_thin, ::map_x_thin]
        else:
            y_thin = lats[::map_y_thin, ::map_x_thin].values
        u_thin = u[::map_y_thin, ::map_x_thin]
        v_thin = v[::map_y_thin, ::map_x_thin]
        ## Assume winds input to here are in m/s, so reduce the barb_increments from 5/10/50 to 2.5/5/25
#       plt.barbs(x_thin, y_thin, u_thin, v_thin, length=5, transform=data_crs, linewidth=0.25,
#           barb_increments={'half':2.5, 'full':5, 'flag':25})
        ax.barbs(x_thin, y_thin, u_thin, v_thin, length=5, transform=data_crs, linewidth=0.25,
            barb_increments={'half':2.5, 'full':5, 'flag':25})

    ## Optional: Add a legend (useful if 2+ sets of markers)
    if lg_text is not None:
        ax.legend(loc=lg_loc, fontsize=lg_fontsize).set_zorder(12)

    ## Save and close the figure
    ## NOTE: This may throw a Shapely TopologicalError with cartopy 0.20.0 & 0.20.1. Fixed with 0.20.2. (see above)
    plt.savefig(fname)
    plt.close()


def mpl_map_stations(fname, suptitle, cart_proj, cart_xlim, cart_ylim, lons, lat,
    borders=None, states=None, oceans=None, lakes=None, water_color='lightblue',
    lat_labels=None, lon_labels=None, suptitle_y=0.95, fontsize=14,
    mark1_lon=None, mark1_lat=None, mark1='o', mark1_size=81, mark1_color='tab:blue',   mark1_edgecolor='black',
    mark2_lon=None, mark2_lat=None, mark2='s', mark2_size=81, mark2_color='tab:orange', mark2_edgecolor='black',
    mark3_lon=None, mark3_lat=None, mark3='v', mark3_size=81, mark3_color='tab:green',  mark3_edgecolor='black',
    mark4_lon=None, mark4_lat=None, mark4='^', mark4_size=81, mark4_color='tab:red',    mark4_edgecolor='black',
    mark5_lon=None, mark5_lat=None, mark5='P', mark5_size=81, mark5_color='tab:purple', mark5_edgecolor='black',
    mark6_lon=None, mark6_lat=None, mark6='p', mark6_size=81, mark6_color='tab:brown',  mark6_edgecolor='black',
    mark7_lon=None, mark7_lat=None, mark7='*', mark7_size=81, mark7_color='tab:pink',   mark7_edgecolor='black',
    mark8_lon=None, mark8_lat=None, mark8='H', mark8_size=81, mark8_color='tab:gray',   mark8_edgecolor='black',
    mark9_lon=None, mark9_lat=None, mark9='X', mark9_size=81, mark9_color='tab:olive',  mark9_edgecolor='black',
    mark10_lon=None,mark10_lat=None,mark10='D',mark10_size=81,mark10_color='tab:cyan',  mark10_edgecolor='black',
    mark11_lon=None,mark11_lat=None,mark11='8',mark11_size=81,mark11_color='black',     mark11_edgecolor='gold',
    lg_text=None, lg_loc='lower left',lg_fontsize=14):

    '''
    -- Procedure to make a map plot with station markers, using matplotlib and Cartopy.
    -- Added by JAL on 5 Jan 2022
    '''

    mpl.rcParams['figure.figsize'] = (10, 8)
    mpl.rcParams['grid.color'] = 'gray'
    mpl.rcParams['grid.linestyle'] = ':'
    mpl.rcParams['font.size'] = fontsize+2
    mpl.rcParams['figure.titlesize'] = fontsize+2
    mpl.rcParams['legend.framealpha'] = 1.0
    mpl.rcParams['savefig.bbox'] = 'tight'
    ll_size = fontsize-2
#    lg_fontsize = fontsize
    data_crs = ccrs.PlateCarree()

    print('Plotting '+str(fname))
    fig = plt.figure()
    ax = plt.subplot(projection=cart_proj)

    ## If cart_xlim and cart_ylim tuples are not provided, then set plot limits from lat/lon data directly
    if cart_xlim is not None and cart_ylim is not None:
        ax.set_xlim(cart_xlim)
        ax.set_ylim(cart_ylim)
    else:
        ax.set_extent([np.min(lons), np.max(lons), np.min(lats), np.max(lats)], crs=cart_proj)

    ## If only 1-D arrays are provided for lons & lats, make them into 2-D arrays
    if lons.ndim == 1 and lats.ndim == 1:
        ll2d = np.meshgrid(lons, lats)
        lons = ll2d[0]
        lats = ll2d[1]

    ## Optional: Add various cartopy features
    if borders != None:
        ax.add_feature(borders, linestyle='-', zorder=3)
    if states != None:
        ax.add_feature(states, linewidth=0.25, edgecolor='black', zorder=4)
    if oceans != None:
        ax.add_feature(oceans, facecolor=water_color, zorder=2)
    if lakes != None:
        ax.add_feature(lakes, facecolor=water_color, linewidth=0.25, edgecolor='black', zorder=5)
    ax.coastlines(zorder=6)

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

    ## Optional: Add a set of markers to the plot
    if mark1_lon is not None and mark1_lat is not None:
        ax.scatter(mark1_lon, mark1_lat, marker=mark1, color=mark1_color, s=mark1_size,
            edgecolors=mark1_edgecolor, transform=data_crs, zorder=10)

    ## Optional: Add a 2nd set of markers to the plot
    if mark2_lon is not None and mark2_lat is not None:
        ax.scatter(mark2_lon, mark2_lat, marker=mark2, color=mark2_color, s=mark2_size,
            edgecolors=mark2_edgecolor, transform=data_crs, zorder=11)

    ## Optional: Add a 3rd set of markers to the plot
    if mark3_lon is not None and mark3_lat is not None:
        ax.scatter(mark3_lon, mark3_lat, marker=mark3, color=mark3_color, s=mark3_size,
            edgecolors=mark3_edgecolor, transform=data_crs, zorder=12)

    ## Optional: Add a 4th set of markers to the plot
    if mark4_lon is not None and mark4_lat is not None:
        ax.scatter(mark4_lon, mark4_lat, marker=mark4, color=mark4_color, s=mark4_size,
            edgecolors=mark4_edgecolor, transform=data_crs, zorder=13)

    ## Optional: Add a 5th set of markers to the plot
    if mark5_lon is not None and mark5_lat is not None:
        ax.scatter(mark5_lon, mark5_lat, marker=mark5, color=mark5_color, s=mark5_size,
            edgecolors=mark5_edgecolor, transform=data_crs, zorder=14)

    ## Optional: Add a 6th set of markers to the plot
    if mark6_lon is not None and mark6_lat is not None:
        ax.scatter(mark6_lon, mark6_lat, marker=mark6, color=mark6_color, s=mark6_size,
            edgecolors=mark6_edgecolor, transform=data_crs, zorder=15)

    ## Optional: Add a 7th set of markers to the plot
    if mark7_lon is not None and mark7_lat is not None:
        ax.scatter(mark7_lon, mark7_lat, marker=mark7, color=mark7_color, s=mark7_size,
            edgecolors=mark7_edgecolor, transform=data_crs, zorder=16)

    ## Optional: Add a 8th set of markers to the plot
    if mark8_lon is not None and mark8_lat is not None:
        ax.scatter(mark8_lon, mark8_lat, marker=mark8, color=mark8_color, s=mark8_size,
            edgecolors=mark8_edgecolor, transform=data_crs, zorder=17)

    ## Optional: Add a 9th set of markers to the plot
    if mark9_lon is not None and mark9_lat is not None:
        ax.scatter(mark9_lon, mark9_lat, marker=mark9, color=mark9_color, s=mark9_size,
            edgecolors=mark9_edgecolor, transform=data_crs, zorder=18)

    ## Optional: Add a 10th set of markers to the plot
    if mark10_lon is not None and mark10_lat is not None:
        ax.scatter(mark10_lon, mark10_lat, marker=mark10, color=mark10_color, s=mark10_size,
            edgecolors=mark10_edgecolor, transform=data_crs, zorder=19)

    ## Optional: Add an 11th set of markers to the plot
    if mark11_lon is not None and mark11_lat is not None:
        ax.scatter(mark11_lon, mark11_lat, marker=mark11, color=mark11_color, s=mark11_size,
            edgecolors=mark11_edgecolor, transform=data_crs, zorder=20)

    ## Optional: Add a legend (useful if 2+ sets of markers)
    if lg_text is not None:
        plt.legend(lg_text, loc=lg_loc, fontsize=lg_fontsize)

    plt.suptitle(suptitle, y=suptitle_y)

    plt.savefig(str(fname))
    plt.close()

    
def make_skewt(fname, title, p1, t1, td1, ws1, wd1, ht1, color1='black', label1=None, isdropsonde1=False,
        p2=None, t2=None, td2=None, ws2=None, wd2=None, ht2=None, color2='magenta', label2=None, isdropsonde2=False,
        lcl_dot=True, fontsize=14):

    '''
    -- Procedure to make a skew-T plot using MetPy. Optionally overlay a second sounding for comparison.
    -- Added by JAL on 2 Nov 2022
    -- Required Positional Inputs:
        - fname: file name (string or Pathlib object)
        - title: plot title (string)
        - p1: pressure [hPa] (1-D numeric array)
        - t1: temperature [degC] (1-D numeric array)
        - td1: dewpoint [degC] (1-D numeric array)
        - ws1: wind speed [m/s] (1-D numeric array)
        - wd1: wind direction [deg] (1-D numeric array)
        - ht1: height/altitude [m] (1-D numeric array)
   -- Optional Inputs:
        - color1: color for 1st skew-T dataset (string)
        - color2: color for 2nd skew-T dataset (string)
        - label1: label for 1st skew-T dataset (string)
        - label2: label for 2nd skew-T dataset (string)
        - isdropsonde1: indicator for whether the 1st skew-T dataset is a dropsonde/aircraft descent (boolean)
        - isdropsonde2: indicator for whether the 2nd skew-T dataset is a dropsonde/aircraft descent (boolean)
        - p2: pressure [hPa] for 2nd skew-T dataset (1-D numeric array)
        - t2: temperature [degC] for 2nd skew-T dataset (1-D numeric array)
        - td2: dewpoint [degC] for 2nd skew-T dataset (1-D numeric array)
        - ws2: wind speed [m/s] for 2nd skew-T dataset (1-D numeric array)
        - wd2: wind direction [deg] for 2nd skew-T dataset (1-D numeric array)
        - ht2: height/altitude [m] for 2nd skew-T dataset (1-D numeric array)
        - lcl_dot: indicator for whether to plot the LCL as a dot on the skew-T (boolean)
        - fontsize: fontsize off of which all other font sizes are based (numeric)
   -- Output: A figure saved to a file.
    '''

    plt.rcParams['figure.figsize'] = (10,8)
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['axes.titlesize'] = fontsize+4
    plt.rcParams['axes.labelsize'] = fontsize+2
    plt.rcParams['xtick.labelsize'] = fontsize+2
    plt.rcParams['ytick.labelsize'] = fontsize+2
    plt.rcParams['savefig.bbox'] = 'tight'

    deg_uni = '\u00B0'

    ## Need to attach units and make these type Quantity to be able to plot a skew-T with MetPy
    pres1 = p1 * mpunits.hPa
    temp1 = t1 * mpunits.degC
    dewp1 = td1 * mpunits.degC
    wspd1 = ws1 * mpunits('m/s')
    wdir1 = wd1 * mpunits.degrees
    alt1 = ht1 * mpunits.meters

    if p2 is not None:
        pres2 = p2 * mpunits.hPa
    if t2 is not None:
        temp2 = t2 * mpunits.degC
    if td2 is not None:
        dewp2 = td2 * mpunits.degC
    if ws2 is not None:
        wspd2 = ws2 * mpunits('m/s')
    if wd2 is not None:
        wdir2 = wd2 * mpunits.degrees
    if ht2 is not None:
        alt2 = ht2 * mpunits.meters

    if p2 is None and t2 is None and td2 is None:
        plot_skewt2 = False
        color_t1 = 'tab:red'
        color_td1 = 'tab:green'
        style_t = '-'
        style_td = '-'
        width_t = 2.0
        width_td = 2.0
    elif p2 is not None and t2 is not None and td2 is not None:
        plot_skewt2 = True
        color_t1 = color1
        color_td1 = color1
        color_t2 = color2
        color_td2 = color2
        style_t = '-'
        style_td = '--'
        width_t = 2.0
        width_td = 1.5
    else:
        print('ERROR: Some but not all of p2, t2, and td2 provided to function make_skewt. Exiting!')
        sys.exit()

    ## Need wind components for MetPy wind barbs
    u1, v1 = mpcalc.wind_components(wspd1, wdir1)
    if ws2 is not None and wd2 is not None:
        u2, v2 = mpcalc.wind_components(wspd2, wdir2)

    print('   Creating plot '+str(fname))
    fig = plt.figure()
    skew = mpplots.SkewT(fig)

    ## Plot 0-degree isotherm
    skew.ax.axvline(0, color='cyan', linestyle='-', linewidth=1)

    ## Plot temperature and dewpoint lines
    skew.plot(pres1, temp1, color=color_t1, linestyle=style_t, linewidth=width_t, label=label1)
    skew.plot(pres1, dewp1, color=color_td1, linestyle=style_td, linewidth=width_td)
    if plot_skewt2:
        skew.plot(pres2, temp2, color=color_t2, linestyle=style_t, linewidth=width_t, label=label2)
        skew.plot(pres2, dewp2, color=color_td2, linestyle=style_td, linewidth=width_td)

    ## Plot wind barbs on the right side of the skew-T
    skew.plot_barbs(pres1, u1, v1, color=color1, y_clip_radius=0.01)

    ## For some reason I am unable to change the color of these lines, so adjust other line properties instead
    skew.plot_dry_adiabats(linewidth=0.5, linestyle=':', alpha=1.0)
    skew.plot_moist_adiabats(linewidth=0.5, linestyle=':', alpha=1.0)
    skew.plot_mixing_lines(linewidth=0.5, linestyle=':', alpha=1.0)
    skew.ax.set_ylim(1000, 100)

    ## Calculate LCL and plot as dot
    ## NOTE: Indexing assumes array element 0 is the bottom of the sounding. Dropsondes will need to be inverted.
    if lcl_dot:
        if isdropsonde1:
            lcl_p1, lcl_t1 = mpcalc.lcl(pres1[-1], temp1[-1], dewp1[-1])
        else:
            lcl_p1, lcl_t1 = mpcalc.lcl(pres1[0], temp1[0], dewp1[0])
        skew.plot(lcl_p1, lcl_t1, color=color1, marker='o', markerfacecolor=color1)
        if plot_skewt2:
            if isdropsonde2:
                lcl_p2, lcl_t2 = mpcalc.lcl(pres2[-1], temp2[-1], dewp2[-1])
            else:
                lcl_p2, lcl_t2 = mpcalc.lcl(pres2[0], temp2[0], dewp2[0])
            skew.plot(lcl_p2, lcl_t2, color=color2, marker='o', markerfacecolor=color2)

    ## Add secondary axis that automatically converts between pressure and height assuming a standard atmosphere.
    ## Source: https://github.com/Unidata/MetPy/issues/262
    ax2 = skew.ax.secondary_yaxis(1.07,
        functions=( lambda pres1: mpcalc.pressure_to_height_std(mpunits.Quantity(pres1, 'hPa')).m_as('km'),
                        lambda alt1: mpcalc.height_to_pressure_std(mpunits.Quantity(alt1, 'km')).m))
    ax2.yaxis.set_major_locator(plt.FixedLocator([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]))
    ax2.yaxis.set_minor_locator(plt.NullLocator())
    ax2.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax2.set_ylabel('Height [km]')

    plt.title(title)
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Temperature ['+deg_uni+'C]')
    if plot_skewt2:
        plt.legend(loc='upper right', framealpha=1.0)
    plt.savefig(fname)
    plt.close()


## NOTE: Commented this out on 9 Feb 2022 after switching to conda-based NCAR Python env. Ngl not installed yet.
#def make_ngl_map_plot_poly(plot_type, plot_file, mpres, poly_lon, poly_lat, polyres,
#                                   poly_lab=[''], text_lon=[np.nan], text_lat=[np.nan], textres=Ngl.Resources(),
#                                   font_ht=0.020, l_str='', l_x=0.25, l_y=0.95, r_str='', r_x=0.78, r_y=0.95):
    '''
   -- Procedure to make an empty map plot with polymarkers and optional text labels.
   -- Added by JAL on 29 Oct 2020
   -- Required Positional Inputs:
      - plot_type: string
      - plot_file: pathlib.Path object or string
      - mpres: Ngl.Resources object
      - poly_lon: 1-D numeric array
      - poly_lat: 1-D numeric array
      - polyres: Ngl.Resources object
    -- Optional Inputs:
      - poly_lab: 1-D string array (default '')
      - text_lat: 1-D numeric array (default np.nan)
      - text_lon: 1-D numeric array (default np.nan)
      - textres: Ngl.Resources object, default Ngl.Resources())
      - font_ht: float (default 0.020), NDC units for font height for Ngl.add_text
      - l_str: string (default ''), left-side plot string (replacement for gsnLeftString in NCL)
      - l_x: float (default 0.25), NDC units for x-position of l_str
      - l_y: float (default 0.95), NDC units for y-position of l_str
      - r_str: string (default ''), right-side plot string (replacement for gsnRightString in NCL)
      - r_x: float (default 0.78), NDC units for x-position of r_str
      - r_y: float (default 0.95), NDC units for y-position of r_str
   -- Output: None.
    '''

#   print('-- Creating plot: '+str(plot_file)+'.'+plot_type)
#   wks = Ngl.open_wks(plot_type, str(plot_file))
#   plot = Ngl.map(wks, mpres)
#   poly = Ngl.add_polymarker(wks, plot, poly_lon, poly_lat, polyres)
#   ## Ngl.add_Text can only work on scalar values, not arrays/lists, so we need to loop
#   for i in range(len(poly_lab)):
#       if not (poly_lab[i] == '' or poly_lab[i] == ' '):
#           text  = Ngl.add_text(wks, plot, poly_lab[i], text_lon[i], text_lat[i], textres)
#           Ngl.draw(plot)

#   ngl_strings(wks, plot, left=l_str, right=r_str, font_ht=font_ht)
#   Ngl.draw(plot)
#   Ngl.frame(wks)
#   Ngl.delete_wks(wks)


def make_ngl_contour_map_plot(plot_type, plot_file, var, mpres,
                                        font_ht=0.020, l_str='', l_x=0.25, l_y=0.95, r_str='', r_x=0.78, r_y=0.95):
    '''
   -- Procedure to make a plot using Ngl.contour_map.
   -- Added by JAL on 13 May 2019
   -- Required Positional Inputs:
      - plot_type: string
      - plot_file: pathlib.Path object or string
      - var: 2-D numeric array
      - mpres: Ngl.Resources object
    -- Optional Inputs:
      - font_ht: float (default 0.020), NDC units for font height for Ngl.add_text
      - l_str: string (default ''), left-side plot string (replacement for gsnLeftString in NCL)
      - l_x: float (default 0.25), NDC units for x-position of l_str
      - l_y: float (default 0.95), NDC units for y-position of l_str
      - r_str: string (default ''), right-side plot string (replacement for gsnRightString in NCL)
      - r_x: float (default 0.78), NDC units for x-position of r_str
      - r_y: float (default 0.95), NDC units for y-position of r_str
   -- Output: None.
    '''

    print('-- Creating plot: '+str(plot_file)+'.'+plot_type)
#   tires   = Ngl.Resources()
#   tires.txFontHeightF = font_ht
    wks = Ngl.open_wks(plot_type, str(plot_file))
    plot    = Ngl.contour_map(wks, var, mpres)
    ngl_strings(wks, plot, left=l_str, right=r_str, font_ht=font_ht)
#   Ngl.text_ndc(wks, l_str, l_x, l_y, tires)
#   Ngl.text_ndc(wks, r_str, r_x, r_y, tires)
    if mpres.nglDraw == False:
        Ngl.draw(plot)
    if mpres.nglFrame == False:
        Ngl.frame(wks)
    Ngl.delete_wks(wks)

def make_ngl_contour_map_plot_vec(plot_type, plot_file, var, mpres, u, v, vcres,
                                             font_ht=0.020, l_str='', l_x=0.25, l_y=0.95, r_str='', r_x=0.78, r_y=0.95):
    '''
   -- Procedure to make a plot using Ngl.contour_map, overlaid with wind vectors using Ngl.contour_map & Ngl.overlay
   -- Added by JAL on 13 May 2019
   -- Required Positional Inputs:
      - plot_type: string
      - plot_file: pathlib.Path object or string
      - var: 2-D numeric array
      - mpres: Ngl.Resources object
      - u: 2-D numeric array (zonal wind component)
      - v: 2-D numeric array (meridional wind component)
      - vcres: Ngl.Resources object
    -- Optional Inputs:
      - font_ht: float (default 0.020), NDC units for font height for Ngl.add_text
      - l_str: string (default ''), left-side plot string (replacement for gsnLeftString in NCL)
      - l_x: float (default 0.25), NDC units for x-position of l_str
      - l_y: float (default 0.95), NDC units for y-position of l_str
      - r_str: string (default ''), right-side plot string (replacement for gsnRightString in NCL)
      - r_x: float (default 0.78), NDC units for x-position of r_str
      - r_y: float (default 0.95), NDC units for y-position of r_str
   -- Output: None.
    '''

    print('-- Creating plot: '+str(plot_file)+'.'+plot_type)
#   tires   = Ngl.Resources()
#   tires.txFontHeightF = font_ht
    wks = Ngl.open_wks(plot_type, str(plot_file))

    ## To overlay wind vectors, ensure nglDraw and nglFrame are False
    vcres.nglDraw   = False
    vcres.nglFrame  = False
    mpres.nglDraw   = False
    mpres.nglFrame  = False

    vcplot  = Ngl.vector_map(wks, u, v, vcres)
    plot        = Ngl.contour_map(wks, var, mapres)
    Ngl.overlay(plot, vcplot)
    ngl_strings(wks, plot, left=l_str, right=r_str, font_ht=font_ht)
#   Ngl.text_ndc(wks, l_str, l_x, l_y, tires)
#   Ngl.text_ndc(wks, r_str, r_x, r_y, tires)
    Ngl.draw(plot)
    Ngl.frame(wks)
    Ngl.delete_wks(wks)

def make_ngl_contour_map_plot_poly(plot_type, plot_file, var, mpres, poly_lon, poly_lat, polyres, poly_lab, text_lon, text_lat, textres,
                                              font_ht=0.020, l_str='', l_x=0.25, l_y=0.95, r_str='', r_x=0.78, r_y=0.95):
    '''
   -- Procedure to make a plot using Ngl.contour_map, with optional polymarkers & text added.
   -- Added by JAL on 13 May 2019
   -- Required Positional Inputs:
      - plot_type: string
      - plot_file: pathlib.Path object or string
      - var: 2-D numeric array
      - mpres: Ngl.Resources object
      - poly_lon: 1-D numeric array
      - poly_lat: 1-D numeric array
      - polyres: Ngl.Resources object
      - poly_lab: 1-D string array
      - text_lat: 1-D numeric array
      - text_lon: 1-D numeric array
      - textres: Ngl.Resources object
    -- Optional Inputs:
      - font_ht: float (default 0.020), NDC units for font height for Ngl.add_text
      - l_str: string (default ''), left-side plot string (replacement for gsnLeftString in NCL)
      - l_x: float (default 0.25), NDC units for x-position of l_str
      - l_y: float (default 0.95), NDC units for y-position of l_str
      - r_str: string (default ''), right-side plot string (replacement for gsnRightString in NCL)
      - r_x: float (default 0.78), NDC units for x-position of r_str
      - r_y: float (default 0.95), NDC units for y-position of r_str
   -- Output: None.
    '''

    print('-- Creating plot: '+str(plot_file)+'.'+plot_type)
#   tires   = Ngl.Resources()
#   tires.txFontHeightF = font_ht
    wks = Ngl.open_wks(plot_type, str(plot_file))

    ## To draw polymarkers and text, ensure nglDraw and nglFrame are False
    mpres.nglDraw   = False
    mpres.nglFrame  = False

    plot    = Ngl.contour_map(wks, var, mpres)
    poly    = Ngl.add_polymarker(wks, plot, poly_lon, poly_lat, polyres)
    ## Ngl.add_Text can only work on scalar values, not arrays/lists, so we need to loop
    for i in range(len(poly_lab)):
        if not (poly_lab[i] == '' or poly_lab[i] == ' '):
            text    = Ngl.add_text(wks, plot, poly_lab[i], text_lon[i], text_lat[i], textres)
            Ngl.draw(plot)

    ngl_strings(wks, plot, left=l_str, right=r_str, font_ht=font_ht)
#   Ngl.text_ndc(wks, l_str, l_x, l_y, tires)
#   Ngl.text_ndc(wks, r_str, r_x, r_y, tires)
    Ngl.draw(plot)
    Ngl.frame(wks)
    Ngl.delete_wks(wks)

def make_ngl_contour_map_plot_poly_vec(plot_type, plot_file, var, mpres, poly_lon, poly_lat, polyres, poly_lab, text_lon, text_lat, textres,
                                                    u, v, vcres,
                                                    font_ht=0.020, l_str='', l_x=0.25, l_y=0.95, r_str='', r_x=0.78, r_y=0.95):
    '''
    make_ngl_contour_map_plot_poly_vec
   -- Procedure to make a plot using Ngl.contour_map, overlaid with wind vectors using
        Ngl.contour_map & Ngl.overlay, with optional polymarkers & text added.
   -- Added by JAL on 13 May 2019
   -- Required Positional Inputs:
      - plot_type: string
      - plot_file: pathlib.Path object or string
      - var: 2-D numeric array
      - mpres: Ngl.Resources object
      - poly_lon: 1-D numeric array
      - poly_lat: 1-D numeric array
      - polyres: Ngl.Resources object
      - poly_lab: 1-D string array
      - text_lat: 1-D numeric array
      - text_lon: 1-D numeric array
      - textres: Ngl.Resources object
      - u: 2-D numeric array (zonal wind component)
      - v: 2-D numeric array (meridional wind component)
      - vcres: Ngl.Resources object
    -- Optional Inputs:
      - font_ht: float (default 0.020), NDC units for font height for Ngl.add_text
      - l_str: string (default ''), left-side plot string (replacement for gsnLeftString in NCL)
      - l_x: float (default 0.25), NDC units for x-position of l_str
      - l_y: float (default 0.95), NDC units for y-position of l_str
      - r_str: string (default ''), right-side plot string (replacement for gsnRightString in NCL)
      - r_x: float (default 0.78), NDC units for x-position of r_str
      - r_y: float (default 0.95), NDC units for y-position of r_str
   -- Output: None.
    '''

    print('-- Creating plot: '+str(plot_file)+'.'+plot_type)
#   tires   = Ngl.Resources()
#   tires.txFontHeightF = font_ht
    wks = Ngl.open_wks(plot_type, str(plot_file))

    ## To draw polymarkers, text, and vectors, ensure nglDraw and nglFrame are False
    vcres.nglDraw   = False
    vcres.nglFrame  = False
    mpres.nglDraw   = False
    mpres.nglFrame  = False

    vcplot  = Ngl.vector(wks, u, v, vcres)
    plot        = Ngl.contour_map(wks, var, mpres)
    Ngl.overlay(plot, vcplot)

    poly    = Ngl.add_polymarker(wks, plot, poly_lon, poly_lat, polyres)
    ## Ngl.add_text can only work on scalar values, not arrays/lists, so we need to loop
    for i in range(len(poly_lab)):
        if not (poly_lab[i] == '' or poly_lab[i] == ' '):
            text    = Ngl.add_text(wks, plot, poly_lab[i], text_lon[i], text_lat[i], textres)
            Ngl.draw(plot)

    ngl_strings(wks, plot, left=l_str, right=r_str, font_ht=font_ht)
#   Ngl.text_ndc(wks, l_str, l_x, l_y, tires)
#   Ngl.text_ndc(wks, r_str, r_x, r_y, tires)
    Ngl.draw(plot)
    Ngl.frame(wks)
    Ngl.delete_wks(wks)

def ngl_strings(wks, plot, left='', center='', right='', font_ht=0.016):
    '''
   -- Adds annotations (left, right, or center) above a PyNGL plot.
      Corresponds to NCL's gsnLeftString, gsnCenterString, gsnRightString.
   -- Adapted from Karin Meier-Fleischer (NCAR/pyngl github issue #11)
   -- Added by JAL on 17 Oct 2019
   -- Required Positional Inputs:
      - wks: Ngl workstation object
      - plot: Ngl plot object
    -- Optional Inputs:
      - left: string (default ''), left string text
      - center: string (default ''), center string text
      - right: string (default ''), right string text
      - font_ht: float, (default 0.016), font height in NDC coordinates
    '''

    assert str(getattr(wks,'__class__')  == "<class 'int'>"), 'ERROR - 1st parameter is not a Ngl wks'
    assert str(getattr(plot,'__class__') == "<class 'ngl.PlotIds'>"), 'ERROR - 2nd parameter is not a Ngl plot'

    vpx = Ngl.get_float(plot,"vpXF")         #-- retrieve value of res.vpXF from plot
    vpy = Ngl.get_float(plot,"vpYF")         #-- retrieve value of res.vpYF from plot
    vpw = Ngl.get_float(plot,"vpWidthF")     #-- retrieve value of res.vpWidthF from plot

    txres = Ngl.Resources()
    txres.txFontHeightF = font_ht            #-- font size for left, center and right string

    y = vpy + 0.025                          #-- y-position

    if(left != ''):
        txres.txJust = "CenterLeft"           #-- text justification
        x = vpx                               #-- x-position
        Ngl.text_ndc(wks, left, x, y, txres)  #-- add text to wks

    if(center != ''):
        txres.txJust = "CenterCenter"         #-- text justification
        x = vpx + vpw/2
        Ngl.text_ndc(wks, center, x, y, txres)#-- add text to wks

    if(right != ''):
        txres.txJust = "CenterRight"          #-- text justification
        x = vpx+vpw                           #-- x-position
        Ngl.text_ndc(wks, right, x, y, txres) #-- add text to wks

def set_ngl_map_res(DS, model, lat, lon):
    '''
   -- Sets numerous basic resources for a contour map plot.
   -- Added by JAL on 4 May 2019
   -- Edited by JAL on 12 May 2019 (changed input arguments)
   -- Inputs:
      - DS: xarray dataset
      - model: string
      - lat: 1-D or 2-D numeric array
      - lon: 1-D or 2-D numeric array
   -- Output: Ngl.Resources object.
    '''

    if model == 'wrf':
        dx              = DS.attrs['DX']
        true_lat1   = DS.attrs['TRUELAT1']
        true_lat2   = DS.attrs['TRUELAT2']
        stand_lon   = DS.attrs['STAND_LON']
        center_lat  = DS.attrs['CEN_LAT']
        center_lon  = DS.attrs['CEN_LON']
        pole_lat        = DS.attrs['POLE_LAT']
        pole_lon        = DS.attrs['POLE_LON']
        map_proj        = DS.attrs['MAP_PROJ']

        n_lat   = len(lat[:,0])
        n_lon   = len(lon[0,:])

        if map_proj  == 0:
            projection  = 'None'
        elif map_proj   == 1:
            projection  = 'Lambert Conformal'
        elif map_proj   == 2:
            projection  = 'Polar'
        elif map_proj   == 3:
            projection  = 'Mercator'
        elif map_proj   == 5:
            projection  = 'Cylindrical'
        elif map_proj   == 6:
            projection  = 'Cassini'
    else:
        projection  = 'Cylindrical'

    mpres = Ngl.Resources()

    ## Core Ngl resources
    mpres.nglDraw           = False # False = require Ngl.draw(plot) be called (allows for other layers)
    mpres.nglFrame          = False # False = require Ngl.frame(wks) be called (allows for other layers)
    mpres.nglMaximize       = True  # True = maximize plot area in the view port
#   mpres.wkPause           = False # don't require a mouse click to clear a workstation

    ## Contour resources
    mpres.cnFillOn                  = True
    mpres.cnFillPalette         = 'NCL_default'
    mpres.cnLinesOn             = False             # no contour lines
    mpres.cnLineLabelsOn            = False             # no contour labels
    mpres.cnInfoLabelOn         = False             # no contour info label
    mpres.cnSpanFillPalette     = True
    mpres.cnFillMode                = 'RasterFill'      # solid fill each grid square (faster for large grids)
    mpres.cnRasterSmoothingOn   = False             # fill each grid square with its value (don't interpolate from neighboring cells)
    mpres.cnLevelSelectionMode  = 'ManualLevels'    # allow for cnMinLevelValF and cnMaxLevelValF to set the contour bounds

    ## Label bar resources
    mpres.lbLabelBarOn          = True
    mpres.lbOrientation         = 'horizontal'
    mpres.lbTitleOn             = True
    mpres.lbTitlePosition       = 'Bottom'      # place title below label bar
    mpres.lbTitleFontHeightF    = 0.015
    mpres.lbLabelFontHeightF    = 0.015

    ## Title resources
    mpres.tiMainFontHeightF     = 0.020
    mpres.tiXAxisFontHeightF    = 0.015
    mpres.tiYAxisFontHeightF    = 0.015

    ## Plot manager resources
    mpres.pmTickMarkDisplayMode = 'Always'  # turn on tickmarks

    ## Tickmark resources
    mpres.tmXBOn                    = True
    mpres.tmXTOn                    = False
    mpres.tmYLOn                    = True
    mpres.tmYROn                    = True
    mpres.tmXBLabelFontHeightF  = 0.008     # default 0.014
    mpres.tmYLLabelFontHeightF  = 0.008

    ## Map resources
    mpres.mpDataBaseVersion                 = 'Ncarg4_1'
    mpres.mpDataSetName                     = 'Earth..4'
    mpres.mpDataResolution                  = 'FinestResolution'
    mpres.mpOutlineBoundarySets         = 'National'
    mpres.mpOutlineSpecifiers               = ['United States:States', 'Canada:Provinces', 'Mexico:States']
    mpres.mpGridAndLimbOn                   = True
    mpres.mpUSStateLineThicknessF           = 1.0
    mpres.mpNationalLineThicknessF      = 2.0
    mpres.mpGeophysicalLineThicknessF   = 2.0
    mpres.mpPerimOn                         = True
    mpres.mpGridLineDashPattern         = 2 # lat/lon lines as dashed

    if projection == 'Cylindrical':
        mpres.mpProjection  = 'CylindricalEquidistant'
        mpres.mpLimitMode       = 'LatLon'
        mpres.mpMinLatF     = np.min(lat)
        mpres.mpMaxLatF     = np.max(lat)
        mpres.mpMinLonF     = np.min(lon)
        mpres.mpMaxLonF     = np.max(lon)

        mpres.tfDoNDCOverlay    = False
        mpres.sfXArray          = lon
        mpres.sfYArray          = lat

    elif projection == 'Lambert Conformal':
        mpres.mpProjection          = 'LambertConformal'
        mpres.mpLambertParallel1F   = true_lat1
        mpres.mpLambertParallel2F   = true_lat2
        mpres.mpLambertMeridianF    = stand_lon
        mpres.mpCenterLatF          = center_lat
        mpres.mpCenterLonF          = center_lon

        mpres.tfDoNDCOverlay        = True

        mpres.mpLimitMode           = 'Corners'
        mpres.mpLeftCornerLatF  = lat[0][0]
        mpres.mpLeftCornerLonF  = lon[0][0]
        mpres.mpRightCornerLatF = lat[n_lat-1][n_lon-1]
        mpres.mpRightCornerLonF = lon[n_lat-1][n_lon-1]

    elif projection == 'Polar':
        print('WARNING: Polar resources mpLeftAngleF, mpRightAngleF, mpTopAngleF, & mpBottomAngleF may need adjustment in functions/map_functions.py.')
        mpres.mpProjection  = 'Stereographic'
        mpres.mpCenterLatF  = center_lat
        mpres.mpCenterLonF  = center_lon

        mpres.tfDoNDCOverlay    = False
        mpres.sfXArray          = lon
        mpres.sfYArray          = lat

        mpres.mpLimitMode       = 'Angles'
        ## The following settings are found by trial and error
            ## Settings for NREL-Alaska project:
        mpres.mpLeftAngleF  = 12.5
        mpres.mpRightAngleF = 12.5
        mpres.mpTopAngleF       = 11.5
        mpres.mpBottomAngleF    = 11.0
            ## Settings for DOE MBL-SE project:
#       mpres.mpLeftAngleF  = 20.0
#       mpres.mpRightAngleF = 20.0
#       mpres.mpTopAngleF       = 19.0
#       mpres.mpBottomAngleF    = 18.0

    elif projection == 'Mercator':
        print('WARNING: Default plotting resources for Mercator projection are currently untested in functions/map_functions.py')
        mpres.mpProjection  = 'Mercator'
        mpres.mpCenterLatF  = 0.0       # This was set by wrf_map_resources in NCL. Are all Mercator-projection wrfout files set to this?? Don't know.
        mpres.mpCenterLonF  = 90.0  # This was set by wrf_map_resources in NCL. Are all Mercator-projection wrfout files set to this?? Don't know.

    else:
        print('WARNING: May need to set some base mapping resources to handle projection '+projection+' in functions/map_functions.py.')


    return mpres

def set_ngl_vcres():
    '''
   -- Sets several basic resources for adding wind barbs on top of a contour map plot.
   -- Added by JAL on 12 May 2019
   -- Inputs: none
   -- Output: Ngl.Resources object.
    '''

    vcres   = Ngl.Resources()

    ## Wind barb overlays
    vcres.nglDraw                       = False         # turn this off for overlaying vector plots
    vcres.nglFrame                      = False         # turn this off for overlaying vector plots
    vcres.tfDoNDCOverlay                = True          # needed to correctly overlay wind vectors on a map plot
    vcres.vcGlyphStyle              = 'WindBarb'    # turn on wind barbs
    vcres.vcWindBarbTickLengthF = 0.4               # default 0.3 in NCL
    vcres.vcWindBarbScaleFactorF    = 1.9438449     # convert [m/s] to [kts] for wind barb plots
    vcres.vcMinDistanceF                = 0.020         # thin vectors (NDC distance between points with wind barbs)
    vcres.vcRefMagnitudeF           = 5.0               # define vector reference magnitude
    vcres.vcRefLengthF              = 0.025         # define length of reference vector
    vcres.vcRefAnnoOrthogonalPosF   = 0.1               # move reference vector down

    return vcres

def set_ngl_polyres():
    '''
   -- Sets several basic resources for adding polymarker(s) at location(s) of interest on a map plot.
   -- Added by JAL on 12 May 2019
   -- Inputs: none
   -- Output: Ngl.Resources object.
    '''

    polyres = Ngl.Resources()

    ## Add a polymarker at the location(s) of the point(s) of interest
    polyres.gsMarkerIndex   = 12            # ploymarker style (12 = star, default)
    polyres.gsMarkerSizeF   = 15.           # polymarker size
    polyres.gsMarkerColor   = 'black'   # polymarker color

    return polyres

def set_ngl_txres():
    '''
   -- Sets a basic resource for adding text at location(s) of interest on a map plot.
   -- Added by JAL on 12 May 2019
   -- Inputs: none
   -- Output: Ngl.Resources object.
    '''

    txres   = Ngl.Resources()

    ## Add text at the location(s) of the point(s) of interest
    txres.txFontHeightF = 0.01  # font height in NDC units (default 0.05 in NCL)

    return txres

## Make a meteogram plot
class Meteogram:
    '''
    Plot a time series of meteorological data from a particular station as a meteogram with standard
    variables to visualize, including thermodynamic, kinematic, and pressure. The functions below
    control the plotting of each variable.
    TO DO: Make the subplot creation dynamic so the number of rows is not static as it is currently.
    Source: https://unidata.github.io/MetPy/latest/examples/meteogram_metpy.html (6 Aug 2021)

    class Meteogram functions (added by JAL on 16 Aug 2021):
   -- plot_ws_wd
      - wind speed and wind direction
   -- plot_t_td
      - temperature and dewpoint
   -- plot_t_td_rh
      - temperature, dewpoint, and relative humidity
   -- plot_rh
      - relative humidity
   -- plot_p
      - mean sea level pressure
   -- plot_ghi_dni
      - global horizontal irradiance and direct normal irradiance
   -- plot_ghi
      - global horizontal irradiance
   -- plot_rain
      - accumulated rainfall
   -- plot_hfx_qfx
      - sensible heat flux and moisture flux
   -- plot_th_e
      - equivalent potential temperature
   -- plot_lcl_lfc_pblh (added 19 Nov 2021)
      - lifted condensation level, level of free convection, atmospheric boundary layer depth
    '''

    def __init__(self, fig, dates, title, time=None, axis=0, n_rows=3, fontsize=24):
        '''
        Required input:
            fig: figure object
            dates: array of dates corresponding to the data
            title: overall plot title
        Optional input:
            time: Time the data is to be plotted
            axis: number that controls the new axis to be plotted (FOR FUTURE)
            n_rows: number of rows of meteogram plots
        '''
        if not time:
            time = dt.datetime.utcnow()
        self.start = dates[0]
        self.fig = fig
        self.end = dates[-1]
        self.axis_num = 0
        self.dates = mpl.dates.date2num(dates)
        self.time = time.strftime('%Y-%m-%d %H:%M UTC')
        self.fontsize = fontsize
        mpl.rcParams['font.size'] = fontsize
        mpl.rcParams['grid.color'] = 'lightgray'
        mpl.rcParams['grid.linestyle'] = '-'
        mpl.rcParams['grid.linewidth'] = 0.5
        mpl.rcParams['lines.linestyle'] = '-'
        mpl.rcParams['lines.linewidth'] = 2.0
        mpl.rcParams['lines.marker'] = 'o'
        mpl.rcParams['lines.markersize'] = 12.0
        mpl.rcParams['legend.fontsize'] = fontsize
        mpl.rcParams['legend.loc'] = 'upper center'
        mpl.rcParams['legend.framealpha'] = 1.0
        mpl.rcParams['savefig.bbox'] = 'tight'
        mpl.rcParams['axes.titlesize'] = fontsize+4
        mpl.rcParams['axes.labelsize'] = fontsize+4
        mpl.rcParams['figure.titlesize'] = fontsize+6

        ## suptitle_y seems to have no effect, and I'm not sure why...
        if n_rows == 1:
            suptitle_y = 0.75
        else:
            suptitle_y = 0.98
        self.fig.suptitle(title, y=suptitle_y)

    ## NOTE: MetPy example has wsmax as a required argument for plotting wind gusts. It's omitted here.
    def plot_ws_wd(self, ws, wd, plot_range=None, n_rows=3, n_cols=1, ind=1):
        '''
        Required input:
            ws: Wind speed (m/s)
            wd: Wind direction (degrees)
        Optional input:
            plot_range: Data range for making figure (list of (min,max,step))
            n_rows: Number of rows of subplots
            n_cols: Number of columns of subplots
            ind: Index of the subplot
        '''
        ## Plot wind speed
        self.ax1 = self.fig.add_subplot(n_rows, n_cols, ind)
        ln_ws = self.ax1.plot(self.dates, ws, label='10-m Wind Speed', color='tab:blue', marker='o', linestyle='-')
#       self.ax1.fill_between(self.dates, ws, 0)
        self.ax1.set_xlim(self.start, self.end)
        self.ax1.set_xlabel('Valid Day/Time [UTC]')
        self.ax1.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        if not plot_range:
            plot_range = [0, 20, 1]
        self.ax1.set_ylabel('Wind Speed [m/s]')
        self.ax1.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        self.ax1.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax1.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)

        ## Plot wind direction
        ax1r = self.ax1.twinx()
        ln_wd = ax1r.plot(self.dates, wd, label='10-m Wind Direction', color='tab:orange', marker='^', linestyle='-')
        ax1r.set_ylabel('Wind Direction ['+u'\N{DEGREE SIGN}]')
        ax1r.set_ylim(0, 360)
#       ax1r.set_yticks(np.arange(45, 405, 90))
#       ax1r.set_yticklabels(['NE', 'SE', 'SW', 'NW'])
        ax1r.set_yticks(np.arange(0, 361, 45))
        ax1r.set_yticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N'])

        lines = ln_ws + ln_wd
        labs = [line.get_label() for line in lines]

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.3
        self.ax1.legend(lines, labs, bbox_to_anchor=(0.5, bbox_y), ncol=2)#, prop={'size':12})

    def plot_t_td(self, t, td, plot_range=None, n_rows=3, n_cols=1, ind=2):
        '''
        Required input:
            t: Temperature (deg C)
            td: Dewpoint (deg C)
        Optional input:
            plot_range: Data range for making figure (list of (min,max,step))
            n_rows: Number of rows of subplots
            n_cols: Number of columns of subplots
            ind: Index of the subplot
        '''
        ## Plot temperature
        if not plot_range:
            plot_range = [-10, 40, 1]
        self.ax2 = self.fig.add_subplot(n_rows, n_cols, ind)#, sharex=self.ax1)
        self.ax2.set_xlim(self.start, self.end)
        ln_t = self.ax2.plot(self.dates, t, label='2-m Temperature', color='tab:red', marker='o', linestyle='-')
#       self.ax2.fill_between(self.dates, t, td, color='red')
        self.ax2.set_ylabel('Temperature ['+u'\N{DEGREE SIGN}C]')
        self.ax2.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax2.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)
        self.ax2.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        self.ax2.set_xlabel('Valid Day/Time [UTC]')
        self.ax2.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))

        ## Plot dewpoint
        ln_td = self.ax2.plot(self.dates, td, label='2-m Dewpoint', color='tab:green', marker='v', linestyle='-')
#       self.ax2.fill_between(self.dates, td, self.ax2.get_ylim()[0], color='green')
        ax2r = self.ax2.twinx()
        ax2r.set_ylim(plot_range[0], plot_range[1], plot_range[2])

        lines = ln_t + ln_td
        labs = [line.get_label() for line in lines]

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:   
            bbox_y = 1.3
        self.ax2.legend(lines, labs, bbox_to_anchor=(0.5, bbox_y), ncol=2)#, prop={'size':12})

    def plot_t_td_rh(self, t, td, rh, plot_range=None, n_rows=3, n_cols=1, ind=2):
        '''
        Required input:
            t: Temperature (deg C)
            td: Dewpoint (deg C)
            rh: Relative Humidity (%)
        Optional input:
            plot_range, n_rows, n_cols, ind: As above
        '''
        ## Plot temperature
        if not plot_range:
            plot_range = [-10, 40, 1]
        self.ax2 = self.fig.add_subplot(n_rows, n_cols, ind)
        self.ax2.set_xlim(self.start, self.end)
        self.ax2.set_xlabel('Valid Day/Time [UTC]')
        self.ax2.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        ln_t = self.ax2.plot(self.dates, t, label='2-m Temperature', color='tab:red', marker='o', linestyle='-')
        self.ax2.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        self.ax2.set_ylabel('Temperature ['+u'\N{DEGREE SIGN}C]')
        self.ax2.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax2.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)

        ## Plot dewpoint
        ln_td = self.ax2.plot(self.dates, td, label='2-m Dewpoint', color='tab:green', marker='v', linestyle='-')

        ## Plot relative humidity
        ax2r = self.ax2.twinx()
        ln_rh = ax2r.plot(self.dates, rh, label='2-m Relative Humidity', color='tab:purple', marker='s', linestyle='-')
        ax2r.set_ylim(0, 100)
        ax2r.set_ylabel('Relative Humidity [%]')

        lines = ln_t + ln_td + ln_rh
        labs = [line.get_label() for line in lines]

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.3
        self.ax2.legend(lines, labs, bbox_to_anchor=(0.5, bbox_y), ncol=3)

    def plot_rh(self, rh, plot_range=None, n_rows=3, n_cols=1, ind=3):
        '''
        Required input:
            rh: Relative humidity (%)
        Optional input:
            plot_range, n_rows, n_cols, ind: As above
        '''
        ## Plot relative humidity
        ## NOTE: It may be worth sticking this onto the temperature/dewpoint plot instead...
        if not plot_range:
            plot_range = [0, 100, 4]
        self.ax3 = self.fig.add_subplot(n_rows, n_cols, ind)#, sharex=self.ax1)
        self.ax3.set_xlim(self.start, self.end)
        self.ax3.plot(self.dates, rh, label='2-m Relative Humidity', color='tab:green', marker='s', linestyle='-')
        self.ax3.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax3.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)
        self.ax3.set_ylim(plot_range[0], plot_range[1], plot_range[2])
#       self.ax3.fill_between(self.dates, rh, self.ax3.get_ylim()[0], color='green')
        self.ax3.set_ylabel('Relative Humidity [%]')
        self.ax3.set_xlabel('Valid Day/Time [UTC]')
        self.ax3.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        ax3r = self.ax3.twinx()
        ax3r.set_ylim(plot_range[0], plot_range[1], plot_range[2])

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.3
        self.ax3.legend(bbox_to_anchor=(0.5, bbox_y))#, prop={'size':12})

    def plot_p(self, p, plot_range=None, n_rows=3, n_cols=1, ind=3):
        '''
        Required input:
            p: Mean sea level pressure (hPa)
        Optional input:
            plot_range, n_rows, n_cols, ind: As above
        '''
        ## Plot pressure
        if not plot_range:
            plot_range = [970, 1030, 2]
        self.ax4 = self.fig.add_subplot(n_rows, n_cols, ind)#, sharex=self.ax1)
        self.ax4.set_xlim(self.start, self.end)
        self.ax4.set_xlabel('Valid Day/Time [UTC]')
        self.ax4.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        self.ax4.plot(self.dates, p, 'm', label='Mean Sea Level Pressure')
        self.ax4.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        self.ax4.set_ylabel('Mean Sea Level Pressure [hPa]')
        self.ax4.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax4.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)

        ax4r = self.ax4.twinx()
        ax4r.set_ylim(plot_range[0], plot_range[1], plot_range[2])

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.3
        self.ax4.legend(bbox_to_anchor=(0.5, bbox_y))#, prop={'size': 12})

    def plot_ghi_dni(self, ghi, dni:None, plot_range=None, n_rows=3, n_cols=1, ind=3):
        '''
        Required input:
            ghi: Global horizontal irradiance (W/m2)
            dni: Direct normal irradiance (W/m2)
        Optional input:
            plot_range, n_rows, n_cols, ind: As above
        '''
        ## Plot GHI
        if not plot_range:
            plot_range = [0, 1100, 25]
        self.ax5 = self.fig.add_subplot(n_rows, n_cols, ind)
        self.ax5.set_xlim(self.start, self.end)
        self.ax5.set_xlabel('Valid Day/Time [UTC]')
        self.ax5.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        ln_ghi = self.ax5.plot(self.dates, ghi, label='Global Horizontal Irradiance', color='tab:olive', marker='o', linestyle='-')
        self.ax5.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        self.ax5.set_ylabel('Irradiance [W m$^{-2}$]')
        self.ax5.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax5.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)

        ## Plot DNI
        ln_dni = self.ax5.plot(self.dates, dni, label='Direct Horizontal Irradiance', color='tab:cyan', marker='v', linestyle='-')

        ax5r = self.ax5.twinx()
        ax5r.set_ylim(plot_range[0], plot_range[1], plot_range[2])

        lines = ln_ghi + ln_dni
        labs = [line.get_label() for line in lines]

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.3
        self.ax5.legend(lines, labs, bbox_to_anchor=(0.5, bbox_y))

    def plot_ghi(self, ghi, plot_range=None, n_rows=3, n_cols=1, ind=3):
        '''
        Required input:
            ghi: Global horizontal irradiance (W/m2)
        Optional input:
            plot_range, n_rows, n_cols, ind: As above
        '''
        if not plot_range:
            plot_range = [0, 1100, 25]
        self.ax5 = self.fig.add_subplot(n_rows, n_cols, ind)
        self.ax5.set_xlim(self.start, self.end)
        self.ax5.set_xlabel('Valid Day/Time [UTC]')
        self.ax5.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        self.ax5.plot(self.dates, ghi, label='Global Horizontal Irradiance', color='tab:olive', marker='o', linestyle='-')
        self.ax5.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        self.ax5.set_ylabel('GHI [W m$^{-2}$]')
        self.ax5.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax5.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.3
        self.ax5.legend(bbox_to_anchor=(0.5, bbox_y))

    def plot_rain(self, rain, plot_range=None, n_rows=3, n_cols=1, ind=3):
        '''
        Required input:
            rain: Total accumulated rainfall (mm)
        Optional input:
            plot_range, n_rows, n_cols, ind: As above
        '''
        ## Plot accumulated rain
        if not plot_range:
            plot_range = [0, 50, 1]
        self.ax6 = self.fig.add_subplot(n_rows, n_cols, ind)
        self.ax6.set_xlim(self.start, self.end)
        self.ax6.set_xlabel('Valid Day/Time [UTC]')
        self.ax6.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        self.ax6.plot(self.dates, rain, label='Accumulated Precipitation', color='tab:green', marker='o', linestyle='-')
        self.ax6.set_ylim(plot_range[0], plot_range[1], plot_range[2])
        self.ax6.set_ylabel('Accum. Precip. [mm]')
        self.ax6.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax6.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.3
        self.ax6.legend(bbox_to_anchor=(0.5, bbox_y))

    def plot_hfx_qfx(self, hfx, qfx, plot_range_hfx=None, plot_range_qfx=None, n_rows=3, n_cols=1, ind=3):
        '''
        Required input:
            hfx: Upward surface sensible heat flux (W/m2)
            qfx: Upward surface moisture flux (kg m-2 s-1)
        Optional input:
            plot_range_hfx, plot_range_qfx: Plot ranges for hfx and qfx, respectively (list of (min,max,step))
            n_rows, n_cols, ind: As above
        '''
        if not plot_range_hfx:
            plot_range_hfx = [-200, 800, 25]
        if not plot_range_qfx:
            plot_range_qfx = [0.0000, 0.0004, 0.00001]
        ## Plot sensible heat flux
        self.ax7 = self.fig.add_subplot(n_rows, n_cols, ind)
        self.ax7.set_xlim(self.start, self.end)
        self.ax7.set_xlabel('Valid Day/Time [UTC]')
        self.ax7.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        ln_hfx = self.ax7.plot(self.dates, hfx, label='Surface Sensible Heat Flux', color='tab:orange', marker='^', linestyle='-')
        self.ax7.set_ylim(plot_range_hfx[0], plot_range_hfx[1], plot_range_hfx[2])
        self.ax7.set_ylabel('Sens. Heat Flux [W m$^{-2}$]')
        self.ax7.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax7.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)

        ## Plot moisture flux
        ax7r = self.ax7.twinx()
        ln_qfx = ax7r.plot(self.dates, qfx, label='Surface Moisture Flux', color='tab:purple', marker='o', linestyle='-')
        ax7r.set_ylim(plot_range_qfx[0], plot_range_qfx[1], plot_range_qfx[2])
        ax7r.set_ylabel('Moist. Flux [kg m$^{-2}$ s$^{-1}$]')

        lines = ln_hfx + ln_qfx
        labs = [line.get_label() for line in lines]

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.3
        self.ax7.legend(lines, labs, bbox_to_anchor=(0.5, bbox_y), ncol=2)

    def plot_th_e(self, th_e, plot_range=None, n_rows=3, n_cols=1, ind=3, k=0):
        '''
        Required input:
            th_e: Equivalent potential temperature (degC)
        Optional input:
            plot_range, n_rows, n_cols, ind: As above
            k: k-level of theta-e (assume 0th [bottom] level)
        '''
        if not plot_range:
            plot_range = [50, 120, 1]
        ## Plot equivalent potential temperature
        self.ax8 = self.fig.add_subplot(n_rows, n_cols, ind)
        self.ax8.set_xlim(self.start, self.end)
        self.ax8.set_xlabel('Valid Day/Time [UTC]')
        self.ax8.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        self.ax8.plot(self.dates, th_e, label='Equivalent Potential Temperature (k='+str(k)+')', color='tab:pink', marker='o', linestyle='-')
        self.ax8.set_ylim(plot_range[0], plot_range[1], plot_range[2])
#       self.ax8.set_ylabel('Equiv. Pot. Temp. ['+u'\N{DEGREE SIGN}C]')
        self.ax8.set_ylabel('Equiv. Pot. Temp. [K]')
        self.ax8.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax8.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.2)

        ax8r = self.ax8.twinx()
        ax8r.set_ylim(plot_range[0], plot_range[1], plot_range[2])

        if n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.35
        self.ax8.legend(bbox_to_anchor=(0.5, bbox_y))

    def plot_lcl_lfc_pblh(self, lcl, lfc, pblh, plot_range_lcl_lfc=None, plot_range_pblh=None, n_rows=3, n_cols=1, ind=3):
        '''
        Required input:
            lcl: Lifted condensation level (m)
            lfc: Level of free convection (m)
            pblh: Atmospheric boundary layer depth (m)
        Optional input:
            plot_range, n_rows, n_cols, ind: As above
        '''
        if not plot_range_lcl_lfc:
            plot_range_lcl_lfc = [0, 5000, 50]
        if not plot_range_pblh:
            plot_range_pblh = [0, 3000, 50]

        ## Plot LCL
        self.ax9 = self.fig.add_subplot(n_rows, n_cols, ind)
        self.ax9.set_xlim(self.start, self.end)
        self.ax9.set_xlabel('Valid Day/Time [UTC]')
        self.ax9.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d/%H%M'))
        ln_lcl = self.ax9.plot(self.dates, lcl, label='Lifted Condensation Level', color='tab:green', marker='o', linestyle='-')
        self.ax9.set_ylim(plot_range_lcl_lfc[0], plot_range_lcl_lfc[1], plot_range_lcl_lfc[2])
        self.ax9.set_ylabel('LCL, LFC [m]')
        self.ax9.grid(b=True, which='major', axis='y', color='black', linestyle='--', linewidth=0.5)
        self.ax9.grid(b=True, which='major', axis='x', color='black', linestyle='--', linewidth=0.5)

        ## Plot LFC
        ln_lfc = self.ax9.plot(self.dates, lfc, label='Level of Free Convection', color='tab:red', marker='v', linestyle='-')

        ## Plot PBL depth
        ax9r = self.ax9.twinx()
        ln_pblh = ax9r.plot(self.dates, pblh, label='Atmospheric Boundary Layer Depth', color='tab:blue', marker='s', linestyle='-')
        ax9r.set_ylim(plot_range_pblh[0], plot_range_pblh[1], plot_range_pblh[2])
        ax9r.set_ylabel('ABL Depth [m]')

        lines = ln_lcl + ln_lfc + ln_pblh
        labs = [line.get_label() for line in lines]

        if n_rows == 1:
            bbox_y = 1.05
        elif n_rows == 2:
            bbox_y = 1.1
        elif n_rows == 3:
            bbox_y = 1.15
        elif n_rows == 4:
            bbox_y = 1.2
        elif n_rows == 5:
            bbox_y = 1.27
        else:
            bbox_y = 1.35

        if n_rows == 1:
#           n_col = 1
#           bbox_y = 1.15
            n_col = 2
            bbox_y = 1.105
        else:
            n_col = 3
        self.ax9.legend(lines, labs, bbox_to_anchor=(0.5, bbox_y), ncol=n_col)
