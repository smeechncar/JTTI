'''
stats_by_epa_region.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 19 Oct 2023

This script calculates and plots bulk statistics within each EPA region.
The observations are from AirNow O3 and PM2.5 stations.
The forecasts are analysis time-series from raw CAMS reanalysis, raw CAMS fcsts, & raw CMAQ fcsts.
Hourly-mean stats are also calculated over the entire time period, domain-wide and for each region.
Daily-mean stats are also calculated as a time series, domain-wide and for each region.
'''

import sys
import pathlib
import argparse
import warnings
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
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
import seaborn as sns
from cycler import cycler
from pandas.plotting import register_matplotlib_converters
from functions import map_functions
from functions import gen_funcs

def plot_stat_xy(fname, data, ymin, ymax, xticklab, xlab, ylab, suptitle, title, colors, dashes, markers, lglab, fontsize=14, lgloc='lower right'):
    plt.rcParams['figure.figsize'] = (8,6)
    plt.rcParams['figure.titlesize'] = fontsize+2
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['axes.prop_cycle'] = cycler(color=colors)+cycler(marker=markers)+cycler(linestyle=dashes)
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['grid.color'] = 'lightgray'
    plt.rcParams['grid.linestyle'] = '-'
    plt.rcParams['grid.linewidth'] = 0.5
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['legend.fontsize'] = fontsize
    plt.rcParams['legend.loc'] = lgloc
    plt.rcParams['legend.framealpha'] = 1.0

    x_dim = data.shape[0]
    x = np.array(range(x_dim))
    x_pos = [i for i, _ in enumerate(x)]
    print('-- Creating plot: '+str(fname))
    fig, ax = plt.subplots()
    ax.plot(data)
    ax.set(xlabel=xlab, ylabel=ylab, title=title)
    ax.set_ylim(bottom=ymin, top=ymax)
    ax.minorticks_on()
    ax.grid(True, which='major', axis='both')
    ax.grid(True, which='minor', axis='both', linestyle='--', color='lightgray', linewidth=0.25)
#    ax.tick_params(axis='x', which='minor', bottom=False)
#    plt.xticks(ticks=x, labels=xticklab)
    plt.xticks(np.arange(min(x), max(x)+1, 4))
    plt.hlines(0.0, x[0], x[-1], color='black')
    plt.figlegend(labels=lglab, ncol=2, bbox_to_anchor=(0.88, 0.11))
#    plt.figlegend(labels=lglab, ncol=2)
    plt.suptitle(suptitle)
    plt.savefig(str(fname))
    plt.close()

def plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr,
        ymin1, ymax1, ymin2, ymax2, ymin3, ymax3,
        xlab, ylab1, ylab2, ylab3, suptitle, title1, title2, title3,
        colors, dashes, markers, lglab, fontsize=14, lgloc='lower center'):
    plt.rcParams['figure.figsize'] = (21,6)
    plt.rcParams['figure.titlesize'] = fontsize+2
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['axes.prop_cycle'] = cycler(color=colors)+cycler(marker=markers)+cycler(linestyle=dashes)
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['grid.color'] = 'lightgray'
    plt.rcParams['grid.linestyle'] = '-'
    plt.rcParams['grid.linewidth'] = 0.5
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['legend.fontsize'] = fontsize
    plt.rcParams['legend.loc'] = lgloc
    plt.rcParams['legend.framealpha'] = 1.0

#    x_dim = data1.shape[0]
#    x = np.array(range(x_dim))
#    x_pos = [i for i, _ in enumerate(x)]
    y_dim1 = data1.shape[1]
    y_dim2 = data2.shape[1]
    y_dim3 = data3.shape[1]

    print('-- Creating plot: '+str(fname))
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.suptitle(suptitle)
    for yy in range(y_dim1):
        ax1.plot(xarr, data1[:,yy])
    ax1.set(xlabel=xlab, ylabel=ylab1, title=title1)

    for yy in range(y_dim2):
        ax2.plot(xarr, data2[:,yy])
    ax2.set(xlabel=xlab, ylabel=ylab2, title=title2)

    for yy in range(y_dim3):
        ax3.plot(xarr, data3[:,yy])
    ax3.set(xlabel=xlab, ylabel=ylab3, title=title3)

    for ax in [ax1, ax2, ax3]:
        ax.grid(True, which='major', axis='both')
        ax.grid(True, which='minor', axis='both', linestyle='--', color='lightgray', linewidth=0.25)
        ax.hlines(0.0, xarr[0]-0.5, xarr[-1]+0.5, color='dimgray')
        ax.xaxis.set_major_locator(mticker.MultipleLocator(3))
        ax.xaxis.set_minor_locator(mticker.MultipleLocator(1))

#    plt.figlegend(labels=lglab, ncol=2, bbox_to_anchor=(0.80, 0.11))
    plt.figlegend(labels=lglab, ncol=6)
    plt.savefig(str(fname))
    plt.close()

def plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, pd_times,
        ymin1, ymax1, ymin2, ymax2, ymin3, ymax3,
        xlab, ylab1, ylab2, ylab3, suptitle, title1, title2, title3,
        colors, dashes, markers, lglab, fontsize=14, lgloc='lower center'):
    plt.rcParams['figure.figsize'] = (15,12)
    plt.rcParams['figure.titlesize'] = fontsize+2
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['axes.prop_cycle'] = cycler(color=colors)+cycler(marker=markers)+cycler(linestyle=dashes)
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['grid.color'] = 'lightgray'
    plt.rcParams['grid.linestyle'] = '-'
    plt.rcParams['grid.linewidth'] = 0.5
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['legend.fontsize'] = fontsize
    plt.rcParams['legend.loc'] = lgloc
    plt.rcParams['legend.framealpha'] = 1.0

    date_fmt = '%b\n%Y'

#    x_dim = data1.shape[0]
    y_dim1 = data1.shape[1]
#    x = np.array(range(x_dim))
#    x_pos = [i for i, _ in enumerate(x)]
    print('-- Creating plot: '+str(fname))
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0.12)
    fig.suptitle(suptitle, y=0.93)
    for yy in range(y_dim1):
        ax1.plot(pd_times, data1[:,yy])
    ax1.set(ylabel=ylab1, title=title1)
#    ax1.set_ylim(bottom=ymin1, top=ymax1)

    y_dim2 = data2.shape[1]
    for yy in range(y_dim2):
        ax2.plot(pd_times, data2[:,yy])
    ax2.set(ylabel=ylab2, title=title2)
#    ax2.set_ylim(bottom=ymin2, top=ymax2)

    y_dim3 = data3.shape[1]
    for yy in range(y_dim3):
        ax3.plot(pd_times, data3[:,yy])
    ax3.set(ylabel=ylab3, title=title3)
#    ax3.set_ylim(bottom=ymin3, top=ymax3)

    for ax in [ax1, ax2, ax3]:
        ax.xaxis.set_major_formatter(mdates.DateFormatter(date_fmt))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
        ax.xaxis.set_minor_locator(mdates.WeekdayLocator(byweekday=mdates.SU))
        ax.set_xlim(left=pd_times[0]-dt.timedelta(days=1), right=pd_times[-1]+dt.timedelta(days=1))
        ax.grid(True, which='major', axis='both')
        ax.grid(True, which='minor', axis='both', linestyle='--', color='lightgray', linewidth=0.25)
        ax.hlines(0.0, pd_times[0]-dt.timedelta(hours=12), pd_times[-1]+dt.timedelta(hours=12), color='dimgray')

#    plt.figlegend(labels=lglab, ncol=2, bbox_to_anchor=(0.80, 0.11))
    plt.figlegend(labels=lglab, ncol=6)
    plt.savefig(str(fname))
    plt.close()


def plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1, xlab2, xlab3, ylab1, ylab2, ylab3, suptitle, title1, title2, title3, colors, fontsize=14, lgloc='lower center', ymin_box=None, yscale_box='linear'):
    plt.rcParams['figure.figsize'] = (20,6)
    plt.rcParams['figure.titlesize'] = fontsize+2
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['axes.prop_cycle'] = cycler(color=colors)
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['grid.color'] = 'lightgray'
    plt.rcParams['grid.linestyle'] = '-'
    plt.rcParams['grid.linewidth'] = 0.5
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['legend.fontsize'] = fontsize
    plt.rcParams['legend.loc'] = lgloc
    plt.rcParams['legend.framealpha'] = 1.0

    print('-- Creating plot: '+str(fname))
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.suptitle(suptitle, y=0.98)

    ## Matplotlib boxplots cannot handle NaNs, so account for this (note: Seaborn boxplots can handle NaNs)
    mask = ~np.isnan(data1)
    data1_filtered = [d[m] for d, m in zip(data1.T, mask.T)]

    box1 = ax1.boxplot(data1_filtered, vert=True, patch_artist=True, labels=xlab1)
    ax1.set_yscale(yscale_box)
    ax1.set_ylabel(ylab1)
    ax1.set_title(title1)
    for patch, color in zip(box1['boxes'], colors):
        patch.set_facecolor(color)
    for median in box1['medians']:
        median.set_color('yellow')
    if ymin_box is not None:
        ax1.set_ylim(bottom=ymin_box)

    ax2.bar(xlab2, data2, color=colors)
    ax2.set_ylabel(ylab2)
    ax2.set_title(title2)

    ax3.bar(xlab3, data3, color=colors)
    ax3.set_ylabel(ylab3)
    ax3.set_title(title3)

    for ax in [ax1, ax2, ax3]:
        ax.axhline(y=0.0, xmin=0.02, xmax=0.98, color='black')
        ax.grid(True, which='major', axis='y')
        ax.grid(True, which='minor', axis='y', linestyle='--', color='lightgray', linewidth=0.25)
        for label in ax.get_xticklabels():
            label.set_rotation(90)

    plt.savefig(str(fname))
    plt.close()

def plot_stat_box_bar_3panel_vert(fname, data1, data2, data3, groups, xlab1, xlab2, xlab3, ylab1, ylab2, ylab3,
        suptitle, title1, title2, title3, colors, fontsize=14, lgloc='lower center', ymin_box=None, yscale_box='linear'):
    plt.rcParams['figure.figsize'] = (15,12)
    plt.rcParams['figure.titlesize'] = fontsize+2
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['axes.prop_cycle'] = cycler(color=colors)
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['grid.color'] = 'lightgray'
    plt.rcParams['grid.linestyle'] = '-'
    plt.rcParams['grid.linewidth'] = 0.5
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['legend.fontsize'] = fontsize-1
    plt.rcParams['legend.loc'] = lgloc
    plt.rcParams['legend.framealpha'] = 0.80

    print('-- Creating plot: '+str(fname))
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    fig.suptitle(suptitle, y=0.98)

    ## Matplotlib boxplots cannot handle NaNs, so account for this (note: Seaborn boxplots can handle NaNs)
#    mask = ~np.isnan(data1)
#    data1_filtered = [d[m] for d, m in zip(data1.T, mask.T)]

#    box1 = ax1.boxplot(data1_filtered, vert=True, patch_artist=True, labels=xlab1)
    
#    ax1.set_yscale(yscale_box)
#    ax1.set_ylabel(ylab1)
#    ax1.set_title(title1)
#    for patch, color in zip(box1['boxes'], colors):
#        patch.set_facecolor(color)
#    for median in box1['medians']:
#        median.set_color('yellow')
#    if ymin_box is not None:
#        ax1.set_ylim(bottom=ymin_box)

    ## Draw a grouped box plot with seaborn and the input pandas dataframe
    ncol = data1['Model'].nunique()
    sns.boxplot(data=data1, x='Region', y='Value', hue='Model', ax=ax1, medianprops={'color':'yellow'})
    sns.move_legend(ax1, 'upper center', title=None, ncol=ncol, bbox_to_anchor=(0.50, 0.95), bbox_transform=fig.transFigure)
    ax1.set_yscale(yscale_box)
    ax1.set_xlabel('')
    ax1.set_ylabel(ylab1)
    ax1.set_title(title1)
    if ymin_box is not None:
        ax1.set_ylim(bottom=ymin_box)
    ax1.tick_params('x', labelbottom=False)

    ## Draw a grouped bar chart with seaborn and the input pandas dataframe
    sns.barplot(data=data2, x='Region', y='Value', hue='Model', ax=ax2)
#    sns.move_legend(ax2, 'upper center', title=None, ncol=5)
    ax2.get_legend().remove()
    ax2.set_xlabel('')
    ax2.set_ylabel(ylab2)
    ax2.set_title(title2)
    ax2.tick_params('x', labelbottom=False)

    sns.barplot(data=data3, x='Region', y='Value', hue='Model', ax=ax3)
    ax3.get_legend().remove()
    ax3.set_xlabel('')
    ax3.set_ylabel(ylab3)
    ax3.set_title(title3)

    for ax in [ax1, ax2, ax3]:
        ax.axhline(y=0.0, xmin=0.00, xmax=1.00, color='black')
        ax.grid(True, which='major', axis='y')
        ax.grid(True, which='minor', axis='y', linestyle='--', color='lightgray', linewidth=0.25)
        for label in ax.get_xticklabels():
            label.set_rotation(30)

    plt.savefig(str(fname))
    plt.close()

def main():
    beg_date_str = '2020-08-01_00'
    end_date_str = '2021-12-31_23'
    plot_type = 'eps'

    calc_stats_all = True
    calc_stats_use = False
    calc_stats_val = False
    write_pairs_files = False
    read_pairs_files = True

    plot_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','plots','anal')
    epa_codes_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI')
    cams_ra_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4')
    cams_fc_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')
    cmaq_fc_raw_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','bcdata','cmaq_hourly','sites_static')
    cams_fc_bm3_dir = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','merge','sites_static')
    cmaq_fc_bm3_dir = cams_fc_bm3_dir

    fmt_yyyymmdd_hhmm = '%Y%m%d_%H%M'
    fmt_yyyymmdd_hh = '%Y%m%d_%H'
    fmt_yyyymmddhh  = '%Y%m%d%H'
    fmt_yyyymmdd = '%Y%m%d'
    fmt_yyyy = '%Y'
    fmt_mm = '%m'
    fmt_hh = '%H'
    fmt_date = '%Y-%m-%d'
    fmt_date_hh = '%Y-%m-%d_%H'
    fmt_date_file = fmt_yyyymmdd_hhmm
    fmt_date_plot = '%d %b %Y/%H%M'

    mpl_o3   = '$\mathregular{O_3}$'
    mpl_pm25 = '$\mathregular{PM_{2.5}}$'
    mpl_ugm3 = '\u03bcg $\mathregular{m^{-3}}$'

    en_dash = u'\u2013'
    em_dash = u'\u2014'

    ## Create datetime array of times to analyze
    dt_beg_valid = dt.datetime.strptime(beg_date_str, fmt_date_hh)
    dt_end_valid = dt.datetime.strptime(end_date_str, fmt_date_hh)
    dt_all_valid = pd.date_range(start=dt_beg_valid, end=dt_end_valid, freq='1H')
    dt_dates_valid = pd.date_range(start=dt_beg_valid, end=dt_end_valid, freq='1D')
    n_valid_times = len(dt_all_valid)
    n_valid_days = (dt_end_valid - dt_beg_valid).days + 1
    date_range_file = dt_dates_valid[0].strftime(fmt_yyyymmdd)+'-'+dt_dates_valid[-1].strftime(fmt_yyyymmdd)

    ## Read in the EPA codes files that Ju-Hye made (PM2.5 obs are misordered & missorted)
#    epa_codes_file_use = epa_codes_dir.joinpath('EPA_CODE_O3_PM25_USE.nc')
#    epa_codes_file_val = epa_codes_dir.joinpath('EPA_CODE_O3_PM25_VAL.nc')
    ## Instead read in the corrected EPA codes files that I made
    epa_codes_file_use = epa_codes_dir.joinpath('epa_code_o3_pm25_use_jared.nc')
    epa_codes_file_val = epa_codes_dir.joinpath('epa_code_o3_pm25_val_jared.nc')

    print('Reading '+str(epa_codes_file_use))
    epa_use_ds = xr.open_dataset(epa_codes_file_use)
#    obs_o3_reg_use = epa_use_ds.epa_code_o3.values.astype(int)
    obs_o3_reg_use = epa_use_ds.obs_o3_reg_use.values
    obs_o3_lat_use = epa_use_ds.obs_o3_lat_use.values
    obs_o3_lon_use = epa_use_ds.obs_o3_lon_use.values
    obs_o3_sid_use = epa_use_ds.obs_o3_sid_use.values
#    obs_pm25_reg_use = epa_use_ds.epa_code_pm25.values.astype(int)
    obs_pm25_reg_use = epa_use_ds.obs_pm25_reg_use.values
    obs_pm25_lat_use = epa_use_ds.obs_pm25_lat_use.values
    obs_pm25_lon_use = epa_use_ds.obs_pm25_lon_use.values
    obs_pm25_sid_use = epa_use_ds.obs_pm25_sid_use.values
    n_obs_o3_use   = len(obs_o3_reg_use)
    n_obs_pm25_use = len(obs_pm25_reg_use)

    print('Reading '+str(epa_codes_file_val))
    epa_val_ds = xr.open_dataset(epa_codes_file_val)
#    obs_o3_reg_val = epa_val_ds.epa_code_o3.values.astype(int)
    obs_o3_reg_val = epa_val_ds.obs_o3_reg_val.values
    obs_o3_lat_val = epa_val_ds.obs_o3_lat_val.values
    obs_o3_lon_val = epa_val_ds.obs_o3_lon_val.values
    obs_o3_sid_val = epa_val_ds.obs_o3_sid_val.values
#    obs_pm25_reg_val = epa_val_ds.epa_code_pm25.values.astype(int)
    obs_pm25_reg_val = epa_val_ds.obs_pm25_reg_val.values
    obs_pm25_lat_val = epa_val_ds.obs_pm25_lat_val.values
    obs_pm25_lon_val = epa_val_ds.obs_pm25_lon_val.values
    obs_pm25_sid_val = epa_val_ds.obs_pm25_sid_val.values
    n_obs_o3_val   = len(obs_o3_reg_val)
    n_obs_pm25_val = len(obs_pm25_reg_val)

    ## Combine the arrays
#    obs_o3_reg_all   = np.append(epa_code_o3_use, epa_code_o3_val)
    obs_o3_reg_all    = np.append(obs_o3_reg_use, obs_o3_reg_val)
    obs_o3_lat_all    = np.append(obs_o3_lat_use, obs_o3_lat_val)
    obs_o3_lon_all    = np.append(obs_o3_lon_use, obs_o3_lon_val)
#    obs_pm25_reg_all = np.append(epa_code_pm25_use, epa_code_pm25_val)
    obs_pm25_reg_all  = np.append(obs_pm25_reg_use, obs_pm25_reg_val)
    obs_pm25_lat_all  = np.append(obs_pm25_lat_use, obs_pm25_lat_val)
    obs_pm25_lon_all  = np.append(obs_pm25_lon_use, obs_pm25_lon_val)
    n_obs_o3_all   = len(obs_o3_reg_all)
    n_obs_pm25_all = len(obs_pm25_reg_all)

#    print(n_obs_o3_use)
#    print(n_obs_o3_val)
#    print(n_obs_o3_all)
#    print(n_obs_pm25_use)
#    print(n_obs_pm25_val)
#    print(n_obs_pm25_all)

    ## Find the indices corresponding to each EPA region
    inds_o3_reg01_all = np.where(obs_o3_reg_all == 1)[0]
    inds_o3_reg02_all = np.where(obs_o3_reg_all == 2)[0]
    inds_o3_reg03_all = np.where(obs_o3_reg_all == 3)[0]
    inds_o3_reg04_all = np.where(obs_o3_reg_all == 4)[0]
    inds_o3_reg05_all = np.where(obs_o3_reg_all == 5)[0]
    inds_o3_reg06_all = np.where(obs_o3_reg_all == 6)[0]
    inds_o3_reg07_all = np.where(obs_o3_reg_all == 7)[0]
    inds_o3_reg08_all = np.where(obs_o3_reg_all == 8)[0]
    inds_o3_reg09_all = np.where(obs_o3_reg_all == 9)[0]
    inds_o3_reg10_all = np.where(obs_o3_reg_all ==10)[0]
    inds_o3_reg99_all = np.where(obs_o3_reg_all ==99)[0]
    inds_pm25_reg01_all = np.where(obs_pm25_reg_all == 1)[0]
    inds_pm25_reg02_all = np.where(obs_pm25_reg_all == 2)[0]
    inds_pm25_reg03_all = np.where(obs_pm25_reg_all == 3)[0]
    inds_pm25_reg04_all = np.where(obs_pm25_reg_all == 4)[0]
    inds_pm25_reg05_all = np.where(obs_pm25_reg_all == 5)[0]
    inds_pm25_reg06_all = np.where(obs_pm25_reg_all == 6)[0]
    inds_pm25_reg07_all = np.where(obs_pm25_reg_all == 7)[0]
    inds_pm25_reg08_all = np.where(obs_pm25_reg_all == 8)[0]
    inds_pm25_reg09_all = np.where(obs_pm25_reg_all == 9)[0]
    inds_pm25_reg10_all = np.where(obs_pm25_reg_all ==10)[0]
    inds_pm25_reg99_all = np.where(obs_pm25_reg_all ==99)[0]

    inds_o3_reg01_val = np.where(obs_o3_reg_val == 1)[0]
    inds_o3_reg02_val = np.where(obs_o3_reg_val == 2)[0]
    inds_o3_reg03_val = np.where(obs_o3_reg_val == 3)[0]
    inds_o3_reg04_val = np.where(obs_o3_reg_val == 4)[0]
    inds_o3_reg05_val = np.where(obs_o3_reg_val == 5)[0]
    inds_o3_reg06_val = np.where(obs_o3_reg_val == 6)[0]
    inds_o3_reg07_val = np.where(obs_o3_reg_val == 7)[0]
    inds_o3_reg08_val = np.where(obs_o3_reg_val == 8)[0]
    inds_o3_reg09_val = np.where(obs_o3_reg_val == 9)[0]
    inds_o3_reg10_val = np.where(obs_o3_reg_val ==10)[0]
    inds_o3_reg99_val = np.where(obs_o3_reg_val ==99)[0]
    inds_pm25_reg01_val = np.where(obs_pm25_reg_val == 1)[0]
    inds_pm25_reg02_val = np.where(obs_pm25_reg_val == 2)[0]
    inds_pm25_reg03_val = np.where(obs_pm25_reg_val == 3)[0]
    inds_pm25_reg04_val = np.where(obs_pm25_reg_val == 4)[0]
    inds_pm25_reg05_val = np.where(obs_pm25_reg_val == 5)[0]
    inds_pm25_reg06_val = np.where(obs_pm25_reg_val == 6)[0]
    inds_pm25_reg07_val = np.where(obs_pm25_reg_val == 7)[0]
    inds_pm25_reg08_val = np.where(obs_pm25_reg_val == 8)[0]
    inds_pm25_reg09_val = np.where(obs_pm25_reg_val == 9)[0]
    inds_pm25_reg10_val = np.where(obs_pm25_reg_val ==10)[0]
    inds_pm25_reg99_val = np.where(obs_pm25_reg_val ==99)[0]

    ## If pairs files are to be written, then all the individual obs/forecast files need to be read in
    if write_pairs_files:
        ## Define arrays in advance rather than appending each time through the time loop
        cams_ra_raw_o3_con_use   = np.full([n_valid_times, n_obs_o3_use], np.nan)
        cams_ra_raw_o3_con_val   = np.full([n_valid_times, n_obs_o3_val], np.nan)
        cams_ra_raw_pm25_con_use = np.full([n_valid_times, n_obs_pm25_use], np.nan)
        cams_ra_raw_pm25_con_val = np.full([n_valid_times, n_obs_pm25_val], np.nan)

        cams_fc_raw_o3_con_use   = np.full([n_valid_times, n_obs_o3_use], np.nan)
        cams_fc_raw_o3_con_val   = np.full([n_valid_times, n_obs_o3_val], np.nan)
        cams_fc_raw_pm25_con_use = np.full([n_valid_times, n_obs_pm25_use], np.nan)
        cams_fc_raw_pm25_con_val = np.full([n_valid_times, n_obs_pm25_val], np.nan)

        cmaq_fc_raw_o3_con_use   = np.full([n_valid_times, n_obs_o3_use], np.nan)
        cmaq_fc_raw_o3_con_val   = np.full([n_valid_times, n_obs_o3_val], np.nan)
        cmaq_fc_raw_pm25_con_use = np.full([n_valid_times, n_obs_pm25_use], np.nan)
        cmaq_fc_raw_pm25_con_val = np.full([n_valid_times, n_obs_pm25_val], np.nan)

        cams_fc_bm3_o3_con_use   = np.full([n_valid_times, n_obs_o3_use], np.nan)
        cams_fc_bm3_o3_con_val   = np.full([n_valid_times, n_obs_o3_val], np.nan)
        cams_fc_bm3_pm25_con_use = np.full([n_valid_times, n_obs_pm25_use], np.nan)
        cams_fc_bm3_pm25_con_val = np.full([n_valid_times, n_obs_pm25_val], np.nan)

        cmaq_fc_bm3_o3_con_use   = np.full([n_valid_times, n_obs_o3_use], np.nan)
        cmaq_fc_bm3_o3_con_val   = np.full([n_valid_times, n_obs_o3_val], np.nan)
        cmaq_fc_bm3_pm25_con_use = np.full([n_valid_times, n_obs_pm25_use], np.nan)
        cmaq_fc_bm3_pm25_con_val = np.full([n_valid_times, n_obs_pm25_val], np.nan)

        obs_o3_con_use   = np.full([n_valid_times, n_obs_o3_use], np.nan)
        obs_o3_con_val   = np.full([n_valid_times, n_obs_o3_val], np.nan)
        obs_pm25_con_use = np.full([n_valid_times, n_obs_pm25_use], np.nan)
        obs_pm25_con_val = np.full([n_valid_times, n_obs_pm25_val], np.nan)
#        print(obs_o3_con_use.shape)
#        print(obs_pm25_con_use.shape)
#        print(obs_o3_con_val.shape)
#        print(obs_pm25_con_val.shape)

        ## Loop over times
        for tt in range(n_valid_times):
            this_dt = dt_all_valid[tt]
            this_yr = this_dt.strftime(fmt_yyyy)
            this_mo = this_dt.strftime(fmt_mm)
            this_hr = this_dt.strftime(fmt_hh)
            this_dt_ymd  = this_dt.strftime(fmt_yyyymmdd)
            this_dt_ymdh = this_dt.strftime(fmt_yyyymmdd_hh)
            this_dt_file = this_dt.strftime(fmt_date_file)
            this_dt_plot = this_dt.strftime(fmt_date_plot)

            print('Reading model data for '+this_dt_plot)

            ## Read in raw CAMS reanalysis
            this_cams_ra_raw_dir = cams_ra_raw_dir.joinpath(this_yr,this_mo)
            cams_ra_raw_use_file = this_cams_ra_raw_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_static_'+this_dt_file+'_use.nc')
            cams_ra_raw_val_file = this_cams_ra_raw_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_static_'+this_dt_file+'_val.nc')

            #print('Reading '+str(cams_ra_raw_use_file))
            cams_ra_raw_use_ds = xr.open_dataset(cams_ra_raw_use_file)
            cams_ra_raw_o3_con_use[tt,:]   = cams_ra_raw_use_ds.cams_o3_con_use.values
            cams_ra_raw_pm25_con_use[tt,:] = cams_ra_raw_use_ds.cams_pm25_con_use.values

            #print('Reading '+str(cams_ra_raw_val_file))
            cams_ra_raw_val_ds = xr.open_dataset(cams_ra_raw_val_file)
            cams_ra_raw_o3_con_val[tt,:]   = cams_ra_raw_val_ds.cams_o3_con_val.values
            cams_ra_raw_pm25_con_val[tt,:] = cams_ra_raw_val_ds.cams_pm25_con_val.values

            ## Read in raw CAMS forecasts
            if int(this_hr) < 12:
                this_cycle = this_dt_ymd+'_00'
            else:
                this_cycle = this_dt_ymd+'_12'
            this_cams_fc_raw_dir = cams_fc_raw_dir.joinpath(this_yr,this_mo,this_cycle)
            cams_fc_raw_use_file = this_cams_fc_raw_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_static_'+this_dt_file+'_use.nc')
            cams_fc_raw_val_file = this_cams_fc_raw_dir.joinpath('cams_regrid_cmaq_airnow_pm2.5_o3_static_'+this_dt_file+'_val.nc')

            #print('Reading '+str(cams_fc_raw_use_file))
            cams_fc_raw_use_ds = xr.open_dataset(cams_fc_raw_use_file)
            cams_fc_raw_o3_con_use[tt,:]   = cams_fc_raw_use_ds.cams_o3_con_use.values
            cams_fc_raw_pm25_con_use[tt,:] = cams_fc_raw_use_ds.cams_pm25_con_use.values

            #print('Reading '+str(cams_fc_raw_val_file))
            cams_fc_raw_val_ds = xr.open_dataset(cams_fc_raw_val_file)
            cams_fc_raw_o3_con_val[tt,:]   = cams_fc_raw_val_ds.cams_o3_con_val.values
            cams_fc_raw_pm25_con_val[tt,:] = cams_fc_raw_val_ds.cams_pm25_con_val.values

            ## Read in raw CMAQ forecasts
            this_cmaq_fc_raw_dir = cmaq_fc_raw_dir.joinpath(this_yr,this_mo)
            cmaq_fc_raw_use_file = this_cmaq_fc_raw_dir.joinpath('cmaq_airnow_pm2.5_o3_static_'+this_dt_file+'_use.nc')
            cmaq_fc_raw_val_file = this_cmaq_fc_raw_dir.joinpath('cmaq_airnow_pm2.5_o3_static_'+this_dt_file+'_val.nc')

            #print('Reading '+str(cmaq_fc_raw_use_file))
            cmaq_fc_raw_use_ds = xr.open_dataset(cmaq_fc_raw_use_file)
            cmaq_fc_raw_o3_con_use[tt,:]   = cmaq_fc_raw_use_ds.cmaq_o3_con_use.values
            cmaq_fc_raw_pm25_con_use[tt,:] = cmaq_fc_raw_use_ds.cmaq_pm25_con_use.values

            #print('Reading '+str(cmaq_fc_raw_val_file))
            cmaq_fc_raw_val_ds = xr.open_dataset(cmaq_fc_raw_val_file)
            cmaq_fc_raw_o3_con_val[tt,:]   = cmaq_fc_raw_val_ds.cmaq_o3_con_val.values
            cmaq_fc_raw_pm25_con_val[tt,:] = cmaq_fc_raw_val_ds.cmaq_pm25_con_val.values

            ## Read in BM3 CAMS & CMAQ forecasts
            ## Also read the observations from these files
            this_cmaq_fc_bm3_dir = cmaq_fc_bm3_dir.joinpath(this_yr,this_mo)
            cmaq_cams_bm3_use_file = this_cmaq_fc_bm3_dir.joinpath('cmaq_bc_cams_bc_airnow_pm2.5_o3_static_'+this_dt_file+'_use.nc')
            cmaq_cams_bm3_val_file = this_cmaq_fc_bm3_dir.joinpath('cmaq_bc_cams_bc_airnow_pm2.5_o3_static_'+this_dt_file+'_val.nc')

            #print('Reading '+str(cmaq_cams_bm3_use_file))
            cmaq_cams_fc_bm3_use_ds = xr.open_dataset(cmaq_cams_bm3_use_file)
            cmaq_fc_bm3_o3_con_use[tt,:]   = cmaq_cams_fc_bm3_use_ds.cmaq_bm3_o3_con_use.values
            cams_fc_bm3_o3_con_use[tt,:]   = cmaq_cams_fc_bm3_use_ds.cams_bm3_o3_con_use.values
            cmaq_fc_bm3_pm25_con_use[tt,:] = cmaq_cams_fc_bm3_use_ds.cmaq_bm3_pm25_con_use.values
            cams_fc_bm3_pm25_con_use[tt,:] = cmaq_cams_fc_bm3_use_ds.cams_bm3_pm25_con_use.values
            obs_o3_con_use[tt,:]   = cmaq_cams_fc_bm3_use_ds.obs_o3_con_use.values
            obs_pm25_con_use[tt,:] = cmaq_cams_fc_bm3_use_ds.obs_pm25_con_use.values

            #print('Reading '+str(cmaq_cams_bm3_val_file))
            cmaq_cams_fc_bm3_val_ds = xr.open_dataset(cmaq_cams_bm3_val_file)
            cmaq_fc_bm3_o3_con_val[tt,:]   = cmaq_cams_fc_bm3_val_ds.cmaq_bm3_o3_con_val.values
            cams_fc_bm3_o3_con_val[tt,:]   = cmaq_cams_fc_bm3_val_ds.cams_bm3_o3_con_val.values
            cmaq_fc_bm3_pm25_con_val[tt,:] = cmaq_cams_fc_bm3_val_ds.cmaq_bm3_pm25_con_val.values
            cams_fc_bm3_pm25_con_val[tt,:] = cmaq_cams_fc_bm3_val_ds.cams_bm3_pm25_con_val.values
            obs_o3_con_val[tt,:]   = cmaq_cams_fc_bm3_val_ds.obs_o3_con_val.values
            obs_pm25_con_val[tt,:] = cmaq_cams_fc_bm3_val_ds.obs_pm25_con_val.values

        ## Build the xarray DataArrays for writing. May as well write out the EPA regions, too.
        n_use_o3 = n_obs_o3_use
        n_val_o3 = n_obs_o3_val
        n_use_pm25 = n_obs_pm25_use
        n_val_pm25 = n_obs_pm25_val

        n_obs_use_o3 = n_obs_o3_use
        n_obs_val_o3 = n_obs_o3_val
        n_obs_use_pm25 = n_obs_pm25_use
        n_obs_val_pm25 = n_obs_pm25_val

        n_times_valid = n_valid_times

        cams_ra_raw_prefix = 'Raw CAMS reanalysis'
        cams_fc_raw_prefix = 'Raw CAMS forecast'
        cmaq_fc_raw_prefix = 'Raw CMAQ forecast'
        cams_fc_bm3_prefix = 'BM3 CAMS forecast'
        cmaq_fc_bm3_prefix = 'BM3 CMAQ forecast'
        o3_use_suffix = 'for the use (training) AirNow O3 sites'
        o3_val_suffix = 'for the val (validation) AirNow O3 sites'
        pm25_use_suffix = 'for the use (training) AirNow PM2.5 sites'
        pm25_val_suffix = 'for the val (validation) AirNow PM2.5 sites'

        obs_o3_lat_use = xr.DataArray(obs_o3_lat_use,
                    coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
                    attrs={'description':'Latitudes '+o3_use_suffix})
        obs_o3_lat_val = xr.DataArray(obs_o3_lat_val,
                    coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
                    attrs={'description':'Latitudes '+o3_use_suffix})
        obs_pm25_lat_use = xr.DataArray(obs_pm25_lat_use,
                    coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
                    attrs={'description':'Latitudes '+pm25_use_suffix})
        obs_pm25_lat_val = xr.DataArray(obs_pm25_lat_val,
                    coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
                    attrs={'description':'Latitudes '+pm25_val_suffix})

        obs_o3_lon_use = xr.DataArray(obs_o3_lon_use,
                    coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
                    attrs={'description':'Longitudes '+o3_use_suffix})
        obs_o3_lon_val = xr.DataArray(obs_o3_lon_val,
                    coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
                    attrs={'description':'Longitudes '+o3_use_suffix})
        obs_pm25_lon_use = xr.DataArray(obs_pm25_lon_use,
                    coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
                    attrs={'description':'Longtudes '+pm25_use_suffix})
        obs_pm25_lon_val = xr.DataArray(obs_pm25_lon_val,
                    coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
                    attrs={'description':'Longitudes '+pm25_val_suffix})

        obs_o3_reg_use = xr.DataArray(obs_o3_reg_use,
                    coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
                    attrs={'description':'EPA region (99=outside U.S.) '+o3_use_suffix})
        obs_o3_reg_val = xr.DataArray(obs_o3_reg_val,
                    coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
                    attrs={'description':'EPA region (99=outside U.S.) '+o3_val_suffix})
        obs_pm25_reg_use = xr.DataArray(obs_pm25_reg_use,
                    coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
                    attrs={'description':'EPA region (99=outside U.S.) '+pm25_use_suffix})
        obs_pm25_reg_val = xr.DataArray(obs_pm25_reg_val,
                    coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
                    attrs={'description':'EPA region (99=outside U.S.) '+pm25_val_suffix})

        obs_o3_sid_use = xr.DataArray(obs_o3_sid_use,
                    coords={'n_obs_use_o3':n_use_o3}, dims=['n_use_o3'],
                    attrs={'description':'Station IDs '+o3_use_suffix})
        obs_o3_sid_val = xr.DataArray(obs_o3_sid_val,
                    coords={'n_obs_val_o3':n_val_o3}, dims=['n_val_o3'],
                    attrs={'description':'Station IDs '+o3_use_suffix})
        obs_pm25_sid_use = xr.DataArray(obs_pm25_sid_use,
                    coords={'n_obs_use_pm25':n_use_pm25}, dims=['n_use_pm25'],
                    attrs={'description':'Station IDs '+pm25_use_suffix})
        obs_pm25_sid_val = xr.DataArray(obs_pm25_sid_val,
                    coords={'n_obs_val_pm25':n_val_pm25}, dims=['n_val_pm25'],
                    attrs={'description':'Station IDs '+pm25_val_suffix})

        obs_o3_con_use = xr.DataArray(obs_o3_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_o3':n_use_o3}, dims=['n_times_valid', 'n_use_o3'],
                    attrs={'description':'Observed values of O3 concentration '+o3_use_suffix})
        obs_o3_con_val = xr.DataArray(obs_o3_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_o3':n_val_o3}, dims=['n_times_valid', 'n_val_o3'],
                    attrs={'description':'Observed values of O3 concentration '+o3_use_suffix})
        obs_pm25_con_use = xr.DataArray(obs_pm25_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_pm25':n_use_pm25}, dims=['n_times_valid', 'n_use_pm25'],
                    attrs={'description':'Observed values of O3 concentration '+pm25_use_suffix})
        obs_pm25_con_val = xr.DataArray(obs_pm25_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_pm25':n_val_pm25}, dims=['n_times_valid', 'n_val_pm25'],
                    attrs={'description':'Observed values of O3 concentration '+pm25_val_suffix})

        cams_ra_raw_o3_con_use = xr.DataArray(cams_ra_raw_o3_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_o3':n_use_o3}, dims=['n_times_valid', 'n_use_o3'],
                    attrs={'description':cams_ra_raw_prefix+' values of O3 concentration '+o3_use_suffix})
        cams_ra_raw_o3_con_val = xr.DataArray(cams_ra_raw_o3_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_o3':n_val_o3}, dims=['n_times_valid', 'n_val_o3'],
                    attrs={'description':cams_ra_raw_prefix+' values of O3 concentration '+o3_val_suffix})
        cams_ra_raw_pm25_con_use = xr.DataArray(cams_ra_raw_pm25_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_pm25':n_use_pm25}, dims=['n_times_valid', 'n_use_pm25'],
                    attrs={'description':cams_ra_raw_prefix+' values of PM2.5 concentration '+pm25_use_suffix})
        cams_ra_raw_pm25_con_val = xr.DataArray(cams_ra_raw_pm25_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_pm25':n_val_pm25}, dims=['n_times_valid', 'n_val_pm25'],
                    attrs={'description':cams_ra_raw_prefix+' values of PM2.5 concentration '+pm25_val_suffix})

        cams_fc_raw_o3_con_use = xr.DataArray(cams_fc_raw_o3_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_o3':n_use_o3}, dims=['n_times_valid', 'n_use_o3'],
                    attrs={'description':cams_fc_raw_prefix+' values of O3 concentration '+o3_use_suffix})
        cams_fc_raw_o3_con_val = xr.DataArray(cams_fc_raw_o3_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_o3':n_val_o3}, dims=['n_times_valid', 'n_val_o3'],
                    attrs={'description':cams_fc_raw_prefix+' values of O3 concentration '+o3_val_suffix})
        cams_fc_raw_pm25_con_use = xr.DataArray(cams_fc_raw_pm25_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_pm25':n_use_pm25}, dims=['n_times_valid', 'n_use_pm25'],
                    attrs={'description':cams_fc_raw_prefix+' values of PM2.5 concentration '+pm25_use_suffix})
        cams_fc_raw_pm25_con_val = xr.DataArray(cams_fc_raw_pm25_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_pm25':n_val_pm25}, dims=['n_times_valid', 'n_val_pm25'],
                    attrs={'description':cams_fc_raw_prefix+' values of PM2.5 concentration '+pm25_val_suffix})

        cmaq_fc_raw_o3_con_use = xr.DataArray(cmaq_fc_raw_o3_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_o3':n_use_o3}, dims=['n_times_valid', 'n_use_o3'],
                    attrs={'description':cmaq_fc_raw_prefix+' values of O3 concentration '+o3_use_suffix})
        cmaq_fc_raw_o3_con_val = xr.DataArray(cmaq_fc_raw_o3_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_o3':n_val_o3}, dims=['n_times_valid', 'n_val_o3'],
                    attrs={'description':cmaq_fc_raw_prefix+' values of O3 concentration '+o3_val_suffix})
        cmaq_fc_raw_pm25_con_use = xr.DataArray(cmaq_fc_raw_pm25_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_pm25':n_use_pm25}, dims=['n_times_valid', 'n_use_pm25'],
                    attrs={'description':cmaq_fc_raw_prefix+' values of PM2.5 concentration '+pm25_use_suffix})
        cmaq_fc_raw_pm25_con_val = xr.DataArray(cmaq_fc_raw_pm25_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_pm25':n_val_pm25}, dims=['n_times_valid', 'n_val_pm25'],
                    attrs={'description':cmaq_fc_raw_prefix+' values of PM2.5 concentration '+pm25_val_suffix})

        cams_fc_bm3_o3_con_use = xr.DataArray(cams_fc_bm3_o3_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_o3':n_use_o3}, dims=['n_times_valid', 'n_use_o3'],
                    attrs={'description':cams_fc_bm3_prefix+' values of O3 concentration '+o3_use_suffix})
        cams_fc_bm3_o3_con_val = xr.DataArray(cams_fc_bm3_o3_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_o3':n_val_o3}, dims=['n_times_valid', 'n_val_o3'],
                    attrs={'description':cams_fc_bm3_prefix+' values of O3 concentration '+o3_val_suffix})
        cams_fc_bm3_pm25_con_use = xr.DataArray(cams_fc_bm3_pm25_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_pm25':n_use_pm25}, dims=['n_times_valid', 'n_use_pm25'],
                    attrs={'description':cams_fc_bm3_prefix+' values of PM2.5 concentration '+pm25_use_suffix})
        cams_fc_bm3_pm25_con_val = xr.DataArray(cams_fc_bm3_pm25_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_pm25':n_val_pm25}, dims=['n_times_valid', 'n_val_pm25'],
                    attrs={'description':cams_fc_bm3_prefix+' values of PM2.5 concentration '+pm25_val_suffix})

        cmaq_fc_bm3_o3_con_use = xr.DataArray(cmaq_fc_bm3_o3_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_o3':n_use_o3}, dims=['n_times_valid', 'n_use_o3'],
                    attrs={'description':cmaq_fc_bm3_prefix+' values of O3 concentration '+o3_use_suffix})
        cmaq_fc_bm3_o3_con_val = xr.DataArray(cmaq_fc_bm3_o3_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_o3':n_val_o3}, dims=['n_times_valid', 'n_val_o3'],
                    attrs={'description':cmaq_fc_bm3_prefix+' values of O3 concentration '+o3_val_suffix})
        cmaq_fc_bm3_pm25_con_use = xr.DataArray(cmaq_fc_bm3_pm25_con_use,
                    coords={'n_valid_times':n_times_valid, 'n_obs_use_pm25':n_use_pm25}, dims=['n_times_valid', 'n_use_pm25'],
                    attrs={'description':cmaq_fc_bm3_prefix+' values of PM2.5 concentration '+pm25_use_suffix})
        cmaq_fc_bm3_pm25_con_val = xr.DataArray(cmaq_fc_bm3_pm25_con_val,
                    coords={'n_valid_times':n_times_valid, 'n_obs_val_pm25':n_val_pm25}, dims=['n_times_valid', 'n_val_pm25'],
                    attrs={'description':cmaq_fc_bm3_prefix+' values of PM2.5 concentration '+pm25_val_suffix})

        dt_valid = xr.DataArray(dt_all_valid,
                    coords={'n_valid_times':n_times_valid}, dims=['n_times_valid'],
                    attrs={'description':'Valid date/time stamps [UTC]'})

        ## Create new xarray Datasets
        ds_use = xr.Dataset(
                    data_vars={'dt_valid':dt_valid,
                        'obs_o3_lat_use':obs_o3_lat_use,
                        'obs_o3_lon_use':obs_o3_lon_use,
                        'obs_o3_reg_use':obs_o3_reg_use,
                        'obs_o3_sid_use':obs_o3_sid_use,
                        'obs_o3_con_use':obs_o3_con_use,
                        'cams_ra_raw_o3_con_use':cams_ra_raw_o3_con_use,
                        'cams_fc_raw_o3_con_use':cams_fc_raw_o3_con_use,
                        'cmaq_fc_raw_o3_con_use':cmaq_fc_raw_o3_con_use,
                        'cams_fc_bm3_o3_con_use':cams_fc_bm3_o3_con_use,
                        'cmaq_fc_bm3_o3_con_use':cmaq_fc_bm3_o3_con_use,
                        'obs_pm25_lat_use':obs_pm25_lat_use,
                        'obs_pm25_lon_use':obs_pm25_lon_use,
                        'obs_pm25_reg_use':obs_pm25_reg_use,
                        'obs_pm25_sid_use':obs_pm25_sid_use,
                        'obs_pm25_con_use':obs_pm25_con_use,
                        'cams_ra_raw_pm25_con_use':cams_ra_raw_pm25_con_use,
                        'cams_fc_raw_pm25_con_use':cams_fc_raw_pm25_con_use,
                        'cmaq_fc_raw_pm25_con_use':cmaq_fc_raw_pm25_con_use,
                        'cams_fc_bm3_pm25_con_use':cams_fc_bm3_pm25_con_use,
                        'cmaq_fc_bm3_pm25_con_use':cmaq_fc_bm3_pm25_con_use},
                    coords={'n_obs_use_o3':n_use_o3, 'n_obs_use_pm25':n_use_pm25, 'n_valid_times':n_times_valid},
                    attrs={'description':'Forecasts & observations for AirNow O3 & PM2.5 stations in the use (training) set'}, )

        ds_val = xr.Dataset(
                    data_vars={'dt_valid':dt_valid,
                        'obs_o3_lat_val':obs_o3_lat_val,
                        'obs_o3_lon_val':obs_o3_lon_val,
                        'obs_o3_reg_use':obs_o3_reg_val,
                        'obs_o3_sid_val':obs_o3_sid_val,
                        'obs_o3_con_val':obs_o3_con_val,
                        'cams_ra_raw_o3_con_val':cams_ra_raw_o3_con_val,
                        'cams_fc_raw_o3_con_val':cams_fc_raw_o3_con_val,
                        'cmaq_fc_raw_o3_con_val':cmaq_fc_raw_o3_con_val,
                        'cams_fc_bm3_o3_con_val':cams_fc_bm3_o3_con_val,
                        'cmaq_fc_bm3_o3_con_val':cmaq_fc_bm3_o3_con_val,
                        'obs_pm25_lat_val':obs_pm25_lat_val,
                        'obs_pm25_lon_val':obs_pm25_lon_val,
                        'obs_pm25_reg_val':obs_pm25_reg_val,
                        'obs_pm25_sid_val':obs_pm25_sid_val,
                        'obs_pm25_con_val':obs_pm25_con_val,
                        'cams_ra_raw_pm25_con_val':cams_ra_raw_pm25_con_val,
                        'cams_fc_raw_pm25_con_val':cams_fc_raw_pm25_con_val,
                        'cmaq_fc_raw_pm25_con_val':cmaq_fc_raw_pm25_con_val,
                        'cams_fc_bm3_pm25_con_val':cams_fc_bm3_pm25_con_val,
                        'cmaq_fc_bm3_pm25_con_val':cmaq_fc_bm3_pm25_con_val},
                    coords={'n_obs_val_o3':n_val_o3, 'n_obs_val_pm25':n_val_pm25, 'n_valid_times':n_times_valid},
                    attrs={'description':'Forecasts & observations for AirNow O3 & PM2.5 stations in the use (training) set'}, )

        ## Set the filenames
        fname_use = plot_dir.joinpath('obs_models_o3_pm25_use_'+date_range_file+'.nc')
        fname_val = plot_dir.joinpath('obs_models_o3_pm25_val_'+date_range_file+'.nc')

        ## Write the datasets to NetCDF
        print('Writing '+str(fname_use))
        ds_use.to_netcdf(fname_use)
        fname_use.chmod(0o644)
        print('Writing '+str(fname_val))
        ds_val.to_netcdf(fname_val)
        fname_val.chmod(0o644)

    ## If the pairs files already exist, then that saves a lot of time
    if read_pairs_files:
        fname_use = plot_dir.joinpath('obs_models_o3_pm25_use_'+date_range_file+'.nc')
        fname_val = plot_dir.joinpath('obs_models_o3_pm25_val_'+date_range_file+'.nc')

        print('Reading '+str(fname_use))
        use_ds = xr.open_dataset(fname_use)
        obs_o3_con_use   = use_ds.obs_o3_con_use
        obs_pm25_con_use = use_ds.obs_pm25_con_use
        cams_ra_raw_o3_con_use = use_ds.cams_ra_raw_o3_con_use
        cams_fc_raw_o3_con_use = use_ds.cams_fc_raw_o3_con_use
        cmaq_fc_raw_o3_con_use = use_ds.cmaq_fc_raw_o3_con_use
        cams_fc_bm3_o3_con_use = use_ds.cams_fc_bm3_o3_con_use
        cmaq_fc_bm3_o3_con_use = use_ds.cmaq_fc_bm3_o3_con_use
        cams_ra_raw_pm25_con_use = use_ds.cams_ra_raw_pm25_con_use
        cams_fc_raw_pm25_con_use = use_ds.cams_fc_raw_pm25_con_use
        cmaq_fc_raw_pm25_con_use = use_ds.cmaq_fc_raw_pm25_con_use
        cams_fc_bm3_pm25_con_use = use_ds.cams_fc_bm3_pm25_con_use
        cmaq_fc_bm3_pm25_con_use = use_ds.cmaq_fc_bm3_pm25_con_use

        print('Reading '+str(fname_val))
        val_ds = xr.open_dataset(fname_val)
        obs_o3_con_val   = val_ds.obs_o3_con_val
        obs_pm25_con_val = val_ds.obs_pm25_con_val
        cams_ra_raw_o3_con_val = val_ds.cams_ra_raw_o3_con_val
        cams_fc_raw_o3_con_val = val_ds.cams_fc_raw_o3_con_val
        cmaq_fc_raw_o3_con_val = val_ds.cmaq_fc_raw_o3_con_val
        cams_fc_bm3_o3_con_val = val_ds.cams_fc_bm3_o3_con_val
        cmaq_fc_bm3_o3_con_val = val_ds.cmaq_fc_bm3_o3_con_val
        cams_ra_raw_pm25_con_val = val_ds.cams_ra_raw_pm25_con_val
        cams_fc_raw_pm25_con_val = val_ds.cams_fc_raw_pm25_con_val
        cmaq_fc_raw_pm25_con_val = val_ds.cmaq_fc_raw_pm25_con_val
        cams_fc_bm3_pm25_con_val = val_ds.cams_fc_bm3_pm25_con_val
        cmaq_fc_bm3_pm25_con_val = val_ds.cmaq_fc_bm3_pm25_con_val

    ## Append the arrays
    cams_ra_raw_o3_con_all = np.append(cams_ra_raw_o3_con_use, cams_ra_raw_o3_con_val, axis=1)
    cams_fc_raw_o3_con_all = np.append(cams_fc_raw_o3_con_use, cams_fc_raw_o3_con_val, axis=1)
    cmaq_fc_raw_o3_con_all = np.append(cmaq_fc_raw_o3_con_use, cmaq_fc_raw_o3_con_val, axis=1)
    cams_fc_bm3_o3_con_all = np.append(cams_fc_bm3_o3_con_use, cams_fc_bm3_o3_con_val, axis=1)
    cmaq_fc_bm3_o3_con_all = np.append(cmaq_fc_bm3_o3_con_use, cmaq_fc_bm3_o3_con_val, axis=1)
    cams_ra_raw_pm25_con_all = np.append(cams_ra_raw_pm25_con_use, cams_ra_raw_pm25_con_val, axis=1)
    cams_fc_raw_pm25_con_all = np.append(cams_fc_raw_pm25_con_use, cams_fc_raw_pm25_con_val, axis=1)
    cmaq_fc_raw_pm25_con_all = np.append(cmaq_fc_raw_pm25_con_use, cmaq_fc_raw_pm25_con_val, axis=1)
    cams_fc_bm3_pm25_con_all = np.append(cams_fc_bm3_pm25_con_use, cams_fc_bm3_pm25_con_val, axis=1)
    cmaq_fc_bm3_pm25_con_all = np.append(cmaq_fc_bm3_pm25_con_use, cmaq_fc_bm3_pm25_con_val, axis=1)
    obs_o3_con_all   = np.append(obs_o3_con_use, obs_o3_con_val, axis=1)
    obs_pm25_con_all = np.append(obs_pm25_con_use, obs_pm25_con_val, axis=1)

    ## Create arrays accumulated by diurnal pattern/hour of day
    cams_ra_raw_o3_con_all_diurnal = np.full([24, n_valid_days, n_obs_o3_all], np.nan)
    cams_fc_raw_o3_con_all_diurnal = np.full([24, n_valid_days, n_obs_o3_all], np.nan)
    cmaq_fc_raw_o3_con_all_diurnal = np.full([24, n_valid_days, n_obs_o3_all], np.nan)
    cams_fc_bm3_o3_con_all_diurnal = np.full([24, n_valid_days, n_obs_o3_all], np.nan)
    cmaq_fc_bm3_o3_con_all_diurnal = np.full([24, n_valid_days, n_obs_o3_all], np.nan)
    cams_ra_raw_pm25_con_all_diurnal = np.full([24, n_valid_days, n_obs_pm25_all], np.nan)
    cams_fc_raw_pm25_con_all_diurnal = np.full([24, n_valid_days, n_obs_pm25_all], np.nan)
    cmaq_fc_raw_pm25_con_all_diurnal = np.full([24, n_valid_days, n_obs_pm25_all], np.nan)
    cams_fc_bm3_pm25_con_all_diurnal = np.full([24, n_valid_days, n_obs_pm25_all], np.nan)
    cmaq_fc_bm3_pm25_con_all_diurnal = np.full([24, n_valid_days, n_obs_pm25_all], np.nan)
    obs_o3_con_all_diurnal   = np.full([24, n_valid_days, n_obs_o3_all], np.nan)
    obs_pm25_con_all_diurnal = np.full([24, n_valid_days, n_obs_pm25_all], np.nan)

    cams_ra_raw_o3_con_val_diurnal = np.full([24, n_valid_days, n_obs_o3_val], np.nan)
    cams_fc_raw_o3_con_val_diurnal = np.full([24, n_valid_days, n_obs_o3_val], np.nan)
    cmaq_fc_raw_o3_con_val_diurnal = np.full([24, n_valid_days, n_obs_o3_val], np.nan)
    cams_fc_bm3_o3_con_val_diurnal = np.full([24, n_valid_days, n_obs_o3_val], np.nan)
    cmaq_fc_bm3_o3_con_val_diurnal = np.full([24, n_valid_days, n_obs_o3_val], np.nan)
    cams_ra_raw_pm25_con_val_diurnal = np.full([24, n_valid_days, n_obs_pm25_val], np.nan)
    cams_fc_raw_pm25_con_val_diurnal = np.full([24, n_valid_days, n_obs_pm25_val], np.nan)
    cmaq_fc_raw_pm25_con_val_diurnal = np.full([24, n_valid_days, n_obs_pm25_val], np.nan)
    cams_fc_bm3_pm25_con_val_diurnal = np.full([24, n_valid_days, n_obs_pm25_val], np.nan)
    cmaq_fc_bm3_pm25_con_val_diurnal = np.full([24, n_valid_days, n_obs_pm25_val], np.nan)
    obs_o3_con_val_diurnal   = np.full([24, n_valid_days, n_obs_o3_val], np.nan)
    obs_pm25_con_val_diurnal = np.full([24, n_valid_days, n_obs_pm25_val], np.nan)

    ## Populate the diurnal arrays
    print('Populating the diurnal arrays')
    for hh in range(24):
        cams_ra_raw_o3_con_all_diurnal[hh,:,:] = cams_ra_raw_o3_con_all[hh::24, :]
        cams_fc_raw_o3_con_all_diurnal[hh,:,:] = cams_fc_raw_o3_con_all[hh::24, :]
        cmaq_fc_raw_o3_con_all_diurnal[hh,:,:] = cmaq_fc_raw_o3_con_all[hh::24, :]
        cams_fc_bm3_o3_con_all_diurnal[hh,:,:] = cams_fc_bm3_o3_con_all[hh::24, :]
        cmaq_fc_bm3_o3_con_all_diurnal[hh,:,:] = cmaq_fc_bm3_o3_con_all[hh::24, :]
        cams_ra_raw_pm25_con_all_diurnal[hh,:,:] = cams_ra_raw_pm25_con_all[hh::24, :]
        cams_fc_raw_pm25_con_all_diurnal[hh,:,:] = cams_fc_raw_pm25_con_all[hh::24, :]
        cmaq_fc_raw_pm25_con_all_diurnal[hh,:,:] = cmaq_fc_raw_pm25_con_all[hh::24, :]
        cams_fc_bm3_pm25_con_all_diurnal[hh,:,:] = cams_fc_bm3_pm25_con_all[hh::24, :]
        cmaq_fc_bm3_pm25_con_all_diurnal[hh,:,:] = cmaq_fc_bm3_pm25_con_all[hh::24, :]
        obs_o3_con_all_diurnal[hh,:,:]   = obs_o3_con_all[hh::24, :]
        obs_pm25_con_all_diurnal[hh,:,:] = obs_pm25_con_all[hh::24, :]

        cams_ra_raw_o3_con_val_diurnal[hh,:,:] = cams_ra_raw_o3_con_val[hh::24, :]
        cams_fc_raw_o3_con_val_diurnal[hh,:,:] = cams_fc_raw_o3_con_val[hh::24, :]
        cmaq_fc_raw_o3_con_val_diurnal[hh,:,:] = cmaq_fc_raw_o3_con_val[hh::24, :]
        cams_fc_bm3_o3_con_val_diurnal[hh,:,:] = cams_fc_bm3_o3_con_val[hh::24, :]
        cmaq_fc_bm3_o3_con_val_diurnal[hh,:,:] = cmaq_fc_bm3_o3_con_val[hh::24, :]
        cams_ra_raw_pm25_con_val_diurnal[hh,:,:] = cams_ra_raw_pm25_con_val[hh::24, :]
        cams_fc_raw_pm25_con_val_diurnal[hh,:,:] = cams_fc_raw_pm25_con_val[hh::24, :]
        cmaq_fc_raw_pm25_con_val_diurnal[hh,:,:] = cmaq_fc_raw_pm25_con_val[hh::24, :]
        cams_fc_bm3_pm25_con_val_diurnal[hh,:,:] = cams_fc_bm3_pm25_con_val[hh::24, :]
        cmaq_fc_bm3_pm25_con_val_diurnal[hh,:,:] = cmaq_fc_bm3_pm25_con_val[hh::24, :]
        obs_o3_con_val_diurnal[hh,:,:]   = obs_o3_con_val[hh::24, :]
        obs_pm25_con_val_diurnal[hh,:,:] = obs_pm25_con_val[hh::24, :]

    ## Create arrays accumulated as daily-average values
    print('Calculating daily averages of model & obs')
    cams_ra_raw_o3_con_all_daily   = np.full([n_valid_days, n_obs_o3_all], np.nan)
    cams_fc_raw_o3_con_all_daily   = np.full([n_valid_days, n_obs_o3_all], np.nan)
    cmaq_fc_raw_o3_con_all_daily   = np.full([n_valid_days, n_obs_o3_all], np.nan)
    cams_fc_bm3_o3_con_all_daily   = np.full([n_valid_days, n_obs_o3_all], np.nan)
    cmaq_fc_bm3_o3_con_all_daily   = np.full([n_valid_days, n_obs_o3_all], np.nan)
    cams_ra_raw_pm25_con_all_daily = np.full([n_valid_days, n_obs_pm25_all], np.nan)
    cams_fc_raw_pm25_con_all_daily = np.full([n_valid_days, n_obs_pm25_all], np.nan)
    cmaq_fc_raw_pm25_con_all_daily = np.full([n_valid_days, n_obs_pm25_all], np.nan)
    cams_fc_bm3_pm25_con_all_daily = np.full([n_valid_days, n_obs_pm25_all], np.nan)
    cmaq_fc_bm3_pm25_con_all_daily = np.full([n_valid_days, n_obs_pm25_all], np.nan)
    obs_o3_con_all_daily   = np.full([n_valid_days, n_obs_o3_all], np.nan)
    obs_pm25_con_all_daily = np.full([n_valid_days, n_obs_pm25_all], np.nan)

    cams_ra_raw_o3_con_val_daily   = np.full([n_valid_days, n_obs_o3_val], np.nan)
    cams_fc_raw_o3_con_val_daily   = np.full([n_valid_days, n_obs_o3_val], np.nan)
    cmaq_fc_raw_o3_con_val_daily   = np.full([n_valid_days, n_obs_o3_val], np.nan)
    cams_fc_bm3_o3_con_val_daily   = np.full([n_valid_days, n_obs_o3_val], np.nan)
    cmaq_fc_bm3_o3_con_val_daily   = np.full([n_valid_days, n_obs_o3_val], np.nan)
    cams_ra_raw_pm25_con_val_daily = np.full([n_valid_days, n_obs_pm25_val], np.nan)
    cams_fc_raw_pm25_con_val_daily = np.full([n_valid_days, n_obs_pm25_val], np.nan)
    cmaq_fc_raw_pm25_con_val_daily = np.full([n_valid_days, n_obs_pm25_val], np.nan)
    cams_fc_bm3_pm25_con_val_daily = np.full([n_valid_days, n_obs_pm25_val], np.nan)
    cmaq_fc_bm3_pm25_con_val_daily = np.full([n_valid_days, n_obs_pm25_val], np.nan)
    obs_o3_con_val_daily   = np.full([n_valid_days, n_obs_o3_val], np.nan)
    obs_pm25_con_val_daily = np.full([n_valid_days, n_obs_pm25_val], np.nan)

    ## Populate the daily-average arrays
    for dd in range(n_valid_days):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            cams_ra_raw_o3_con_all_daily[dd,:] = np.nanmean(cams_ra_raw_o3_con_all[dd*24:dd*24+24, :], axis=0)
            cams_fc_raw_o3_con_all_daily[dd,:] = np.nanmean(cams_fc_raw_o3_con_all[dd*24:dd*24+24, :], axis=0)
            cmaq_fc_raw_o3_con_all_daily[dd,:] = np.nanmean(cmaq_fc_raw_o3_con_all[dd*24:dd*24+24, :], axis=0)
            cams_fc_bm3_o3_con_all_daily[dd,:] = np.nanmean(cams_fc_bm3_o3_con_all[dd*24:dd*24+24, :], axis=0)
            cmaq_fc_bm3_o3_con_all_daily[dd,:] = np.nanmean(cmaq_fc_bm3_o3_con_all[dd*24:dd*24+24, :], axis=0)
            cams_ra_raw_pm25_con_all_daily[dd,:] = np.nanmean(cams_ra_raw_pm25_con_all[dd*24:dd*24+24, :], axis=0)
            cams_fc_raw_pm25_con_all_daily[dd,:] = np.nanmean(cams_fc_raw_pm25_con_all[dd*24:dd*24+24, :], axis=0)
            cmaq_fc_raw_pm25_con_all_daily[dd,:] = np.nanmean(cmaq_fc_raw_pm25_con_all[dd*24:dd*24+24, :], axis=0)
            cams_fc_bm3_pm25_con_all_daily[dd,:] = np.nanmean(cams_fc_bm3_pm25_con_all[dd*24:dd*24+24, :], axis=0)
            cmaq_fc_bm3_pm25_con_all_daily[dd,:] = np.nanmean(cmaq_fc_bm3_pm25_con_all[dd*24:dd*24+24, :], axis=0)
            obs_o3_con_all_daily[dd,:]   = np.nanmean(obs_o3_con_all[dd*24:dd*24+24, :], axis=0)
            obs_pm25_con_all_daily[dd,:] = np.nanmean(obs_pm25_con_all[dd*24:dd*24+24, :], axis=0)

            cams_ra_raw_o3_con_val_daily[dd,:] = np.nanmean(cams_ra_raw_o3_con_val[dd*24:dd*24+24, :], axis=0)
            cams_fc_raw_o3_con_val_daily[dd,:] = np.nanmean(cams_fc_raw_o3_con_val[dd*24:dd*24+24, :], axis=0)
            cmaq_fc_raw_o3_con_val_daily[dd,:] = np.nanmean(cmaq_fc_raw_o3_con_val[dd*24:dd*24+24, :], axis=0)
            cams_fc_bm3_o3_con_val_daily[dd,:] = np.nanmean(cams_fc_bm3_o3_con_val[dd*24:dd*24+24, :], axis=0)
            cmaq_fc_bm3_o3_con_val_daily[dd,:] = np.nanmean(cmaq_fc_bm3_o3_con_val[dd*24:dd*24+24, :], axis=0)
            cams_ra_raw_pm25_con_val_daily[dd,:] = np.nanmean(cams_ra_raw_pm25_con_val[dd*24:dd*24+24, :], axis=0)
            cams_fc_raw_pm25_con_val_daily[dd,:] = np.nanmean(cams_fc_raw_pm25_con_val[dd*24:dd*24+24, :], axis=0)
            cmaq_fc_raw_pm25_con_val_daily[dd,:] = np.nanmean(cmaq_fc_raw_pm25_con_val[dd*24:dd*24+24, :], axis=0)
            cams_fc_bm3_pm25_con_val_daily[dd,:] = np.nanmean(cams_fc_bm3_pm25_con_val[dd*24:dd*24+24, :], axis=0)
            cmaq_fc_bm3_pm25_con_val_daily[dd,:] = np.nanmean(cmaq_fc_bm3_pm25_con_val[dd*24:dd*24+24, :], axis=0)
            obs_o3_con_val_daily[dd,:]   = np.nanmean(obs_o3_con_val[dd*24:dd*24+24, :], axis=0)
            obs_pm25_con_val_daily[dd,:] = np.nanmean(obs_pm25_con_val[dd*24:dd*24+24, :], axis=0)

    ## Calculate statistics for each region (overall=0, 1-10, outside U.S.=11)
    ## But first, define the arrays
    print('Calculating diurnal statistics')
    cams_ra_raw_o3_con_all_diurnal_mean = np.full([24, 12], np.nan)
    cams_fc_raw_o3_con_all_diurnal_mean = np.full([24, 12], np.nan)
    cmaq_fc_raw_o3_con_all_diurnal_mean = np.full([24, 12], np.nan)
    cams_fc_bm3_o3_con_all_diurnal_mean = np.full([24, 12], np.nan)
    cmaq_fc_bm3_o3_con_all_diurnal_mean = np.full([24, 12], np.nan)
    cams_ra_raw_pm25_con_all_diurnal_mean = np.full([24, 12], np.nan)
    cams_fc_raw_pm25_con_all_diurnal_mean = np.full([24, 12], np.nan)
    cmaq_fc_raw_pm25_con_all_diurnal_mean = np.full([24, 12], np.nan)
    cams_fc_bm3_pm25_con_all_diurnal_mean = np.full([24, 12], np.nan)
    cmaq_fc_bm3_pm25_con_all_diurnal_mean = np.full([24, 12], np.nan)
    obs_o3_con_all_diurnal_mean   = np.full([24, 12], np.nan)
    obs_pm25_con_all_diurnal_mean = np.full([24, 12], np.nan)

    cams_ra_raw_o3_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    cams_fc_raw_o3_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    cmaq_fc_raw_o3_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    cams_fc_bm3_o3_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    cmaq_fc_bm3_o3_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    cams_ra_raw_pm25_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    cams_fc_raw_pm25_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    cmaq_fc_raw_pm25_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    cams_fc_bm3_pm25_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    cmaq_fc_bm3_pm25_con_all_diurnal_sdev = np.full([24, 12], np.nan)
    obs_o3_con_all_diurnal_sdev   = np.full([24, 12], np.nan)
    obs_pm25_con_all_diurnal_sdev = np.full([24, 12], np.nan)

    cams_ra_raw_o3_con_all_diurnal_rmse = np.full([24, 12], np.nan)
    cams_fc_raw_o3_con_all_diurnal_rmse = np.full([24, 12], np.nan)
    cmaq_fc_raw_o3_con_all_diurnal_rmse = np.full([24, 12], np.nan)
    cams_fc_bm3_o3_con_all_diurnal_rmse = np.full([24, 12], np.nan)
    cmaq_fc_bm3_o3_con_all_diurnal_rmse = np.full([24, 12], np.nan)
    cams_ra_raw_pm25_con_all_diurnal_rmse = np.full([24, 12], np.nan)
    cams_fc_raw_pm25_con_all_diurnal_rmse = np.full([24, 12], np.nan)
    cmaq_fc_raw_pm25_con_all_diurnal_rmse = np.full([24, 12], np.nan)
    cams_fc_bm3_pm25_con_all_diurnal_rmse = np.full([24, 12], np.nan)
    cmaq_fc_bm3_pm25_con_all_diurnal_rmse = np.full([24, 12], np.nan)

    cams_ra_raw_o3_con_all_diurnal_mbe = np.full([24, 12], np.nan)
    cams_fc_raw_o3_con_all_diurnal_mbe = np.full([24, 12], np.nan)
    cmaq_fc_raw_o3_con_all_diurnal_mbe = np.full([24, 12], np.nan)
    cams_fc_bm3_o3_con_all_diurnal_mbe = np.full([24, 12], np.nan)
    cmaq_fc_bm3_o3_con_all_diurnal_mbe = np.full([24, 12], np.nan)
    cams_ra_raw_pm25_con_all_diurnal_mbe = np.full([24, 12], np.nan)
    cams_fc_raw_pm25_con_all_diurnal_mbe = np.full([24, 12], np.nan)
    cmaq_fc_raw_pm25_con_all_diurnal_mbe = np.full([24, 12], np.nan)
    cams_fc_bm3_pm25_con_all_diurnal_mbe = np.full([24, 12], np.nan)
    cmaq_fc_bm3_pm25_con_all_diurnal_mbe = np.full([24, 12], np.nan)

    ## Overall
    axis = (1,2)
    obs = obs_o3_con_all_diurnal[:,:,:]
    obs_o3_con_all_diurnal_mean[:, 0] = np.nanmean(obs, axis=axis)
    obs_o3_con_all_diurnal_sdev[:, 0] = np.nanstd( obs, axis=axis)
    fcst = cams_ra_raw_o3_con_all_diurnal[:,:,:]
    cams_ra_raw_o3_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_ra_raw_o3_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_ra_raw_o3_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_ra_raw_o3_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_raw_o3_con_all_diurnal[:,:,:]
    cams_fc_raw_o3_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_raw_o3_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_raw_o3_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_raw_o3_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_raw_o3_con_all_diurnal[:,:,:]
    cmaq_fc_raw_o3_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_raw_o3_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cmaq_fc_raw_o3_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cmaq_fc_raw_o3_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_bm3_o3_con_all_diurnal[:,:,:]
    cams_fc_bm3_o3_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_o3_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_bm3_o3_con_all_diurnal[:,:,:]
    cmaq_fc_bm3_o3_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_bm3_o3_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_o3_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)

    obs = obs_pm25_con_all_diurnal[:,:,:]
    obs_pm25_con_all_diurnal_mean[:, 0] = np.nanmean(obs, axis=axis)
    obs_pm25_con_all_diurnal_sdev[:, 0] = np.nanstd( obs, axis=axis)
    fcst = cams_ra_raw_pm25_con_all_diurnal[:,:,:]
    cams_ra_raw_pm25_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_ra_raw_pm25_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_ra_raw_pm25_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_ra_raw_pm25_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_raw_pm25_con_all_diurnal[:,:,:]
    cams_fc_raw_pm25_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_raw_pm25_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_raw_pm25_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_raw_pm25_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_raw_pm25_con_all_diurnal[:,:,:]
    cmaq_fc_raw_pm25_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_raw_pm25_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cmaq_fc_raw_pm25_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cmaq_fc_raw_pm25_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_bm3_pm25_con_all_diurnal[:,:,:]
    cams_fc_bm3_pm25_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_pm25_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_bm3_pm25_con_all_diurnal[:,:,:]
    cmaq_fc_bm3_pm25_con_all_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_bm3_pm25_con_all_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_all_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_pm25_con_all_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)

    ## Regions 1-10, 99
    for rr in range(1,12):
        if rr == 1:
            inds_o3_reg = inds_o3_reg01_all
            inds_pm25_reg = inds_pm25_reg01_all
        elif rr == 2:
            inds_o3_reg = inds_o3_reg02_all
            inds_pm25_reg = inds_pm25_reg02_all
        elif rr == 3:
            inds_o3_reg = inds_o3_reg03_all
            inds_pm25_reg = inds_pm25_reg03_all
        elif rr == 4:
            inds_o3_reg = inds_o3_reg04_all
            inds_pm25_reg = inds_pm25_reg04_all
        elif rr == 5:
            inds_o3_reg = inds_o3_reg05_all
            inds_pm25_reg = inds_pm25_reg05_all
        elif rr == 6:
            inds_o3_reg = inds_o3_reg06_all
            inds_pm25_reg = inds_pm25_reg06_all
        elif rr == 7:
            inds_o3_reg = inds_o3_reg07_all
            inds_pm25_reg = inds_pm25_reg07_all
        elif rr == 8:
            inds_o3_reg = inds_o3_reg08_all
            inds_pm25_reg = inds_pm25_reg08_all
        elif rr == 9:
            inds_o3_reg = inds_o3_reg09_all
            inds_pm25_reg = inds_pm25_reg09_all
        elif rr == 10:
            inds_o3_reg = inds_o3_reg10_all
            inds_pm25_reg = inds_pm25_reg10_all
        elif rr == 11:
            inds_o3_reg = inds_o3_reg99_all
            inds_pm25_reg = inds_pm25_reg99_all

        obs = obs_o3_con_all_diurnal[:,:,inds_o3_reg]
        obs_o3_con_all_diurnal_mean[:, rr] = np.nanmean(obs, axis=axis)
        obs_o3_con_all_diurnal_sdev[:, rr] = np.nanstd( obs, axis=axis)
        fcst = cams_ra_raw_o3_con_all_diurnal[:,:,inds_o3_reg]
        cams_ra_raw_o3_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_ra_raw_o3_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_ra_raw_o3_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_ra_raw_o3_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_raw_o3_con_all_diurnal[:,:,inds_o3_reg]
        cams_fc_raw_o3_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_raw_o3_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_raw_o3_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_raw_o3_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_raw_o3_con_all_diurnal[:,:,inds_o3_reg]
        cmaq_fc_raw_o3_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_raw_o3_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_raw_o3_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_raw_o3_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_bm3_o3_con_all_diurnal[:,:,inds_o3_reg]
        cams_fc_bm3_o3_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_bm3_o3_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_bm3_o3_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_bm3_o3_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_bm3_o3_con_all_diurnal[:,:,inds_o3_reg]
        cmaq_fc_bm3_o3_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_bm3_o3_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_bm3_o3_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_bm3_o3_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)

        obs = obs_pm25_con_all_diurnal[:,:,inds_pm25_reg]
        obs_pm25_con_all_diurnal_mean[:, rr] = np.nanmean(obs, axis=axis)
        obs_pm25_con_all_diurnal_sdev[:, rr] = np.nanstd( obs, axis=axis)
        fcst = cams_ra_raw_pm25_con_all_diurnal[:,:,inds_pm25_reg]
        cams_ra_raw_pm25_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_ra_raw_pm25_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_ra_raw_pm25_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_ra_raw_pm25_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_raw_pm25_con_all_diurnal[:,:,inds_pm25_reg]
        cams_fc_raw_pm25_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_raw_pm25_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_raw_pm25_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_raw_pm25_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_raw_pm25_con_all_diurnal[:,:,inds_pm25_reg]
        cmaq_fc_raw_pm25_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_raw_pm25_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_raw_pm25_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_raw_pm25_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_bm3_pm25_con_all_diurnal[:,:,inds_pm25_reg]
        cams_fc_bm3_pm25_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_bm3_pm25_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_bm3_pm25_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_bm3_pm25_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_bm3_pm25_con_all_diurnal[:,:,inds_pm25_reg]
        cmaq_fc_bm3_pm25_con_all_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_bm3_pm25_con_all_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_bm3_pm25_con_all_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_bm3_pm25_con_all_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)

    ## Now do the same but only over the validation stations
    cams_ra_raw_o3_con_val_diurnal_mean = np.full([24, 12], np.nan)
    cams_fc_raw_o3_con_val_diurnal_mean = np.full([24, 12], np.nan)
    cmaq_fc_raw_o3_con_val_diurnal_mean = np.full([24, 12], np.nan)
    cams_fc_bm3_o3_con_val_diurnal_mean = np.full([24, 12], np.nan)
    cmaq_fc_bm3_o3_con_val_diurnal_mean = np.full([24, 12], np.nan)
    cams_ra_raw_pm25_con_val_diurnal_mean = np.full([24, 12], np.nan)
    cams_fc_raw_pm25_con_val_diurnal_mean = np.full([24, 12], np.nan)
    cmaq_fc_raw_pm25_con_val_diurnal_mean = np.full([24, 12], np.nan)
    cams_fc_bm3_pm25_con_val_diurnal_mean = np.full([24, 12], np.nan)
    cmaq_fc_bm3_pm25_con_val_diurnal_mean = np.full([24, 12], np.nan)
    obs_o3_con_val_diurnal_mean   = np.full([24, 12], np.nan)
    obs_pm25_con_val_diurnal_mean = np.full([24, 12], np.nan)

    cams_ra_raw_o3_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    cams_fc_raw_o3_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    cmaq_fc_raw_o3_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    cams_fc_bm3_o3_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    cmaq_fc_bm3_o3_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    cams_ra_raw_pm25_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    cams_fc_raw_pm25_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    cmaq_fc_raw_pm25_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    cams_fc_bm3_pm25_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    cmaq_fc_bm3_pm25_con_val_diurnal_sdev = np.full([24, 12], np.nan)
    obs_o3_con_val_diurnal_sdev   = np.full([24, 12], np.nan)
    obs_pm25_con_val_diurnal_sdev = np.full([24, 12], np.nan)

    cams_ra_raw_o3_con_val_diurnal_rmse = np.full([24, 12], np.nan)
    cams_fc_raw_o3_con_val_diurnal_rmse = np.full([24, 12], np.nan)
    cmaq_fc_raw_o3_con_val_diurnal_rmse = np.full([24, 12], np.nan)
    cams_fc_bm3_o3_con_val_diurnal_rmse = np.full([24, 12], np.nan)
    cmaq_fc_bm3_o3_con_val_diurnal_rmse = np.full([24, 12], np.nan)
    cams_ra_raw_pm25_con_val_diurnal_rmse = np.full([24, 12], np.nan)
    cams_fc_raw_pm25_con_val_diurnal_rmse = np.full([24, 12], np.nan)
    cmaq_fc_raw_pm25_con_val_diurnal_rmse = np.full([24, 12], np.nan)
    cams_fc_bm3_pm25_con_val_diurnal_rmse = np.full([24, 12], np.nan)
    cmaq_fc_bm3_pm25_con_val_diurnal_rmse = np.full([24, 12], np.nan)

    cams_ra_raw_o3_con_val_diurnal_mbe = np.full([24, 12], np.nan)
    cams_fc_raw_o3_con_val_diurnal_mbe = np.full([24, 12], np.nan)
    cmaq_fc_raw_o3_con_val_diurnal_mbe = np.full([24, 12], np.nan)
    cams_fc_bm3_o3_con_val_diurnal_mbe = np.full([24, 12], np.nan)
    cmaq_fc_bm3_o3_con_val_diurnal_mbe = np.full([24, 12], np.nan)
    cams_ra_raw_pm25_con_val_diurnal_mbe = np.full([24, 12], np.nan)
    cams_fc_raw_pm25_con_val_diurnal_mbe = np.full([24, 12], np.nan)
    cmaq_fc_raw_pm25_con_val_diurnal_mbe = np.full([24, 12], np.nan)
    cams_fc_bm3_pm25_con_val_diurnal_mbe = np.full([24, 12], np.nan)
    cmaq_fc_bm3_pm25_con_val_diurnal_mbe = np.full([24, 12], np.nan)

    ## Overall
    axis = (1,2)
    obs = obs_o3_con_val_diurnal[:,:,:]
    obs_o3_con_val_diurnal_mean[:, 0] = np.nanmean(obs, axis=axis)
    obs_o3_con_val_diurnal_sdev[:, 0] = np.nanstd( obs, axis=axis)
    fcst = cams_ra_raw_o3_con_val_diurnal[:,:,:]
    cams_ra_raw_o3_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_ra_raw_o3_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_ra_raw_o3_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_ra_raw_o3_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_raw_o3_con_val_diurnal[:,:,:]
    cams_fc_raw_o3_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_raw_o3_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_raw_o3_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_raw_o3_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_raw_o3_con_val_diurnal[:,:,:]
    cmaq_fc_raw_o3_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_raw_o3_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cmaq_fc_raw_o3_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cmaq_fc_raw_o3_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_bm3_o3_con_val_diurnal[:,:,:]
    cams_fc_bm3_o3_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_o3_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_bm3_o3_con_val_diurnal[:,:,:]
    cmaq_fc_bm3_o3_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_bm3_o3_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_o3_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)

    obs = obs_pm25_con_val_diurnal[:,:,:]
    obs_pm25_con_val_diurnal_mean[:, 0] = np.nanmean(obs, axis=axis)
    obs_pm25_con_val_diurnal_sdev[:, 0] = np.nanstd( obs, axis=axis)
    fcst = cams_ra_raw_pm25_con_val_diurnal[:,:,:]
    cams_ra_raw_pm25_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_ra_raw_pm25_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_ra_raw_pm25_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_ra_raw_pm25_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_raw_pm25_con_val_diurnal[:,:,:]
    cams_fc_raw_pm25_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_raw_pm25_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_raw_pm25_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_raw_pm25_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_raw_pm25_con_val_diurnal[:,:,:]
    cmaq_fc_raw_pm25_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_raw_pm25_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cmaq_fc_raw_pm25_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cmaq_fc_raw_pm25_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_bm3_pm25_con_val_diurnal[:,:,:]
    cams_fc_bm3_pm25_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_pm25_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_bm3_pm25_con_val_diurnal[:,:,:]
    cmaq_fc_bm3_pm25_con_val_diurnal_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_bm3_pm25_con_val_diurnal_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_val_diurnal_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_pm25_con_val_diurnal_sdev[:, 0] = np.nanstd( fcst, axis=axis)

    ## Regions 1-10, 99
    for rr in range(1,12):
        if rr == 1:
            inds_o3_reg = inds_o3_reg01_val
            inds_pm25_reg = inds_pm25_reg01_val
        elif rr == 2:
            inds_o3_reg = inds_o3_reg02_val
            inds_pm25_reg = inds_pm25_reg02_val
        elif rr == 3:
            inds_o3_reg = inds_o3_reg03_val
            inds_pm25_reg = inds_pm25_reg03_val
        elif rr == 4:
            inds_o3_reg = inds_o3_reg04_val
            inds_pm25_reg = inds_pm25_reg04_val
        elif rr == 5:
            inds_o3_reg = inds_o3_reg05_val
            inds_pm25_reg = inds_pm25_reg05_val
        elif rr == 6:
            inds_o3_reg = inds_o3_reg06_val
            inds_pm25_reg = inds_pm25_reg06_val
        elif rr == 7:
            inds_o3_reg = inds_o3_reg07_val
            inds_pm25_reg = inds_pm25_reg07_val
        elif rr == 8:
            inds_o3_reg = inds_o3_reg08_val
            inds_pm25_reg = inds_pm25_reg08_val
        elif rr == 9:
            inds_o3_reg = inds_o3_reg09_val
            inds_pm25_reg = inds_pm25_reg09_val
        elif rr == 10:
            inds_o3_reg = inds_o3_reg10_val
            inds_pm25_reg = inds_pm25_reg10_val
        elif rr == 11:
            inds_o3_reg = inds_o3_reg99_val
            inds_pm25_reg = inds_pm25_reg99_val

        obs = obs_o3_con_val_diurnal[:,:,inds_o3_reg]
        obs_o3_con_val_diurnal_mean[:, rr] = np.nanmean(obs, axis=axis)
        obs_o3_con_val_diurnal_sdev[:, rr] = np.nanstd( obs, axis=axis)
        fcst = cams_ra_raw_o3_con_val_diurnal[:,:,inds_o3_reg]
        cams_ra_raw_o3_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_ra_raw_o3_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_ra_raw_o3_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_ra_raw_o3_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_raw_o3_con_val_diurnal[:,:,inds_o3_reg]
        cams_fc_raw_o3_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_raw_o3_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_raw_o3_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_raw_o3_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_raw_o3_con_val_diurnal[:,:,inds_o3_reg]
        cmaq_fc_raw_o3_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_raw_o3_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_raw_o3_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_raw_o3_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_bm3_o3_con_val_diurnal[:,:,inds_o3_reg]
        cams_fc_bm3_o3_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_bm3_o3_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_bm3_o3_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_bm3_o3_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_bm3_o3_con_val_diurnal[:,:,inds_o3_reg]
        cmaq_fc_bm3_o3_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_bm3_o3_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_bm3_o3_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_bm3_o3_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)

        obs = obs_pm25_con_val_diurnal[:,:,inds_pm25_reg]
        obs_pm25_con_val_diurnal_mean[:, rr] = np.nanmean(obs, axis=axis)
        obs_pm25_con_val_diurnal_sdev[:, rr] = np.nanstd( obs, axis=axis)
        fcst = cams_ra_raw_pm25_con_val_diurnal[:,:,inds_pm25_reg]
        cams_ra_raw_pm25_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_ra_raw_pm25_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_ra_raw_pm25_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_ra_raw_pm25_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_raw_pm25_con_val_diurnal[:,:,inds_pm25_reg]
        cams_fc_raw_pm25_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_raw_pm25_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_raw_pm25_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_raw_pm25_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_raw_pm25_con_val_diurnal[:,:,inds_pm25_reg]
        cmaq_fc_raw_pm25_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_raw_pm25_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_raw_pm25_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_raw_pm25_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_bm3_pm25_con_val_diurnal[:,:,inds_pm25_reg]
        cams_fc_bm3_pm25_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_bm3_pm25_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_bm3_pm25_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_bm3_pm25_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_bm3_pm25_con_val_diurnal[:,:,inds_pm25_reg]
        cmaq_fc_bm3_pm25_con_val_diurnal_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_bm3_pm25_con_val_diurnal_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_bm3_pm25_con_val_diurnal_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_bm3_pm25_con_val_diurnal_sdev[:, rr] = np.nanstd( fcst, axis=axis)

    ## Daily-mean stats
    ## Define arrays first
    print('Calculating daily-mean statistics')
    cams_ra_raw_o3_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_o3_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_o3_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_o3_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_o3_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_ra_raw_pm25_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_pm25_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_pm25_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_pm25_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_pm25_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)
    obs_o3_con_all_daily_mean   = np.full([n_valid_days, 12], np.nan)
    obs_pm25_con_all_daily_mean = np.full([n_valid_days, 12], np.nan)

    cams_ra_raw_o3_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_o3_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_o3_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_o3_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_o3_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_ra_raw_pm25_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_pm25_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_pm25_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_pm25_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_pm25_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)
    obs_o3_con_all_daily_sdev   = np.full([n_valid_days, 12], np.nan)
    obs_pm25_con_all_daily_sdev = np.full([n_valid_days, 12], np.nan)

    cams_ra_raw_o3_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_o3_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_o3_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_o3_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_o3_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_ra_raw_pm25_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_pm25_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_pm25_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_pm25_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_pm25_con_all_daily_rmse = np.full([n_valid_days, 12], np.nan)

    cams_ra_raw_o3_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_o3_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_o3_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_o3_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_o3_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_ra_raw_pm25_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_pm25_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_pm25_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_pm25_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_pm25_con_all_daily_mbe = np.full([n_valid_days, 12], np.nan)

    ## Overall
    axis = 1
    obs = obs_o3_con_all_daily[:,:]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        obs_o3_con_all_daily_mean[:, 0] = np.nanmean(obs, axis=axis)
        obs_o3_con_all_daily_sdev[:, 0] = np.nanstd( obs, axis=axis)
    fcst = cams_ra_raw_o3_con_all_daily[:,:]
    cams_ra_raw_o3_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_ra_raw_o3_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_ra_raw_o3_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_ra_raw_o3_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_raw_o3_con_all_daily[:,:]
    cams_fc_raw_o3_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_raw_o3_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_raw_o3_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_raw_o3_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_raw_o3_con_all_daily[:,:]
    cmaq_fc_raw_o3_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_raw_o3_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cmaq_fc_raw_o3_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cmaq_fc_raw_o3_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_bm3_o3_con_all_daily[:,:]
    cams_fc_bm3_o3_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_o3_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_bm3_o3_con_all_daily[:,:]
    cmaq_fc_bm3_o3_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_bm3_o3_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_o3_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)

    obs = obs_pm25_con_all_daily[:,:]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        obs_pm25_con_all_daily_mean[:, 0] = np.nanmean(obs, axis=axis)
        obs_pm25_con_all_daily_sdev[:, 0] = np.nanstd( obs, axis=axis)
    fcst = cams_ra_raw_pm25_con_all_daily[:,:]
    cams_ra_raw_pm25_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_ra_raw_pm25_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_ra_raw_pm25_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_ra_raw_pm25_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_raw_pm25_con_all_daily[:,:]
    cams_fc_raw_pm25_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_raw_pm25_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_raw_pm25_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_raw_pm25_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_raw_pm25_con_all_daily[:,:]
    cmaq_fc_raw_pm25_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_raw_pm25_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cmaq_fc_raw_pm25_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cmaq_fc_raw_pm25_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_bm3_pm25_con_all_daily[:,:]
    cams_fc_bm3_pm25_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_pm25_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_bm3_pm25_con_all_daily[:,:]
    cmaq_fc_bm3_pm25_con_all_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_bm3_pm25_con_all_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_all_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_pm25_con_all_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)

    ## Regions 1-10, 99
    for rr in range(1,12):
        if rr == 1:
            inds_o3_reg = inds_o3_reg01_all
            inds_pm25_reg = inds_pm25_reg01_all
        elif rr == 2:
            inds_o3_reg = inds_o3_reg02_all
            inds_pm25_reg = inds_pm25_reg02_all
        elif rr == 3:
            inds_o3_reg = inds_o3_reg03_all
            inds_pm25_reg = inds_pm25_reg03_all
        elif rr == 4:
            inds_o3_reg = inds_o3_reg04_all
            inds_pm25_reg = inds_pm25_reg04_all
        elif rr == 5:
            inds_o3_reg = inds_o3_reg05_all
            inds_pm25_reg = inds_pm25_reg05_all
        elif rr == 6:
            inds_o3_reg = inds_o3_reg06_all
            inds_pm25_reg = inds_pm25_reg06_all
        elif rr == 7:
            inds_o3_reg = inds_o3_reg07_all
            inds_pm25_reg = inds_pm25_reg07_all
        elif rr == 8:
            inds_o3_reg = inds_o3_reg08_all
            inds_pm25_reg = inds_pm25_reg08_all
        elif rr == 9:
            inds_o3_reg = inds_o3_reg09_all
            inds_pm25_reg = inds_pm25_reg09_all
        elif rr == 10:
            inds_o3_reg = inds_o3_reg10_all
            inds_pm25_reg = inds_pm25_reg10_all
        elif rr == 11:
            inds_o3_reg = inds_o3_reg99_all
            inds_pm25_reg = inds_pm25_reg99_all

        obs = obs_o3_con_all_daily[:,inds_o3_reg]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            obs_o3_con_all_daily_mean[:, rr] = np.nanmean(obs, axis=axis)
            obs_o3_con_all_daily_sdev[:, rr] = np.nanstd( obs, axis=axis)
        fcst = cams_ra_raw_o3_con_all_daily[:,inds_o3_reg]
        cams_ra_raw_o3_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_ra_raw_o3_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_ra_raw_o3_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_ra_raw_o3_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_raw_o3_con_all_daily[:,inds_o3_reg]
        cams_fc_raw_o3_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_raw_o3_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_raw_o3_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_raw_o3_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_raw_o3_con_all_daily[:,inds_o3_reg]
        cmaq_fc_raw_o3_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_raw_o3_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_raw_o3_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_raw_o3_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_bm3_o3_con_all_daily[:,inds_o3_reg]
        cams_fc_bm3_o3_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_bm3_o3_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_bm3_o3_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_bm3_o3_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_bm3_o3_con_all_daily[:,inds_o3_reg]
        cmaq_fc_bm3_o3_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_bm3_o3_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_bm3_o3_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_bm3_o3_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)

        obs = obs_pm25_con_all_daily[:,inds_pm25_reg]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            obs_pm25_con_all_daily_mean[:, rr] = np.nanmean(obs, axis=axis)
            obs_pm25_con_all_daily_sdev[:, rr] = np.nanstd( obs, axis=axis)
        fcst = cams_ra_raw_pm25_con_all_daily[:,inds_pm25_reg]
        cams_ra_raw_pm25_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_ra_raw_pm25_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_ra_raw_pm25_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_ra_raw_pm25_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_raw_pm25_con_all_daily[:,inds_pm25_reg]
        cams_fc_raw_pm25_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_raw_pm25_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_raw_pm25_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_raw_pm25_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_raw_pm25_con_all_daily[:,inds_pm25_reg]
        cmaq_fc_raw_pm25_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_raw_pm25_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_raw_pm25_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_raw_pm25_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_bm3_pm25_con_all_daily[:,inds_pm25_reg]
        cams_fc_bm3_pm25_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_bm3_pm25_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_bm3_pm25_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_bm3_pm25_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_bm3_pm25_con_all_daily[:,inds_pm25_reg]
        cmaq_fc_bm3_pm25_con_all_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_bm3_pm25_con_all_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_bm3_pm25_con_all_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_bm3_pm25_con_all_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)

    ## Now do the same thing for the validation stations
    cams_ra_raw_o3_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_o3_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_o3_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_o3_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_o3_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_ra_raw_pm25_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_pm25_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_pm25_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_pm25_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_pm25_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)
    obs_o3_con_val_daily_mean   = np.full([n_valid_days, 12], np.nan)
    obs_pm25_con_val_daily_mean = np.full([n_valid_days, 12], np.nan)

    cams_ra_raw_o3_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_o3_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_o3_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_o3_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_o3_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_ra_raw_pm25_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_pm25_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_pm25_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_pm25_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_pm25_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)
    obs_o3_con_val_daily_sdev   = np.full([n_valid_days, 12], np.nan)
    obs_pm25_con_val_daily_sdev = np.full([n_valid_days, 12], np.nan)

    cams_ra_raw_o3_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_o3_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_o3_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_o3_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_o3_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_ra_raw_pm25_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_pm25_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_pm25_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_pm25_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_pm25_con_val_daily_rmse = np.full([n_valid_days, 12], np.nan)

    cams_ra_raw_o3_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_o3_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_o3_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_o3_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_o3_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_ra_raw_pm25_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_fc_raw_pm25_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_raw_pm25_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cams_fc_bm3_pm25_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)
    cmaq_fc_bm3_pm25_con_val_daily_mbe = np.full([n_valid_days, 12], np.nan)

    ## Overall
    axis = 1
    obs = obs_o3_con_val_daily[:,:]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        obs_o3_con_val_daily_mean[:, 0] = np.nanmean(obs, axis=axis)
        obs_o3_con_val_daily_sdev[:, 0] = np.nanstd( obs, axis=axis)
    fcst = cams_ra_raw_o3_con_val_daily[:,:]
    cams_ra_raw_o3_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_ra_raw_o3_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_ra_raw_o3_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_ra_raw_o3_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_raw_o3_con_val_daily[:,:]
    cams_fc_raw_o3_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_raw_o3_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_raw_o3_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_raw_o3_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_raw_o3_con_val_daily[:,:]
    cmaq_fc_raw_o3_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_raw_o3_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cmaq_fc_raw_o3_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cmaq_fc_raw_o3_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_bm3_o3_con_val_daily[:,:]
    cams_fc_bm3_o3_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_o3_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_bm3_o3_con_val_daily[:,:]
    cmaq_fc_bm3_o3_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_bm3_o3_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_o3_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_o3_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)

    obs = obs_pm25_con_val_daily[:,:]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        obs_pm25_con_val_daily_mean[:, 0] = np.nanmean(obs, axis=axis)
        obs_pm25_con_val_daily_sdev[:, 0] = np.nanstd( obs, axis=axis)
    fcst = cams_ra_raw_pm25_con_val_daily[:,:]
    cams_ra_raw_pm25_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_ra_raw_pm25_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_ra_raw_pm25_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_ra_raw_pm25_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_raw_pm25_con_val_daily[:,:]
    cams_fc_raw_pm25_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_raw_pm25_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_raw_pm25_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_raw_pm25_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_raw_pm25_con_val_daily[:,:]
    cmaq_fc_raw_pm25_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_raw_pm25_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cmaq_fc_raw_pm25_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cmaq_fc_raw_pm25_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cams_fc_bm3_pm25_con_val_daily[:,:]
    cams_fc_bm3_pm25_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_pm25_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)
    fcst = cmaq_fc_bm3_pm25_con_val_daily[:,:]
    cmaq_fc_bm3_pm25_con_val_daily_mbe[ :, 0] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
    cmaq_fc_bm3_pm25_con_val_daily_rmse[:, 0] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
    cams_fc_bm3_pm25_con_val_daily_mean[:, 0] = np.nanmean(fcst, axis=axis)
    cams_fc_bm3_pm25_con_val_daily_sdev[:, 0] = np.nanstd( fcst, axis=axis)

    ## Regions 1-10, 99
    for rr in range(1,12):
        if rr == 1:
            inds_o3_reg = inds_o3_reg01_val
            inds_pm25_reg = inds_pm25_reg01_val
        elif rr == 2:
            inds_o3_reg = inds_o3_reg02_val
            inds_pm25_reg = inds_pm25_reg02_val
        elif rr == 3:
            inds_o3_reg = inds_o3_reg03_val
            inds_pm25_reg = inds_pm25_reg03_val
        elif rr == 4:
            inds_o3_reg = inds_o3_reg04_val
            inds_pm25_reg = inds_pm25_reg04_val
        elif rr == 5:
            inds_o3_reg = inds_o3_reg05_val
            inds_pm25_reg = inds_pm25_reg05_val
        elif rr == 6:
            inds_o3_reg = inds_o3_reg06_val
            inds_pm25_reg = inds_pm25_reg06_val
        elif rr == 7:
            inds_o3_reg = inds_o3_reg07_val
            inds_pm25_reg = inds_pm25_reg07_val
        elif rr == 8:
            inds_o3_reg = inds_o3_reg08_val
            inds_pm25_reg = inds_pm25_reg08_val
        elif rr == 9:
            inds_o3_reg = inds_o3_reg09_val
            inds_pm25_reg = inds_pm25_reg09_val
        elif rr == 10:
            inds_o3_reg = inds_o3_reg10_val
            inds_pm25_reg = inds_pm25_reg10_val
        elif rr == 11:
            inds_o3_reg = inds_o3_reg99_val
            inds_pm25_reg = inds_pm25_reg99_val

        obs = obs_o3_con_val_daily[:,inds_o3_reg]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            obs_o3_con_val_daily_mean[:, rr] = np.nanmean(obs, axis=axis)
            obs_o3_con_val_daily_sdev[:, rr] = np.nanstd( obs, axis=axis)
        fcst = cams_ra_raw_o3_con_val_daily[:,inds_o3_reg]
        cams_ra_raw_o3_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_ra_raw_o3_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_ra_raw_o3_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_ra_raw_o3_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_raw_o3_con_val_daily[:,inds_o3_reg]
        cams_fc_raw_o3_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_raw_o3_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_raw_o3_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_raw_o3_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_raw_o3_con_val_daily[:,inds_o3_reg]
        cmaq_fc_raw_o3_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_raw_o3_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_raw_o3_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_raw_o3_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_bm3_o3_con_val_daily[:,inds_o3_reg]
        cams_fc_bm3_o3_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_bm3_o3_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_bm3_o3_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_bm3_o3_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_bm3_o3_con_val_daily[:,inds_o3_reg]
        cmaq_fc_bm3_o3_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_bm3_o3_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_bm3_o3_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_bm3_o3_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)

        obs = obs_pm25_con_val_daily[:,inds_pm25_reg]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            obs_pm25_con_val_daily_mean[:, rr] = np.nanmean(obs, axis=axis)
            obs_pm25_con_val_daily_sdev[:, rr] = np.nanstd( obs, axis=axis)
        fcst = cams_ra_raw_pm25_con_val_daily[:,inds_pm25_reg]
        cams_ra_raw_pm25_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_ra_raw_pm25_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_ra_raw_pm25_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_ra_raw_pm25_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_raw_pm25_con_val_daily[:,inds_pm25_reg]
        cams_fc_raw_pm25_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_raw_pm25_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_raw_pm25_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_raw_pm25_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_raw_pm25_con_val_daily[:,inds_pm25_reg]
        cmaq_fc_raw_pm25_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_raw_pm25_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_raw_pm25_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_raw_pm25_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cams_fc_bm3_pm25_con_val_daily[:,inds_pm25_reg]
        cams_fc_bm3_pm25_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cams_fc_bm3_pm25_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cams_fc_bm3_pm25_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cams_fc_bm3_pm25_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)
        fcst = cmaq_fc_bm3_pm25_con_val_daily[:,inds_pm25_reg]
        cmaq_fc_bm3_pm25_con_val_daily_mbe[ :, rr] = gen_funcs.calc_mbe( fcst, obs, axis=axis)
        cmaq_fc_bm3_pm25_con_val_daily_rmse[:, rr] = gen_funcs.calc_rmse(fcst, obs, axis=axis)
        cmaq_fc_bm3_pm25_con_val_daily_mean[:, rr] = np.nanmean(fcst, axis=axis)
        cmaq_fc_bm3_pm25_con_val_daily_sdev[:, rr] = np.nanstd( fcst, axis=axis)

    ## Calculate overall stats for each region
    ## Create arrays to hold the overall mean, sdev, rmse, mbe
    cams_ra_raw_o3_con_all_overall_mean = np.full([12], np.nan)
    cams_fc_raw_o3_con_all_overall_mean = np.full([12], np.nan)
    cmaq_fc_raw_o3_con_all_overall_mean = np.full([12], np.nan)
    cams_fc_bm3_o3_con_all_overall_mean = np.full([12], np.nan)
    cmaq_fc_bm3_o3_con_all_overall_mean = np.full([12], np.nan)
    cams_ra_raw_pm25_con_all_overall_mean = np.full([12], np.nan)
    cams_fc_raw_pm25_con_all_overall_mean = np.full([12], np.nan)
    cmaq_fc_raw_pm25_con_all_overall_mean = np.full([12], np.nan)
    cams_fc_bm3_pm25_con_all_overall_mean = np.full([12], np.nan)
    cmaq_fc_bm3_pm25_con_all_overall_mean = np.full([12], np.nan)
    obs_o3_con_all_overall_mean   = np.full([12], np.nan)
    obs_pm25_con_all_overall_mean = np.full([12], np.nan)

    cams_ra_raw_o3_con_all_overall_sdev = np.full([12], np.nan)
    cams_fc_raw_o3_con_all_overall_sdev = np.full([12], np.nan)
    cmaq_fc_raw_o3_con_all_overall_sdev = np.full([12], np.nan)
    cams_fc_bm3_o3_con_all_overall_sdev = np.full([12], np.nan)
    cmaq_fc_bm3_o3_con_all_overall_sdev = np.full([12], np.nan)
    cams_ra_raw_pm25_con_all_overall_sdev = np.full([12], np.nan)
    cams_fc_raw_pm25_con_all_overall_sdev = np.full([12], np.nan)
    cmaq_fc_raw_pm25_con_all_overall_sdev = np.full([12], np.nan)
    cams_fc_bm3_pm25_con_all_overall_sdev = np.full([12], np.nan)
    cmaq_fc_bm3_pm25_con_all_overall_sdev = np.full([12], np.nan)
    obs_o3_con_all_overall_sdev   = np.full([12], np.nan)
    obs_pm25_con_all_overall_sdev = np.full([12], np.nan)

    cams_ra_raw_o3_con_all_overall_rmse = np.full([12], np.nan)
    cams_fc_raw_o3_con_all_overall_rmse = np.full([12], np.nan)
    cmaq_fc_raw_o3_con_all_overall_rmse = np.full([12], np.nan)
    cams_fc_bm3_o3_con_all_overall_rmse = np.full([12], np.nan)
    cmaq_fc_bm3_o3_con_all_overall_rmse = np.full([12], np.nan)
    cams_ra_raw_pm25_con_all_overall_rmse = np.full([12], np.nan)
    cams_fc_raw_pm25_con_all_overall_rmse = np.full([12], np.nan)
    cmaq_fc_raw_pm25_con_all_overall_rmse = np.full([12], np.nan)
    cams_fc_bm3_pm25_con_all_overall_rmse = np.full([12], np.nan)
    cmaq_fc_bm3_pm25_con_all_overall_rmse = np.full([12], np.nan)

    cams_ra_raw_o3_con_all_overall_mbe = np.full([12], np.nan)
    cams_fc_raw_o3_con_all_overall_mbe = np.full([12], np.nan)
    cmaq_fc_raw_o3_con_all_overall_mbe = np.full([12], np.nan)
    cams_fc_bm3_o3_con_all_overall_mbe = np.full([12], np.nan)
    cmaq_fc_bm3_o3_con_all_overall_mbe = np.full([12], np.nan)
    cams_ra_raw_pm25_con_all_overall_mbe = np.full([12], np.nan)
    cams_fc_raw_pm25_con_all_overall_mbe = np.full([12], np.nan)
    cmaq_fc_raw_pm25_con_all_overall_mbe = np.full([12], np.nan)
    cams_fc_bm3_pm25_con_all_overall_mbe = np.full([12], np.nan)
    cmaq_fc_bm3_pm25_con_all_overall_mbe = np.full([12], np.nan)

    ## CONUS (all regions)
    obs = obs_o3_con_all[:,:]
    obs_o3_con_all_overall_mean[0] = np.nanmean(obs)
    obs_o3_con_all_overall_sdev[0] = np.nanstd( obs)
    fcst = cams_ra_raw_o3_con_all[:,:]
    cams_ra_raw_o3_con_all_overall_mean[0] = np.nanmean(fcst)
    cams_ra_raw_o3_con_all_overall_sdev[0] = np.nanstd( fcst)
    cams_ra_raw_o3_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_ra_raw_o3_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cams_fc_raw_o3_con_all[:,:]
    cams_fc_raw_o3_con_all_overall_mean[0] = np.nanmean(fcst)
    cams_fc_raw_o3_con_all_overall_sdev[0] = np.nanstd( fcst)
    cams_fc_raw_o3_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_fc_raw_o3_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cmaq_fc_raw_o3_con_all[:,:]
    cmaq_fc_raw_o3_con_all_overall_mean[0] = np.nanmean(fcst)
    cmaq_fc_raw_o3_con_all_overall_sdev[0] = np.nanstd( fcst)
    cmaq_fc_raw_o3_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cmaq_fc_raw_o3_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cams_fc_bm3_o3_con_all[:,:]
    cams_fc_bm3_o3_con_all_overall_mean[0] = np.nanmean(fcst)
    cams_fc_bm3_o3_con_all_overall_sdev[0] = np.nanstd( fcst)
    cams_fc_bm3_o3_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_fc_bm3_o3_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cmaq_fc_bm3_o3_con_all[:,:]
    cmaq_fc_bm3_o3_con_all_overall_mean[0] = np.nanmean(fcst)
    cmaq_fc_bm3_o3_con_all_overall_sdev[0] = np.nanstd( fcst)
    cmaq_fc_bm3_o3_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cmaq_fc_bm3_o3_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)

    obs = obs_pm25_con_all[:,:]
    obs_pm25_con_all_overall_mean[0] = np.nanmean(obs)
    obs_pm25_con_all_overall_sdev[0] = np.nanstd( obs)
    fcst = cams_ra_raw_pm25_con_all[:,:]
    cams_ra_raw_pm25_con_all_overall_mean[0] = np.nanmean(fcst)
    cams_ra_raw_pm25_con_all_overall_sdev[0] = np.nanstd( fcst)
    cams_ra_raw_pm25_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_ra_raw_pm25_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cams_fc_raw_pm25_con_all[:,:]
    cams_fc_raw_pm25_con_all_overall_mean[0] = np.nanmean(fcst)
    cams_fc_raw_pm25_con_all_overall_sdev[0] = np.nanstd( fcst)
    cams_fc_raw_pm25_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_fc_raw_pm25_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cmaq_fc_raw_pm25_con_all[:,:]
    cmaq_fc_raw_pm25_con_all_overall_mean[0] = np.nanmean(fcst)
    cmaq_fc_raw_pm25_con_all_overall_sdev[0] = np.nanstd( fcst)
    cmaq_fc_raw_pm25_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cmaq_fc_raw_pm25_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cams_fc_bm3_pm25_con_all[:,:]
    cams_fc_bm3_pm25_con_all_overall_mean[0] = np.nanmean(fcst)
    cams_fc_bm3_pm25_con_all_overall_sdev[0] = np.nanstd( fcst)
    cams_fc_bm3_pm25_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_fc_bm3_pm25_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cmaq_fc_bm3_pm25_con_all[:,:]
    cmaq_fc_bm3_pm25_con_all_overall_mean[0] = np.nanmean(fcst)
    cmaq_fc_bm3_pm25_con_all_overall_sdev[0] = np.nanstd( fcst)
    cmaq_fc_bm3_pm25_con_all_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cmaq_fc_bm3_pm25_con_all_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)

    ## Regions 1-10, 99
    for rr in range(1,12):
        if rr == 1:
            inds_o3_reg = inds_o3_reg01_all
            inds_pm25_reg = inds_pm25_reg01_all
        elif rr == 2:
            inds_o3_reg = inds_o3_reg02_all
            inds_pm25_reg = inds_pm25_reg02_all
        elif rr == 3:
            inds_o3_reg = inds_o3_reg03_all
            inds_pm25_reg = inds_pm25_reg03_all
        elif rr == 4:
            inds_o3_reg = inds_o3_reg04_all
            inds_pm25_reg = inds_pm25_reg04_all
        elif rr == 5:
            inds_o3_reg = inds_o3_reg05_all
            inds_pm25_reg = inds_pm25_reg05_all
        elif rr == 6:
            inds_o3_reg = inds_o3_reg06_all
            inds_pm25_reg = inds_pm25_reg06_all
        elif rr == 7:
            inds_o3_reg = inds_o3_reg07_all
            inds_pm25_reg = inds_pm25_reg07_all
        elif rr == 8:
            inds_o3_reg = inds_o3_reg08_all
            inds_pm25_reg = inds_pm25_reg08_all
        elif rr == 9:
            inds_o3_reg = inds_o3_reg09_all
            inds_pm25_reg = inds_pm25_reg09_all
        elif rr == 10:
            inds_o3_reg = inds_o3_reg10_all
            inds_pm25_reg = inds_pm25_reg10_all
        elif rr == 11:
            inds_o3_reg = inds_o3_reg99_all
            inds_pm25_reg = inds_pm25_reg99_all

        obs = obs_o3_con_all[:,inds_o3_reg]
        obs_o3_con_all_overall_mean[rr] = np.nanmean(obs)
        obs_o3_con_all_overall_sdev[rr] = np.nanstd( obs)
        fcst = cams_ra_raw_o3_con_all[:,inds_o3_reg]
        cams_ra_raw_o3_con_all_overall_mean[rr] = np.nanmean(fcst)
        cams_ra_raw_o3_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cams_ra_raw_o3_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_ra_raw_o3_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cams_fc_raw_o3_con_all[:,inds_o3_reg]
        cams_fc_raw_o3_con_all_overall_mean[rr] = np.nanmean(fcst)
        cams_fc_raw_o3_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cams_fc_raw_o3_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_fc_raw_o3_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cmaq_fc_raw_o3_con_all[:,inds_o3_reg]
        cmaq_fc_raw_o3_con_all_overall_mean[rr] = np.nanmean(fcst)
        cmaq_fc_raw_o3_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cmaq_fc_raw_o3_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cmaq_fc_raw_o3_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cams_fc_bm3_o3_con_all[:,inds_o3_reg]
        cams_fc_bm3_o3_con_all_overall_mean[rr] = np.nanmean(fcst)
        cams_fc_bm3_o3_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cams_fc_bm3_o3_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_fc_bm3_o3_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cmaq_fc_bm3_o3_con_all[:,inds_o3_reg]
        cmaq_fc_bm3_o3_con_all_overall_mean[rr] = np.nanmean(fcst)
        cmaq_fc_bm3_o3_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cmaq_fc_bm3_o3_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cmaq_fc_bm3_o3_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)

        obs = obs_pm25_con_all[:,inds_pm25_reg]
        obs_pm25_con_all_overall_mean[rr] = np.nanmean(obs)
        obs_pm25_con_all_overall_sdev[rr] = np.nanstd( obs)
        fcst = cams_ra_raw_pm25_con_all[:,inds_pm25_reg]
        cams_ra_raw_pm25_con_all_overall_mean[rr] = np.nanmean(fcst)
        cams_ra_raw_pm25_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cams_ra_raw_pm25_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_ra_raw_pm25_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cams_fc_raw_pm25_con_all[:,inds_pm25_reg]
        cams_fc_raw_pm25_con_all_overall_mean[rr] = np.nanmean(fcst)
        cams_fc_raw_pm25_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cams_fc_raw_pm25_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_fc_raw_pm25_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cmaq_fc_raw_pm25_con_all[:,inds_pm25_reg]
        cmaq_fc_raw_pm25_con_all_overall_mean[rr] = np.nanmean(fcst)
        cmaq_fc_raw_pm25_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cmaq_fc_raw_pm25_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cmaq_fc_raw_pm25_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cams_fc_bm3_pm25_con_all[:,inds_pm25_reg]
        cams_fc_bm3_pm25_con_all_overall_mean[rr] = np.nanmean(fcst)
        cams_fc_bm3_pm25_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cams_fc_bm3_pm25_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_fc_bm3_pm25_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cmaq_fc_bm3_pm25_con_all[:,inds_pm25_reg]
        cmaq_fc_bm3_pm25_con_all_overall_mean[rr] = np.nanmean(fcst)
        cmaq_fc_bm3_pm25_con_all_overall_sdev[rr] = np.nanstd( fcst)
        cmaq_fc_bm3_pm25_con_all_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cmaq_fc_bm3_pm25_con_all_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)

    ## Now do the same for the validation stations
    cams_ra_raw_o3_con_val_overall_mean = np.full([12], np.nan)
    cams_fc_raw_o3_con_val_overall_mean = np.full([12], np.nan)
    cmaq_fc_raw_o3_con_val_overall_mean = np.full([12], np.nan)
    cams_fc_bm3_o3_con_val_overall_mean = np.full([12], np.nan)
    cmaq_fc_bm3_o3_con_val_overall_mean = np.full([12], np.nan)
    cams_ra_raw_pm25_con_val_overall_mean = np.full([12], np.nan)
    cams_fc_raw_pm25_con_val_overall_mean = np.full([12], np.nan)
    cmaq_fc_raw_pm25_con_val_overall_mean = np.full([12], np.nan)
    cams_fc_bm3_pm25_con_val_overall_mean = np.full([12], np.nan)
    cmaq_fc_bm3_pm25_con_val_overall_mean = np.full([12], np.nan)
    obs_o3_con_val_overall_mean   = np.full([12], np.nan)
    obs_pm25_con_val_overall_mean = np.full([12], np.nan)

    cams_ra_raw_o3_con_val_overall_sdev = np.full([12], np.nan)
    cams_fc_raw_o3_con_val_overall_sdev = np.full([12], np.nan)
    cmaq_fc_raw_o3_con_val_overall_sdev = np.full([12], np.nan)
    cams_fc_bm3_o3_con_val_overall_sdev = np.full([12], np.nan)
    cmaq_fc_bm3_o3_con_val_overall_sdev = np.full([12], np.nan)
    cams_ra_raw_pm25_con_val_overall_sdev = np.full([12], np.nan)
    cams_fc_raw_pm25_con_val_overall_sdev = np.full([12], np.nan)
    cmaq_fc_raw_pm25_con_val_overall_sdev = np.full([12], np.nan)
    cams_fc_bm3_pm25_con_val_overall_sdev = np.full([12], np.nan)
    cmaq_fc_bm3_pm25_con_val_overall_sdev = np.full([12], np.nan)
    obs_o3_con_val_overall_sdev   = np.full([12], np.nan)
    obs_pm25_con_val_overall_sdev = np.full([12], np.nan)

    cams_ra_raw_o3_con_val_overall_rmse = np.full([12], np.nan)
    cams_fc_raw_o3_con_val_overall_rmse = np.full([12], np.nan)
    cmaq_fc_raw_o3_con_val_overall_rmse = np.full([12], np.nan)
    cams_fc_bm3_o3_con_val_overall_rmse = np.full([12], np.nan)
    cmaq_fc_bm3_o3_con_val_overall_rmse = np.full([12], np.nan)
    cams_ra_raw_pm25_con_val_overall_rmse = np.full([12], np.nan)
    cams_fc_raw_pm25_con_val_overall_rmse = np.full([12], np.nan)
    cmaq_fc_raw_pm25_con_val_overall_rmse = np.full([12], np.nan)
    cams_fc_bm3_pm25_con_val_overall_rmse = np.full([12], np.nan)
    cmaq_fc_bm3_pm25_con_val_overall_rmse = np.full([12], np.nan)

    cams_ra_raw_o3_con_val_overall_mbe = np.full([12], np.nan)
    cams_fc_raw_o3_con_val_overall_mbe = np.full([12], np.nan)
    cmaq_fc_raw_o3_con_val_overall_mbe = np.full([12], np.nan)
    cams_fc_bm3_o3_con_val_overall_mbe = np.full([12], np.nan)
    cmaq_fc_bm3_o3_con_val_overall_mbe = np.full([12], np.nan)
    cams_ra_raw_pm25_con_val_overall_mbe = np.full([12], np.nan)
    cams_fc_raw_pm25_con_val_overall_mbe = np.full([12], np.nan)
    cmaq_fc_raw_pm25_con_val_overall_mbe = np.full([12], np.nan)
    cams_fc_bm3_pm25_con_val_overall_mbe = np.full([12], np.nan)
    cmaq_fc_bm3_pm25_con_val_overall_mbe = np.full([12], np.nan)

    ## CONUS (all regions)
    obs = obs_o3_con_val[:,:]
    obs_o3_con_val_overall_mean[0] = np.nanmean(obs)
    obs_o3_con_val_overall_sdev[0] = np.nanstd( obs)
    fcst = cams_ra_raw_o3_con_val[:,:]
    cams_ra_raw_o3_con_val_overall_mean[0] = np.nanmean(fcst)
    cams_ra_raw_o3_con_val_overall_sdev[0] = np.nanstd( fcst)
    cams_ra_raw_o3_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_ra_raw_o3_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cams_fc_raw_o3_con_val[:,:]
    cams_fc_raw_o3_con_val_overall_mean[0] = np.nanmean(fcst)
    cams_fc_raw_o3_con_val_overall_sdev[0] = np.nanstd( fcst)
    cams_fc_raw_o3_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_fc_raw_o3_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cmaq_fc_raw_o3_con_val[:,:]
    cmaq_fc_raw_o3_con_val_overall_mean[0] = np.nanmean(fcst)
    cmaq_fc_raw_o3_con_val_overall_sdev[0] = np.nanstd( fcst)
    cmaq_fc_raw_o3_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cmaq_fc_raw_o3_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cams_fc_bm3_o3_con_val[:,:]
    cams_fc_bm3_o3_con_val_overall_mean[0] = np.nanmean(fcst)
    cams_fc_bm3_o3_con_val_overall_sdev[0] = np.nanstd( fcst)
    cams_fc_bm3_o3_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_fc_bm3_o3_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cmaq_fc_bm3_o3_con_val[:,:]
    cmaq_fc_bm3_o3_con_val_overall_mean[0] = np.nanmean(fcst)
    cmaq_fc_bm3_o3_con_val_overall_sdev[0] = np.nanstd( fcst)
    cmaq_fc_bm3_o3_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cmaq_fc_bm3_o3_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)

    obs = obs_pm25_con_val[:,:]
    obs_pm25_con_val_overall_mean[0] = np.nanmean(obs)
    obs_pm25_con_val_overall_sdev[0] = np.nanstd( obs)
    fcst = cams_ra_raw_pm25_con_val[:,:]
    cams_ra_raw_pm25_con_val_overall_mean[0] = np.nanmean(fcst)
    cams_ra_raw_pm25_con_val_overall_sdev[0] = np.nanstd( fcst)
    cams_ra_raw_pm25_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_ra_raw_pm25_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cams_fc_raw_pm25_con_val[:,:]
    cams_fc_raw_pm25_con_val_overall_mean[0] = np.nanmean(fcst)
    cams_fc_raw_pm25_con_val_overall_sdev[0] = np.nanstd( fcst)
    cams_fc_raw_pm25_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_fc_raw_pm25_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cmaq_fc_raw_pm25_con_val[:,:]
    cmaq_fc_raw_pm25_con_val_overall_mean[0] = np.nanmean(fcst)
    cmaq_fc_raw_pm25_con_val_overall_sdev[0] = np.nanstd( fcst)
    cmaq_fc_raw_pm25_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cmaq_fc_raw_pm25_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cams_fc_bm3_pm25_con_val[:,:]
    cams_fc_bm3_pm25_con_val_overall_mean[0] = np.nanmean(fcst)
    cams_fc_bm3_pm25_con_val_overall_sdev[0] = np.nanstd( fcst)
    cams_fc_bm3_pm25_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cams_fc_bm3_pm25_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)
    fcst = cmaq_fc_bm3_pm25_con_val[:,:]
    cmaq_fc_bm3_pm25_con_val_overall_mean[0] = np.nanmean(fcst)
    cmaq_fc_bm3_pm25_con_val_overall_sdev[0] = np.nanstd( fcst)
    cmaq_fc_bm3_pm25_con_val_overall_rmse[0] = gen_funcs.calc_rmse(fcst, obs)
    cmaq_fc_bm3_pm25_con_val_overall_mbe[ 0] = gen_funcs.calc_mbe( fcst, obs)

    ## Regions 1-10, 99
    for rr in range(1,12):
        if rr == 1:
            inds_o3_reg = inds_o3_reg01_val
            inds_pm25_reg = inds_pm25_reg01_val
        elif rr == 2:
            inds_o3_reg = inds_o3_reg02_val
            inds_pm25_reg = inds_pm25_reg02_val
        elif rr == 3:
            inds_o3_reg = inds_o3_reg03_val
            inds_pm25_reg = inds_pm25_reg03_val
        elif rr == 4:
            inds_o3_reg = inds_o3_reg04_val
            inds_pm25_reg = inds_pm25_reg04_val
        elif rr == 5:
            inds_o3_reg = inds_o3_reg05_val
            inds_pm25_reg = inds_pm25_reg05_val
        elif rr == 6:
            inds_o3_reg = inds_o3_reg06_val
            inds_pm25_reg = inds_pm25_reg06_val
        elif rr == 7:
            inds_o3_reg = inds_o3_reg07_val
            inds_pm25_reg = inds_pm25_reg07_val
        elif rr == 8:
            inds_o3_reg = inds_o3_reg08_val
            inds_pm25_reg = inds_pm25_reg08_val
        elif rr == 9:
            inds_o3_reg = inds_o3_reg09_val
            inds_pm25_reg = inds_pm25_reg09_val
        elif rr == 10:
            inds_o3_reg = inds_o3_reg10_val
            inds_pm25_reg = inds_pm25_reg10_val
        elif rr == 11:
            inds_o3_reg = inds_o3_reg99_val
            inds_pm25_reg = inds_pm25_reg99_val

        obs = obs_o3_con_val[:,inds_o3_reg]
        obs_o3_con_val_overall_mean[rr] = np.nanmean(obs)
        obs_o3_con_val_overall_sdev[rr] = np.nanstd( obs)
        fcst = cams_ra_raw_o3_con_val[:,inds_o3_reg]
        cams_ra_raw_o3_con_val_overall_mean[rr] = np.nanmean(fcst)
        cams_ra_raw_o3_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cams_ra_raw_o3_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_ra_raw_o3_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cams_fc_raw_o3_con_val[:,inds_o3_reg]
        cams_fc_raw_o3_con_val_overall_mean[rr] = np.nanmean(fcst)
        cams_fc_raw_o3_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cams_fc_raw_o3_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_fc_raw_o3_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cmaq_fc_raw_o3_con_val[:,inds_o3_reg]
        cmaq_fc_raw_o3_con_val_overall_mean[rr] = np.nanmean(fcst)
        cmaq_fc_raw_o3_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cmaq_fc_raw_o3_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cmaq_fc_raw_o3_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cams_fc_bm3_o3_con_val[:,inds_o3_reg]
        cams_fc_bm3_o3_con_val_overall_mean[rr] = np.nanmean(fcst)
        cams_fc_bm3_o3_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cams_fc_bm3_o3_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_fc_bm3_o3_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cmaq_fc_bm3_o3_con_val[:,inds_o3_reg]
        cmaq_fc_bm3_o3_con_val_overall_mean[rr] = np.nanmean(fcst)
        cmaq_fc_bm3_o3_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cmaq_fc_bm3_o3_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cmaq_fc_bm3_o3_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)

        obs = obs_pm25_con_val[:,inds_pm25_reg]
        obs_pm25_con_val_overall_mean[rr] = np.nanmean(obs)
        obs_pm25_con_val_overall_sdev[rr] = np.nanstd( obs)
        fcst = cams_ra_raw_pm25_con_val[:,inds_pm25_reg]
        cams_ra_raw_pm25_con_val_overall_mean[rr] = np.nanmean(fcst)
        cams_ra_raw_pm25_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cams_ra_raw_pm25_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_ra_raw_pm25_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cams_fc_raw_pm25_con_val[:,inds_pm25_reg]
        cams_fc_raw_pm25_con_val_overall_mean[rr] = np.nanmean(fcst)
        cams_fc_raw_pm25_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cams_fc_raw_pm25_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_fc_raw_pm25_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cmaq_fc_raw_pm25_con_val[:,inds_pm25_reg]
        cmaq_fc_raw_pm25_con_val_overall_mean[rr] = np.nanmean(fcst)
        cmaq_fc_raw_pm25_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cmaq_fc_raw_pm25_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cmaq_fc_raw_pm25_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cams_fc_bm3_pm25_con_val[:,inds_pm25_reg]
        cams_fc_bm3_pm25_con_val_overall_mean[rr] = np.nanmean(fcst)
        cams_fc_bm3_pm25_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cams_fc_bm3_pm25_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cams_fc_bm3_pm25_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)
        fcst = cmaq_fc_bm3_pm25_con_val[:,inds_pm25_reg]
        cmaq_fc_bm3_pm25_con_val_overall_mean[rr] = np.nanmean(fcst)
        cmaq_fc_bm3_pm25_con_val_overall_sdev[rr] = np.nanstd( fcst)
        cmaq_fc_bm3_pm25_con_val_overall_rmse[rr] = gen_funcs.calc_rmse(fcst, obs)
        cmaq_fc_bm3_pm25_con_val_overall_mbe[ rr] = gen_funcs.calc_mbe( fcst, obs)

    ## ****************
    ## PLOTTING SECTION
    ## ****************

    colors_raw_bm3  = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'black']
    dashes_raw_bm3  = ['--','--','--','-','-','-']
    markers_raw_bm3 = ['v','^','<','>','s','o']
    lglab_raw_bm3_obs = ['CAMS RA Raw', 'CAMS FC Raw', 'CMAQ FC Raw', 'CAMS FC BC', 'CMAQ FC BC', 'AirNow Obs']
    lglab_raw_bm3 = ['CAMS RA Raw', 'CAMS FC Raw', 'CMAQ FC Raw', 'CAMS FC BC', 'CMAQ FC BC']

    colors_raw  = ['tab:blue', 'tab:orange', 'tab:green', 'black']
    dashes_raw  = ['--','--','--','-']
    markers_raw = ['v','^','<','o']
    lglab_raw_obs = ['CAMS RA Raw', 'CAMS FC Raw', 'CMAQ FC Raw', 'AirNow Obs']
    lglab_raw = ['CAMS RA Raw', 'CAMS FC Raw', 'CMAQ FC Raw']

    regions = ['CONUS', 'Region 1', 'Region 2', 'Region 3', 'Region 4', 'Region 5',
                'Region 6', 'Region 7', 'Region 8', 'Region 9', 'Region 10', 'Non-U.S.']

    xarr_diurnal = np.arange(0,24)
    xarr_daily = dt_dates_valid
    xlab_diurnal = 'Hour of Day [UTC]\n\n '
    xlab_daily = 'Date'
    ylab_o3_conc = mpl_o3+' Concentration [ppbv]'
    ylab_o3_mean = mpl_o3+' Mean Concentration [ppbv]'
    ylab_o3_rmse = mpl_o3+' RMSE [ppbv]'
    ylab_o3_mbe  = mpl_o3+' MBE [ppbv]'
    ylab_pm25_conc = mpl_pm25+' Concentration ['+mpl_ugm3+']'
    ylab_pm25_mean = mpl_pm25+' Mean Concentration ['+mpl_ugm3+']'
    ylab_pm25_rmse = mpl_pm25+' RMSE ['+mpl_ugm3+']'
    ylab_pm25_mbe  = mpl_pm25+' MBE ['+mpl_ugm3+']'

    ymin_o3_mean = 0.0
    ymax_o3_mean = 60.0
    ymin_o3_rmse = 0.0
    ymax_o3_rmse = 25.0
    ymin_o3_mbe  = -5.0
    ymax_o3_mbe  = 20.0
    ymin_o3_box  = 0.5

    ymin_pm25_mean = 0.0
    ymax_pm25_mean = 40.0
    ymin_pm25_rmse = 0.0
    ymax_pm25_rmse = 30.0
    ymin_pm25_mbe  = -20.0
    ymax_pm25_mbe  = 20.0
    ymin_pm25_box  = 0.5

    title1 = 'Concentration'
    title2 = 'Root Mean Squared Error'
    title3 = 'Mean Bias Error'

    ## Make a multi-panel plot of the overall data, both using all and validation stations
    ## Left/top panel: box plot of O3/PM2.5 distributions;
    ## Middle panel: Bar chart or RMSE;
    ## Right/bottom panel: Bar chart of MBE
    xlab1_raw_bm3 = lglab_raw_bm3_obs
    xlab2_raw_bm3 = lglab_raw_bm3
    xlab3_raw_bm3 = lglab_raw_bm3
    xlab1_raw = lglab_raw_obs
    xlab2_raw = lglab_raw
    xlab3_raw = lglab_raw

    title1_box = 'Distribution of Values'
    yscale_box_pm25 = 'log' # needed because of how large the outliers are compared to the IQR
    yscale_box_o3 = 'log'

    ## Make overall stats box/bar plots with all regions in one plot
    ## Start by making the plots against all stations and include the bias-corrected forecasts
    fname = plot_dir.joinpath('overall_o3_allreg_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle = 'Overall Stats By Region for '+mpl_o3+' at All Stations'

    data1_vals = np.full([n_valid_times*n_obs_o3_all, 4], np.nan)
    data1_vals[:,0] = np.ravel(cams_ra_raw_o3_con_all[:,:])
    data1_vals[:,1] = np.ravel(cams_fc_raw_o3_con_all[:,:])
    data1_vals[:,2] = np.ravel(cmaq_fc_raw_o3_con_all[:,:])
    data1_vals[:,3] = np.ravel(obs_o3_con_all[:,:])
    regions1 = np.full([n_valid_times*n_obs_o3_all, 4], regions[0], dtype=object)
    models1 = np.full([n_valid_times*n_obs_o3_all, 4], '', dtype=object)
    models1[:,0] = lglab_raw_obs[0]
    models1[:,1] = lglab_raw_obs[1]
    models1[:,2] = lglab_raw_obs[2]
    models1[:,3] = lglab_raw_obs[3]
    ## To conserve memory usage, and to ensure that only distributions of obs-model pairs are used,
    ## find rows where any value is NaN and drop that entire row. Obs will have more NaNs than models.
    for mm in range(4):
        if np.isnan(data1_vals[:,mm]).any():
            inds_nan = np.where(np.isnan(data1_vals[:,mm]))[0]
            data1_vals = np.delete(data1_vals, inds_nan, axis=0)
            regions1 = np.delete(regions1, inds_nan, axis=0)
            models1 = np.delete(models1, inds_nan, axis=0)
            
    data1_vals = np.ravel(data1_vals)
    regions1 = np.ravel(regions1)
    models1 = np.ravel(models1)
    ## Loop over the regions to keep building the very long dataframe
    for rr in range(1,12):
        if rr == 1:
            inds = inds_o3_reg01_all
        elif rr == 2:
            inds = inds_o3_reg02_all
        elif rr == 3:
            inds = inds_o3_reg03_all
        elif rr == 4:
            inds = inds_o3_reg04_all
        elif rr == 5:
            inds = inds_o3_reg05_all
        elif rr == 6:
            inds = inds_o3_reg06_all
        elif rr == 7:
            inds = inds_o3_reg07_all
        elif rr == 8:
            inds = inds_o3_reg08_all
        elif rr == 9:
            inds = inds_o3_reg09_all
        elif rr == 10:
            inds = inds_o3_reg10_all
        elif rr == 11:
            inds = inds_o3_reg99_all
        n_obs = len(inds)
        data1_vals_this = np.full([n_valid_times*n_obs, 4], np.nan)
        data1_vals_this[:,0] = np.ravel(cams_ra_raw_o3_con_all[:,inds])
        data1_vals_this[:,1] = np.ravel(cams_fc_raw_o3_con_all[:,inds])
        data1_vals_this[:,2] = np.ravel(cmaq_fc_raw_o3_con_all[:,inds])
        data1_vals_this[:,3] = np.ravel(obs_o3_con_all[:,inds])
        regions1_this = np.full([n_valid_times*n_obs, 4], regions[rr], dtype=object)
        models1_this = np.full([n_valid_times*n_obs, 4], '', dtype=object)
        models1_this[:,0] = lglab_raw_obs[0]
        models1_this[:,1] = lglab_raw_obs[1]
        models1_this[:,2] = lglab_raw_obs[2]
        models1_this[:,3] = lglab_raw_obs[3]

        ## To conserve memory usage, and to ensure that only distributions of valid obs-model pairs are used,
        ## find rows where any value is NaN and drop that entire row. Obs will have more NaNs than models.
        for mm in range(4):
            if np.isnan(data1_vals_this[:,mm]).any():
                inds_nan = np.where(np.isnan(data1_vals_this[:,mm]))[0]
                data1_vals_this = np.delete(data1_vals_this, inds_nan, axis=0)
                regions1_this = np.delete(regions1_this, inds_nan, axis=0)
                models1_this = np.delete(models1_this, inds_nan, axis=0)

        data1_vals_this = np.ravel(data1_vals_this)
        regions1_this = np.ravel(regions1_this)
        models1_this = np.ravel(models1_this)

        data1_vals = np.append(data1_vals, data1_vals_this)
        regions1 = np.append(regions1, regions1_this)
        models1 = np.append(models1, models1_this)

    data1 = pd.DataFrame({
                'Region': regions1,
                'Value': data1_vals,
                'Model': models1})

    data2_vals = np.full([3*12], np.nan)
    regions2 = np.full([3*12], '', dtype=object)
    models2 = np.full([3*12], '', dtype=object)
    for xx in range(12):
        data2_vals[xx*3+0] = cams_ra_raw_o3_con_all_overall_rmse[xx]
        data2_vals[xx*3+1] = cams_fc_raw_o3_con_all_overall_rmse[xx]
        data2_vals[xx*3+2] = cmaq_fc_raw_o3_con_all_overall_rmse[xx]
        regions2[xx*3+0:xx*3+3] = regions[xx]
        models2[xx*3+0] = lglab_raw[0]
        models2[xx*3+1] = lglab_raw[1]
        models2[xx*3+2] = lglab_raw[2]
    data2 = pd.DataFrame({
                'Value': data2_vals,
                'Region': regions2,
                'Model': models2})

    data3_vals = np.full([3*12], np.nan)
    regions3 = regions2
    models3 = models2
    for xx in range(12):
        data3_vals[xx*3+0] = cams_ra_raw_o3_con_all_overall_mbe[xx]
        data3_vals[xx*3+1] = cams_fc_raw_o3_con_all_overall_mbe[xx]
        data3_vals[xx*3+2] = cmaq_fc_raw_o3_con_all_overall_mbe[xx]
    data3 = pd.DataFrame({
                'Value': data3_vals,
                'Region': regions3,
                'Model': models3})
    plot_stat_box_bar_3panel_vert(fname, data1, data2, data3, regions, xlab1_raw, xlab2_raw, xlab3_raw,
        ylab_o3_conc, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1_box, title2, title3, colors_raw,
        ymin_box=ymin_o3_box, yscale_box=yscale_box_o3)


    fname = plot_dir.joinpath('overall_pm25_allreg_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle = 'Overall Stats By Region for '+mpl_pm25+' at All Stations'

    data1_vals = np.full([n_valid_times*n_obs_pm25_all, 4], np.nan)
    data1_vals[:,0] = np.ravel(cams_ra_raw_pm25_con_all[:,:])
    data1_vals[:,1] = np.ravel(cams_fc_raw_pm25_con_all[:,:])
    data1_vals[:,2] = np.ravel(cmaq_fc_raw_pm25_con_all[:,:])
    data1_vals[:,3] = np.ravel(obs_pm25_con_all[:,:])
    regions1 = np.full([n_valid_times*n_obs_pm25_all, 4], regions[0], dtype=object)
    models1 = np.full([n_valid_times*n_obs_pm25_all, 4], '', dtype=object)
    models1[:,0] = lglab_raw_obs[0]
    models1[:,1] = lglab_raw_obs[1]
    models1[:,2] = lglab_raw_obs[2]
    models1[:,3] = lglab_raw_obs[3]
    ## To conserve memory usage, and to ensure that only distributions of obs-model pairs are used,
    ## find rows where any value is NaN and drop that entire row. Obs will have more NaNs than models.
    for mm in range(4):
        if np.isnan(data1_vals[:,mm]).any():
            inds_nan = np.where(np.isnan(data1_vals[:,mm]))[0]
            data1_vals = np.delete(data1_vals, inds_nan, axis=0)
            regions1 = np.delete(regions1, inds_nan, axis=0)
            models1 = np.delete(models1, inds_nan, axis=0)
            
    data1_vals = np.ravel(data1_vals)
    regions1 = np.ravel(regions1)
    models1 = np.ravel(models1)
    ## Loop over the regions to keep building the very long dataframe
    for rr in range(1,12):
        if rr == 1:
            inds = inds_pm25_reg01_all
        elif rr == 2:
            inds = inds_pm25_reg02_all
        elif rr == 3:
            inds = inds_pm25_reg03_all
        elif rr == 4:
            inds = inds_pm25_reg04_all
        elif rr == 5:
            inds = inds_pm25_reg05_all
        elif rr == 6:
            inds = inds_pm25_reg06_all
        elif rr == 7:
            inds = inds_pm25_reg07_all
        elif rr == 8:
            inds = inds_pm25_reg08_all
        elif rr == 9:
            inds = inds_pm25_reg09_all
        elif rr == 10:
            inds = inds_pm25_reg10_all
        elif rr == 11:
            inds = inds_pm25_reg99_all
        n_obs = len(inds)
        data1_vals_this = np.full([n_valid_times*n_obs, 4], np.nan)
        data1_vals_this[:,0] = np.ravel(cams_ra_raw_pm25_con_all[:,inds])
        data1_vals_this[:,1] = np.ravel(cams_fc_raw_pm25_con_all[:,inds])
        data1_vals_this[:,2] = np.ravel(cmaq_fc_raw_pm25_con_all[:,inds])
        data1_vals_this[:,3] = np.ravel(obs_pm25_con_all[:,inds])
        regions1_this = np.full([n_valid_times*n_obs, 4], regions[rr], dtype=object)
        models1_this = np.full([n_valid_times*n_obs, 4], '', dtype=object)
        models1_this[:,0] = lglab_raw_obs[0]
        models1_this[:,1] = lglab_raw_obs[1]
        models1_this[:,2] = lglab_raw_obs[2]
        models1_this[:,3] = lglab_raw_obs[3]

        ## To conserve memory usage, and to ensure that only distributions of valid obs-model pairs are used,
        ## find rows where any value is NaN and drop that entire row. Obs will have more NaNs than models.
        for mm in range(4):
            if np.isnan(data1_vals_this[:,mm]).any():
                inds_nan = np.where(np.isnan(data1_vals_this[:,mm]))[0]
                data1_vals_this = np.delete(data1_vals_this, inds_nan, axis=0)
                regions1_this = np.delete(regions1_this, inds_nan, axis=0)
                models1_this = np.delete(models1_this, inds_nan, axis=0)

        data1_vals_this = np.ravel(data1_vals_this)
        regions1_this = np.ravel(regions1_this)
        models1_this = np.ravel(models1_this)

        data1_vals = np.append(data1_vals, data1_vals_this)
        regions1 = np.append(regions1, regions1_this)
        models1 = np.append(models1, models1_this)

    data1 = pd.DataFrame({
                'Value': data1_vals,
                'Region': regions1,
                'Model': models1})

    data2_vals = np.full([3*12], np.nan)
    regions2 = np.full([3*12], '', dtype=object)
    models2 = np.full([3*12], '', dtype=object)
    for xx in range(12):
        data2_vals[xx*3+0] = cams_ra_raw_pm25_con_all_overall_rmse[xx]
        data2_vals[xx*3+1] = cams_fc_raw_pm25_con_all_overall_rmse[xx]
        data2_vals[xx*3+2] = cmaq_fc_raw_pm25_con_all_overall_rmse[xx]
        regions2[xx*3+0:xx*3+3] = regions[xx]
        models2[xx*3+0] = lglab_raw[0]
        models2[xx*3+1] = lglab_raw[1]
        models2[xx*3+2] = lglab_raw[2]
    data2 = pd.DataFrame({
                'Value': data2_vals,
                'Region': regions2,
                'Model': models2})

    data3_vals = np.full([3*12], np.nan)
    regions3 = regions2
    models3 = models2
    for xx in range(12):
        data3_vals[xx*3+0] = cams_ra_raw_pm25_con_all_overall_mbe[xx]
        data3_vals[xx*3+1] = cams_fc_raw_pm25_con_all_overall_mbe[xx]
        data3_vals[xx*3+2] = cmaq_fc_raw_pm25_con_all_overall_mbe[xx]
    data3 = pd.DataFrame({
                'Value': data3_vals,
                'Region': regions3,
                'Model': models3})
    plot_stat_box_bar_3panel_vert(fname, data1, data2, data3, regions, xlab1_raw, xlab2_raw, xlab3_raw,
        ylab_pm25_conc, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1_box, title2, title3, colors_raw,
        ymin_box=ymin_pm25_box, yscale_box=yscale_box_pm25)

    ## Now make the plots against validation stations and include the bias-corrected forecasts,
    ## one region at a time per 3-panel plot
    fname = plot_dir.joinpath('overall_o3_allreg_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle = 'Overall Stats By Region for '+mpl_o3+' at Validation Stations'

    data1_vals = np.full([n_valid_times*n_obs_o3_val, 6], np.nan)
    data1_vals[:,0] = np.ravel(cams_ra_raw_o3_con_val[:,:])
    data1_vals[:,1] = np.ravel(cams_fc_raw_o3_con_val[:,:])
    data1_vals[:,2] = np.ravel(cmaq_fc_raw_o3_con_val[:,:])
    data1_vals[:,3] = np.ravel(cams_fc_bm3_o3_con_val[:,:])
    data1_vals[:,4] = np.ravel(cmaq_fc_bm3_o3_con_val[:,:])
    data1_vals[:,5] = np.ravel(obs_o3_con_val[:,:])
    regions1 = np.full([n_valid_times*n_obs_o3_val, 6], regions[0], dtype=object)
    models1 = np.full([n_valid_times*n_obs_o3_val, 6], '', dtype=object)
    models1[:,0] = lglab_raw_bm3_obs[0]
    models1[:,1] = lglab_raw_bm3_obs[1]
    models1[:,2] = lglab_raw_bm3_obs[2]
    models1[:,3] = lglab_raw_bm3_obs[3]
    models1[:,4] = lglab_raw_bm3_obs[4]
    models1[:,5] = lglab_raw_bm3_obs[5]
    ## To conserve memory usage, and to ensure that only distributions of obs-model pairs are used,
    ## find rows where any value is NaN and drop that entire row. Obs will have more NaNs than models.
    for mm in range(6):
        if np.isnan(data1_vals[:,mm]).any():
            inds_nan = np.where(np.isnan(data1_vals[:,mm]))[0]
            data1_vals = np.delete(data1_vals, inds_nan, axis=0)
            regions1 = np.delete(regions1, inds_nan, axis=0)
            models1 = np.delete(models1, inds_nan, axis=0)
            
    data1_vals = np.ravel(data1_vals)
    regions1 = np.ravel(regions1)
    models1 = np.ravel(models1)
    ## Loop over the regions to keep building the very long dataframe
    for rr in range(1,12):
        if rr == 1:
            inds = inds_o3_reg01_val
        elif rr == 2:
            inds = inds_o3_reg02_val
        elif rr == 3:
            inds = inds_o3_reg03_val
        elif rr == 4:
            inds = inds_o3_reg04_val
        elif rr == 5:
            inds = inds_o3_reg05_val
        elif rr == 6:
            inds = inds_o3_reg06_val
        elif rr == 7:
            inds = inds_o3_reg07_val
        elif rr == 8:
            inds = inds_o3_reg08_val
        elif rr == 9:
            inds = inds_o3_reg09_val
        elif rr == 10:
            inds = inds_o3_reg10_val
        elif rr == 11:
            inds = inds_o3_reg99_val
        n_obs = len(inds)
        data1_vals_this = np.full([n_valid_times*n_obs, 6], np.nan)
        data1_vals_this[:,0] = np.ravel(cams_ra_raw_o3_con_val[:,inds])
        data1_vals_this[:,1] = np.ravel(cams_fc_raw_o3_con_val[:,inds])
        data1_vals_this[:,2] = np.ravel(cmaq_fc_raw_o3_con_val[:,inds])
        data1_vals_this[:,3] = np.ravel(cams_fc_bm3_o3_con_val[:,inds])
        data1_vals_this[:,4] = np.ravel(cmaq_fc_bm3_o3_con_val[:,inds])
        data1_vals_this[:,5] = np.ravel(obs_o3_con_val[:,inds])
        regions1_this = np.full([n_valid_times*n_obs, 6], regions[rr], dtype=object)
        models1_this = np.full([n_valid_times*n_obs, 6], '', dtype=object)
        models1_this[:,0] = lglab_raw_bm3_obs[0]
        models1_this[:,1] = lglab_raw_bm3_obs[1]
        models1_this[:,2] = lglab_raw_bm3_obs[2]
        models1_this[:,3] = lglab_raw_bm3_obs[3]
        models1_this[:,4] = lglab_raw_bm3_obs[4]
        models1_this[:,5] = lglab_raw_bm3_obs[5]

        ## To conserve memory usage, and to ensure that only distributions of valid obs-model pairs are used,
        ## find rows where any value is NaN and drop that entire row. Obs will have more NaNs than models.
        for mm in range(6):
            if np.isnan(data1_vals_this[:,mm]).any():
                inds_nan = np.where(np.isnan(data1_vals_this[:,mm]))[0]
                data1_vals_this = np.delete(data1_vals_this, inds_nan, axis=0)
                regions1_this = np.delete(regions1_this, inds_nan, axis=0)
                models1_this = np.delete(models1_this, inds_nan, axis=0)

        data1_vals_this = np.ravel(data1_vals_this)
        regions1_this = np.ravel(regions1_this)
        models1_this = np.ravel(models1_this)

        data1_vals = np.append(data1_vals, data1_vals_this)
        regions1 = np.append(regions1, regions1_this)
        models1 = np.append(models1, models1_this)

    data1 = pd.DataFrame({
                'Region': regions1,
                'Value': data1_vals,
                'Model': models1})

    data2_vals = np.full([5*12], np.nan)
    regions2 = np.full([5*12], '', dtype=object)
    models2 = np.full([5*12], '', dtype=object)
    for xx in range(12):
        data2_vals[xx*5+0] = cams_ra_raw_o3_con_val_overall_rmse[xx]
        data2_vals[xx*5+1] = cams_fc_raw_o3_con_val_overall_rmse[xx]
        data2_vals[xx*5+2] = cmaq_fc_raw_o3_con_val_overall_rmse[xx]
        data2_vals[xx*5+3] = cams_fc_bm3_o3_con_val_overall_rmse[xx]
        data2_vals[xx*5+4] = cmaq_fc_bm3_o3_con_val_overall_rmse[xx]
        regions2[xx*5+0:xx*5+5] = regions[xx]
        models2[xx*5+0] = lglab_raw_bm3[0]
        models2[xx*5+1] = lglab_raw_bm3[1]
        models2[xx*5+2] = lglab_raw_bm3[2]
        models2[xx*5+3] = lglab_raw_bm3[3]
        models2[xx*5+4] = lglab_raw_bm3[4]
    data2 = pd.DataFrame({
                'Value': data2_vals,
                'Region': regions2,
                'Model': models2})

    data3_vals = np.full([5*12], np.nan)
    regions3 = regions2
    models3 = models2
    for xx in range(12):
        data3_vals[xx*5+0] = cams_ra_raw_o3_con_val_overall_mbe[xx]
        data3_vals[xx*5+1] = cams_fc_raw_o3_con_val_overall_mbe[xx]
        data3_vals[xx*5+2] = cmaq_fc_raw_o3_con_val_overall_mbe[xx]
        data3_vals[xx*5+3] = cams_fc_bm3_o3_con_val_overall_mbe[xx]
        data3_vals[xx*5+4] = cmaq_fc_bm3_o3_con_val_overall_mbe[xx]
    data3 = pd.DataFrame({
                'Value': data3_vals,
                'Region': regions3,
                'Model': models3})
    plot_stat_box_bar_3panel_vert(fname, data1, data2, data3, regions, xlab1_raw_bm3, xlab2_raw_bm3, xlab3_raw_bm3,
        ylab_o3_conc, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1_box, title2, title3, colors_raw_bm3,
        ymin_box=ymin_o3_box, yscale_box=yscale_box_o3)


    fname = plot_dir.joinpath('overall_pm25_allreg_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle = 'Overall Stats By Region for '+mpl_pm25+' at Validation Stations'

    data1_vals = np.full([n_valid_times*n_obs_pm25_val, 6], np.nan)
    data1_vals[:,0] = np.ravel(cams_ra_raw_pm25_con_val[:,:])
    data1_vals[:,1] = np.ravel(cams_fc_raw_pm25_con_val[:,:])
    data1_vals[:,2] = np.ravel(cmaq_fc_raw_pm25_con_val[:,:])
    data1_vals[:,3] = np.ravel(cams_fc_bm3_pm25_con_val[:,:])
    data1_vals[:,4] = np.ravel(cmaq_fc_bm3_pm25_con_val[:,:])
    data1_vals[:,5] = np.ravel(obs_pm25_con_val[:,:])
    regions1 = np.full([n_valid_times*n_obs_pm25_val, 6], regions[0], dtype=object)
    models1 = np.full([n_valid_times*n_obs_pm25_val, 6], '', dtype=object)
    models1[:,0] = lglab_raw_bm3_obs[0]
    models1[:,1] = lglab_raw_bm3_obs[1]
    models1[:,2] = lglab_raw_bm3_obs[2]
    models1[:,3] = lglab_raw_bm3_obs[3]
    models1[:,4] = lglab_raw_bm3_obs[4]
    models1[:,5] = lglab_raw_bm3_obs[5]
    ## To conserve memory usage, and to ensure that only distributions of obs-model pairs are used,
    ## find rows where any value is NaN and drop that entire row. Obs will have more NaNs than models.
    for mm in range(6):
        if np.isnan(data1_vals[:,mm]).any():
            inds_nan = np.where(np.isnan(data1_vals[:,mm]))[0]
            data1_vals = np.delete(data1_vals, inds_nan, axis=0)
            regions1 = np.delete(regions1, inds_nan, axis=0)
            models1 = np.delete(models1, inds_nan, axis=0)
            
    data1_vals = np.ravel(data1_vals)
    regions1 = np.ravel(regions1)
    models1 = np.ravel(models1)
    ## Loop over the regions to keep building the very long dataframe
    for rr in range(1,12):
        if rr == 1:
            inds = inds_pm25_reg01_val
        elif rr == 2:
            inds = inds_pm25_reg02_val
        elif rr == 3:
            inds = inds_pm25_reg03_val
        elif rr == 4:
            inds = inds_pm25_reg04_val
        elif rr == 5:
            inds = inds_pm25_reg05_val
        elif rr == 6:
            inds = inds_pm25_reg06_val
        elif rr == 7:
            inds = inds_pm25_reg07_val
        elif rr == 8:
            inds = inds_pm25_reg08_val
        elif rr == 9:
            inds = inds_pm25_reg09_val
        elif rr == 10:
            inds = inds_pm25_reg10_val
        elif rr == 11:
            inds = inds_pm25_reg99_val
        n_obs = len(inds)
        data1_vals_this = np.full([n_valid_times*n_obs, 6], np.nan)
        data1_vals_this[:,0] = np.ravel(cams_ra_raw_pm25_con_val[:,inds])
        data1_vals_this[:,1] = np.ravel(cams_fc_raw_pm25_con_val[:,inds])
        data1_vals_this[:,2] = np.ravel(cmaq_fc_raw_pm25_con_val[:,inds])
        data1_vals_this[:,3] = np.ravel(cams_fc_bm3_pm25_con_val[:,inds])
        data1_vals_this[:,4] = np.ravel(cmaq_fc_bm3_pm25_con_val[:,inds])
        data1_vals_this[:,5] = np.ravel(obs_pm25_con_val[:,inds])
        regions1_this = np.full([n_valid_times*n_obs, 6], regions[rr], dtype=object)
        models1_this = np.full([n_valid_times*n_obs, 6], '', dtype=object)
        models1_this[:,0] = lglab_raw_bm3_obs[0]
        models1_this[:,1] = lglab_raw_bm3_obs[1]
        models1_this[:,2] = lglab_raw_bm3_obs[2]
        models1_this[:,3] = lglab_raw_bm3_obs[3]
        models1_this[:,4] = lglab_raw_bm3_obs[4]
        models1_this[:,5] = lglab_raw_bm3_obs[5]

        ## To conserve memory usage, and to ensure that only distributions of valid obs-model pairs are used,
        ## find rows where any value is NaN and drop that entire row. Obs will have more NaNs than models.
        for mm in range(6):
            if np.isnan(data1_vals_this[:,mm]).any():
                inds_nan = np.where(np.isnan(data1_vals_this[:,mm]))[0]
                data1_vals_this = np.delete(data1_vals_this, inds_nan, axis=0)
                regions1_this = np.delete(regions1_this, inds_nan, axis=0)
                models1_this = np.delete(models1_this, inds_nan, axis=0)

        data1_vals_this = np.ravel(data1_vals_this)
        regions1_this = np.ravel(regions1_this)
        models1_this = np.ravel(models1_this)

        data1_vals = np.append(data1_vals, data1_vals_this)
        regions1 = np.append(regions1, regions1_this)
        models1 = np.append(models1, models1_this)

    data1 = pd.DataFrame({
                'Value': data1_vals,
                'Region': regions1,
                'Model': models1})

    data2_vals = np.full([5*12], np.nan)
    regions2 = np.full([5*12], '', dtype=object)
    models2 = np.full([5*12], '', dtype=object)
    for xx in range(12):
        data2_vals[xx*5+0] = cams_ra_raw_pm25_con_val_overall_rmse[xx]
        data2_vals[xx*5+1] = cams_fc_raw_pm25_con_val_overall_rmse[xx]
        data2_vals[xx*5+2] = cmaq_fc_raw_pm25_con_val_overall_rmse[xx]
        data2_vals[xx*5+3] = cams_fc_bm3_pm25_con_val_overall_rmse[xx]
        data2_vals[xx*5+4] = cmaq_fc_bm3_pm25_con_val_overall_rmse[xx]
        regions2[xx*5+0:xx*5+5] = regions[xx]
        models2[xx*5+0] = lglab_raw_bm3[0]
        models2[xx*5+1] = lglab_raw_bm3[1]
        models2[xx*5+2] = lglab_raw_bm3[2]
        models2[xx*5+3] = lglab_raw_bm3[3]
        models2[xx*5+4] = lglab_raw_bm3[4]
    data2 = pd.DataFrame({
                'Value': data2_vals,
                'Region': regions2,
                'Model': models2})

    data3_vals = np.full([5*12], np.nan)
    regions3 = regions2
    models3 = models2
    for xx in range(12):
        data3_vals[xx*5+0] = cams_ra_raw_pm25_con_val_overall_mbe[xx]
        data3_vals[xx*5+1] = cams_fc_raw_pm25_con_val_overall_mbe[xx]
        data3_vals[xx*5+2] = cmaq_fc_raw_pm25_con_val_overall_mbe[xx]
        data3_vals[xx*5+3] = cams_fc_bm3_pm25_con_val_overall_mbe[xx]
        data3_vals[xx*5+4] = cmaq_fc_bm3_pm25_con_val_overall_mbe[xx]
    data3 = pd.DataFrame({
                'Value': data3_vals,
                'Region': regions3,
                'Model': models3})
    plot_stat_box_bar_3panel_vert(fname, data1, data2, data3, regions, xlab1_raw_bm3, xlab2_raw_bm3, xlab3_raw_bm3,
        ylab_pm25_conc, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1_box, title2, title3, colors_raw_bm3,
        ymin_box=ymin_pm25_box, yscale_box=yscale_box_pm25)


    ## Start by making the plots using all stations and excluding the bias-corrected forecasts
    fname = plot_dir.joinpath('overall_o3_conus_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Overall Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'All Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_times*n_obs_o3_all, 4], np.nan)
    data1[:,0] = np.ravel(cams_ra_raw_o3_con_all[:,:])
    data1[:,1] = np.ravel(cams_fc_raw_o3_con_all[:,:])
    data1[:,2] = np.ravel(cmaq_fc_raw_o3_con_all[:,:])
    data1[:,3] = np.ravel(obs_o3_con_all[:,:])
    data2 = np.full([3], np.nan)
    data2[0] = cams_ra_raw_o3_con_all_overall_rmse[0]
    data2[1] = cams_fc_raw_o3_con_all_overall_rmse[0]
    data2[2] = cmaq_fc_raw_o3_con_all_overall_rmse[0]
    data3 = np.full([3], np.nan)
    data3[0] = cams_ra_raw_o3_con_all_overall_mbe[0]
    data3[1] = cams_fc_raw_o3_con_all_overall_mbe[0]
    data3[2] = cmaq_fc_raw_o3_con_all_overall_mbe[0]
    plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw, xlab2_raw, xlab3_raw,
        ylab_o3_conc, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1_box, title2, title3, colors_raw)

    fname = plot_dir.joinpath('overall_pm25_conus_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Overall Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'All Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_times*n_obs_pm25_all, 4], np.nan)
    data1[:,0] = np.ravel(cams_ra_raw_pm25_con_all[:,:])
    data1[:,1] = np.ravel(cams_fc_raw_pm25_con_all[:,:])
    data1[:,2] = np.ravel(cmaq_fc_raw_pm25_con_all[:,:])
    data1[:,3] = np.ravel(obs_pm25_con_all[:,:])
    data2 = np.full([3], np.nan)
    data2[0] = cams_ra_raw_pm25_con_all_overall_rmse[0]
    data2[1] = cams_fc_raw_pm25_con_all_overall_rmse[0]
    data2[2] = cmaq_fc_raw_pm25_con_all_overall_rmse[0]
    data3 = np.full([3], np.nan)
    data3[0] = cams_ra_raw_pm25_con_all_overall_mbe[0]
    data3[1] = cams_fc_raw_pm25_con_all_overall_mbe[0]
    data3[2] = cmaq_fc_raw_pm25_con_all_overall_mbe[0]
    plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw, xlab2_raw, xlab3_raw,
        ylab_pm25_conc, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1_box, title2, title3, colors_raw,
        ymin_box=ymin_pm25_box, yscale_box=yscale_box_pm25)

    ## Regions 1-10
    for rr in range(1,11):
        reg_str = 'reg'+str(rr).zfill(2)
        if rr == 1:
            inds_o3_reg = inds_o3_reg01_all
            inds_pm25_reg = inds_pm25_reg01_all
        elif rr == 2:
            inds_o3_reg = inds_o3_reg02_all
            inds_pm25_reg = inds_pm25_reg02_all
        elif rr == 3:
            inds_o3_reg = inds_o3_reg03_all
            inds_pm25_reg = inds_pm25_reg03_all
        elif rr == 4:
            inds_o3_reg = inds_o3_reg04_all
            inds_pm25_reg = inds_pm25_reg04_all
        elif rr == 5:
            inds_o3_reg = inds_o3_reg05_all
            inds_pm25_reg = inds_pm25_reg05_all
        elif rr == 6:
            inds_o3_reg = inds_o3_reg06_all
            inds_pm25_reg = inds_pm25_reg06_all
        elif rr == 7:
            inds_o3_reg = inds_o3_reg07_all
            inds_pm25_reg = inds_pm25_reg07_all
        elif rr == 8:
            inds_o3_reg = inds_o3_reg08_all
            inds_pm25_reg = inds_pm25_reg08_all
        elif rr == 9:
            inds_o3_reg = inds_o3_reg09_all
            inds_pm25_reg = inds_pm25_reg09_all
        elif rr == 10:
            inds_o3_reg = inds_o3_reg10_all
            inds_pm25_reg = inds_pm25_reg10_all
        n_obs_o3_reg = len(inds_o3_reg)
        n_obs_pm25_reg = len(inds_pm25_reg)

        fname = plot_dir.joinpath('overall_o3_'+reg_str+'_3panel_'+date_range_file+'_all.'+plot_type)
        suptitle_prefix = 'Overall Stats for '+mpl_o3+' '+em_dash+' '
        suptitle_suffix = 'All Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([n_valid_times*n_obs_o3_reg, 4], np.nan)
        data1[:,0] = np.ravel(cams_ra_raw_o3_con_all[:,inds_o3_reg])
        data1[:,1] = np.ravel(cams_fc_raw_o3_con_all[:,inds_o3_reg])
        data1[:,2] = np.ravel(cmaq_fc_raw_o3_con_all[:,inds_o3_reg])
        data1[:,3] = np.ravel(obs_o3_con_all[:,inds_o3_reg])
        data2 = np.full([3], np.nan)
        data2[0] = cams_ra_raw_o3_con_all_overall_rmse[rr]
        data2[1] = cams_fc_raw_o3_con_all_overall_rmse[rr]
        data2[2] = cmaq_fc_raw_o3_con_all_overall_rmse[rr]
        data3 = np.full([3], np.nan)
        data3[0] = cams_ra_raw_o3_con_all_overall_mbe[rr]
        data3[1] = cams_fc_raw_o3_con_all_overall_mbe[rr]
        data3[2] = cmaq_fc_raw_o3_con_all_overall_mbe[rr]
        plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw, xlab2_raw, xlab3_raw,
            ylab_o3_conc, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1_box, title2, title3, colors_raw)

        fname = plot_dir.joinpath('overall_pm25_'+reg_str+'_3panel_'+date_range_file+'_all.'+plot_type)
        suptitle_prefix = 'Overall Stats for '+mpl_pm25+' '+em_dash+' '
        suptitle_suffix = 'All Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([n_valid_times*n_obs_pm25_reg, 4], np.nan)
        data1[:,0] = np.ravel(cams_ra_raw_pm25_con_all[:,inds_pm25_reg])
        data1[:,1] = np.ravel(cams_fc_raw_pm25_con_all[:,inds_pm25_reg])
        data1[:,2] = np.ravel(cmaq_fc_raw_pm25_con_all[:,inds_pm25_reg])
        data1[:,3] = np.ravel(obs_pm25_con_all[:,inds_pm25_reg])
        data2 = np.full([3], np.nan)
        data2[0] = cams_ra_raw_pm25_con_all_overall_rmse[rr]
        data2[1] = cams_fc_raw_pm25_con_all_overall_rmse[rr]
        data2[2] = cmaq_fc_raw_pm25_con_all_overall_rmse[rr]
        data3 = np.full([3], np.nan)
        data3[0] = cams_ra_raw_pm25_con_all_overall_mbe[rr]
        data3[1] = cams_fc_raw_pm25_con_all_overall_mbe[rr]
        data3[2] = cmaq_fc_raw_pm25_con_all_overall_mbe[rr]
        plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw, xlab2_raw, xlab3_raw,
            ylab_pm25_conc, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1_box, title2, title3, colors_raw,
            ymin_box=ymin_pm25_box, yscale_box=yscale_box_pm25)

    ## Region 99 (outside U.S.)
    n_obs_o3_reg99_all = len(inds_o3_reg99_all)
    n_obs_pm25_reg99_all = len(inds_pm25_reg99_all)
    fname = plot_dir.joinpath('overall_o3_canmex_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Overall Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'All Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_times*n_obs_o3_reg99_all, 4], np.nan)
    data1[:,0] = np.ravel(cams_ra_raw_o3_con_all[:,inds_o3_reg99_all])
    data1[:,1] = np.ravel(cams_fc_raw_o3_con_all[:,inds_o3_reg99_all])
    data1[:,2] = np.ravel(cmaq_fc_raw_o3_con_all[:,inds_o3_reg99_all])
    data1[:,3] = np.ravel(obs_o3_con_all[:,inds_o3_reg99_all])
    data2 = np.full([3], np.nan)
    data2[0] = cams_ra_raw_o3_con_all_overall_rmse[11]
    data2[1] = cams_fc_raw_o3_con_all_overall_rmse[11]
    data2[2] = cmaq_fc_raw_o3_con_all_overall_rmse[11]
    data3 = np.full([3], np.nan)
    data3[0] = cams_ra_raw_o3_con_all_overall_mbe[11]
    data3[1] = cams_fc_raw_o3_con_all_overall_mbe[11]
    data3[2] = cmaq_fc_raw_o3_con_all_overall_mbe[11]
    plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw, xlab2_raw, xlab3_raw,
        ylab_o3_conc, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1_box, title2, title3, colors_raw)

    fname = plot_dir.joinpath('overall_pm25_canmex_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Overall Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'All Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_times*n_obs_pm25_reg99_all, 4], np.nan)
    data1[:,0] = np.ravel(cams_ra_raw_pm25_con_all[:,inds_pm25_reg99_all])
    data1[:,1] = np.ravel(cams_fc_raw_pm25_con_all[:,inds_pm25_reg99_all])
    data1[:,2] = np.ravel(cmaq_fc_raw_pm25_con_all[:,inds_pm25_reg99_all])
    data1[:,3] = np.ravel(obs_pm25_con_all[:,inds_pm25_reg99_all])
    data2 = np.full([3], np.nan)
    data2[0] = cams_ra_raw_pm25_con_all_overall_rmse[11]
    data2[1] = cams_fc_raw_pm25_con_all_overall_rmse[11]
    data2[2] = cmaq_fc_raw_pm25_con_all_overall_rmse[11]
    data3 = np.full([3], np.nan)
    data3[0] = cams_ra_raw_pm25_con_all_overall_mbe[11]
    data3[1] = cams_fc_raw_pm25_con_all_overall_mbe[11]
    data3[2] = cmaq_fc_raw_pm25_con_all_overall_mbe[11]
    plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw, xlab2_raw, xlab3_raw,
        ylab_pm25_conc, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1_box, title2, title3, colors_raw,
        ymin_box=ymin_pm25_box, yscale_box=yscale_box_pm25)

    ## Now make the plots against validation stations and include the bias-corrected forecasts
    fname = plot_dir.joinpath('overall_o3_conus_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Overall Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_times*n_obs_o3_val, 6], np.nan)
    data1[:,0] = np.ravel(cams_ra_raw_o3_con_val[:,:])
    data1[:,1] = np.ravel(cams_fc_raw_o3_con_val[:,:])
    data1[:,2] = np.ravel(cmaq_fc_raw_o3_con_val[:,:])
    data1[:,3] = np.ravel(cams_fc_bm3_o3_con_val[:,:])
    data1[:,4] = np.ravel(cmaq_fc_bm3_o3_con_val[:,:])
    data1[:,5] = np.ravel(obs_o3_con_val[:,:])
    data2 = np.full([5], np.nan)
    data2[0] = cams_ra_raw_o3_con_val_overall_rmse[0]
    data2[1] = cams_fc_raw_o3_con_val_overall_rmse[0]
    data2[2] = cmaq_fc_raw_o3_con_val_overall_rmse[0]
    data2[3] = cams_fc_bm3_o3_con_val_overall_rmse[0]
    data2[4] = cmaq_fc_bm3_o3_con_val_overall_rmse[0]
    data3 = np.full([5], np.nan)
    data3[0] = cams_ra_raw_o3_con_val_overall_mbe[0]
    data3[1] = cams_fc_raw_o3_con_val_overall_mbe[0]
    data3[2] = cmaq_fc_raw_o3_con_val_overall_mbe[0]
    data3[3] = cams_fc_bm3_o3_con_val_overall_mbe[0]
    data3[4] = cmaq_fc_bm3_o3_con_val_overall_mbe[0]
    plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw_bm3, xlab2_raw_bm3, xlab3_raw_bm3,
        ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1_box, title2, title3, colors_raw_bm3)

    fname = plot_dir.joinpath('overall_pm25_conus_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Overall Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_times*n_obs_pm25_val, 6], np.nan)
    data1[:,0] = np.ravel(cams_ra_raw_pm25_con_val[:,:])
    data1[:,1] = np.ravel(cams_fc_raw_pm25_con_val[:,:])
    data1[:,2] = np.ravel(cmaq_fc_raw_pm25_con_val[:,:])
    data1[:,3] = np.ravel(cams_fc_bm3_pm25_con_val[:,:])
    data1[:,4] = np.ravel(cmaq_fc_bm3_pm25_con_val[:,:])
    data1[:,5] = np.ravel(obs_pm25_con_val[:,:])
    data2 = np.full([5], np.nan)
    data2[0] = cams_ra_raw_pm25_con_val_overall_rmse[0]
    data2[1] = cams_fc_raw_pm25_con_val_overall_rmse[0]
    data2[2] = cmaq_fc_raw_pm25_con_val_overall_rmse[0]
    data2[3] = cams_fc_bm3_pm25_con_val_overall_rmse[0]
    data2[4] = cmaq_fc_bm3_pm25_con_val_overall_rmse[0]
    data3 = np.full([5], np.nan)
    data3[0] = cams_ra_raw_pm25_con_val_overall_mbe[0]
    data3[1] = cams_fc_raw_pm25_con_val_overall_mbe[0]
    data3[2] = cmaq_fc_raw_pm25_con_val_overall_mbe[0]
    data3[3] = cams_fc_bm3_pm25_con_val_overall_mbe[0]
    data3[4] = cmaq_fc_bm3_pm25_con_val_overall_mbe[0]
    plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw_bm3, xlab2_raw_bm3, xlab3_raw_bm3,
        ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1_box, title2, title3, colors_raw_bm3,
        ymin_box=ymin_pm25_box, yscale_box=yscale_box_pm25)

    ## Regions 1-10
    for rr in range(1,11):
        reg_str = 'reg'+str(rr).zfill(2)
        if rr == 1:
            inds_o3_reg = inds_o3_reg01_val
            inds_pm25_reg = inds_pm25_reg01_val
        elif rr == 2:
            inds_o3_reg = inds_o3_reg02_val
            inds_pm25_reg = inds_pm25_reg02_val
        elif rr == 3:
            inds_o3_reg = inds_o3_reg03_val
            inds_pm25_reg = inds_pm25_reg03_val
        elif rr == 4:
            inds_o3_reg = inds_o3_reg04_val
            inds_pm25_reg = inds_pm25_reg04_val
        elif rr == 5:
            inds_o3_reg = inds_o3_reg05_val
            inds_pm25_reg = inds_pm25_reg05_val
        elif rr == 6:
            inds_o3_reg = inds_o3_reg06_val
            inds_pm25_reg = inds_pm25_reg06_val
        elif rr == 7:
            inds_o3_reg = inds_o3_reg07_val
            inds_pm25_reg = inds_pm25_reg07_val
        elif rr == 8:
            inds_o3_reg = inds_o3_reg08_val
            inds_pm25_reg = inds_pm25_reg08_val
        elif rr == 9:
            inds_o3_reg = inds_o3_reg09_val
            inds_pm25_reg = inds_pm25_reg09_val
        elif rr == 10:
            inds_o3_reg = inds_o3_reg10_val
            inds_pm25_reg = inds_pm25_reg10_val
        n_obs_o3_reg = len(inds_o3_reg)
        n_obs_pm25_reg = len(inds_pm25_reg)

        fname = plot_dir.joinpath('overall_o3_'+reg_str+'_3panel_'+date_range_file+'_val.'+plot_type)
        suptitle_prefix = 'Overall Stats for '+mpl_o3+' '+em_dash+' '
        suptitle_suffix = 'Validation Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([n_valid_times*n_obs_o3_reg, 6], np.nan)
        data1[:,0] = np.ravel(cams_ra_raw_o3_con_val[:,inds_o3_reg])
        data1[:,1] = np.ravel(cams_fc_raw_o3_con_val[:,inds_o3_reg])
        data1[:,2] = np.ravel(cmaq_fc_raw_o3_con_val[:,inds_o3_reg])
        data1[:,3] = np.ravel(cams_fc_bm3_o3_con_val[:,inds_o3_reg])
        data1[:,4] = np.ravel(cmaq_fc_bm3_o3_con_val[:,inds_o3_reg])
        data1[:,5] = np.ravel(obs_o3_con_val[:,inds_o3_reg])
        data2 = np.full([5], np.nan)
        data2[0] = cams_ra_raw_o3_con_val_overall_rmse[rr]
        data2[1] = cams_fc_raw_o3_con_val_overall_rmse[rr]
        data2[2] = cmaq_fc_raw_o3_con_val_overall_rmse[rr]
        data2[3] = cams_fc_bm3_o3_con_val_overall_rmse[rr]
        data2[4] = cmaq_fc_bm3_o3_con_val_overall_rmse[rr]
        data3 = np.full([5], np.nan)
        data3[0] = cams_ra_raw_o3_con_val_overall_mbe[rr]
        data3[1] = cams_fc_raw_o3_con_val_overall_mbe[rr]
        data3[2] = cmaq_fc_raw_o3_con_val_overall_mbe[rr]
        data3[3] = cams_fc_bm3_o3_con_val_overall_mbe[rr]
        data3[4] = cmaq_fc_bm3_o3_con_val_overall_mbe[rr]
        plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw_bm3, xlab2_raw_bm3, xlab3_raw_bm3,
            ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1_box, title2, title3, colors_raw_bm3)

        fname = plot_dir.joinpath('overall_pm25_'+reg_str+'_3panel_'+date_range_file+'_val.'+plot_type)
        suptitle_prefix = 'Overall Stats for '+mpl_pm25+' '+em_dash+' '
        suptitle_suffix = 'Validation Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([n_valid_times*n_obs_pm25_reg, 6], np.nan)
        data1[:,0] = np.ravel(cams_ra_raw_pm25_con_val[:,inds_pm25_reg])
        data1[:,1] = np.ravel(cams_fc_raw_pm25_con_val[:,inds_pm25_reg])
        data1[:,2] = np.ravel(cmaq_fc_raw_pm25_con_val[:,inds_pm25_reg])
        data1[:,3] = np.ravel(cams_fc_bm3_pm25_con_val[:,inds_pm25_reg])
        data1[:,4] = np.ravel(cmaq_fc_bm3_pm25_con_val[:,inds_pm25_reg])
        data1[:,5] = np.ravel(obs_pm25_con_val[:,inds_pm25_reg])
        data2 = np.full([5], np.nan)
        data2[0] = cams_ra_raw_pm25_con_val_overall_rmse[rr]
        data2[1] = cams_fc_raw_pm25_con_val_overall_rmse[rr]
        data2[2] = cmaq_fc_raw_pm25_con_val_overall_rmse[rr]
        data2[3] = cams_fc_bm3_pm25_con_val_overall_rmse[rr]
        data2[4] = cmaq_fc_bm3_pm25_con_val_overall_rmse[rr]
        data3 = np.full([5], np.nan)
        data3[0] = cams_ra_raw_pm25_con_val_overall_mbe[rr]
        data3[1] = cams_fc_raw_pm25_con_val_overall_mbe[rr]
        data3[2] = cmaq_fc_raw_pm25_con_val_overall_mbe[rr]
        data3[3] = cams_fc_bm3_pm25_con_val_overall_mbe[rr]
        data3[4] = cmaq_fc_bm3_pm25_con_val_overall_mbe[rr]
        plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw_bm3, xlab2_raw_bm3, xlab3_raw_bm3,
            ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1_box, title2, title3, colors_raw_bm3,
            ymin_box=ymin_pm25_box, yscale_box=yscale_box_pm25)

    ## Region 99 (outside U.S.)
    n_obs_o3_reg99_val = len(inds_o3_reg99_val)
    n_obs_pm25_reg99_val = len(inds_pm25_reg99_val)
    fname = plot_dir.joinpath('overall_o3_canmex_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Overall Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_times*n_obs_o3_reg99_val, 6], np.nan)
    data1[:,0] = np.ravel(cams_ra_raw_o3_con_val[:,inds_o3_reg99_val])
    data1[:,1] = np.ravel(cams_fc_raw_o3_con_val[:,inds_o3_reg99_val])
    data1[:,2] = np.ravel(cmaq_fc_raw_o3_con_val[:,inds_o3_reg99_val])
    data1[:,3] = np.ravel(cams_fc_bm3_o3_con_val[:,inds_o3_reg99_val])
    data1[:,4] = np.ravel(cmaq_fc_bm3_o3_con_val[:,inds_o3_reg99_val])
    data1[:,5] = np.ravel(obs_o3_con_val[:,inds_o3_reg99_val])
    data2 = np.full([5], np.nan)
    data2[0] = cams_ra_raw_o3_con_val_overall_rmse[11]
    data2[1] = cams_fc_raw_o3_con_val_overall_rmse[11]
    data2[2] = cmaq_fc_raw_o3_con_val_overall_rmse[11]
    data2[3] = cams_fc_bm3_o3_con_val_overall_rmse[11]
    data2[4] = cmaq_fc_bm3_o3_con_val_overall_rmse[11]
    data3 = np.full([5], np.nan)
    data3[0] = cams_ra_raw_o3_con_val_overall_mbe[11]
    data3[1] = cams_fc_raw_o3_con_val_overall_mbe[11]
    data3[2] = cmaq_fc_raw_o3_con_val_overall_mbe[11]
    data3[3] = cams_fc_bm3_o3_con_val_overall_mbe[11]
    data3[4] = cmaq_fc_bm3_o3_con_val_overall_mbe[11]
    plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw_bm3, xlab2_raw_bm3, xlab3_raw_bm3,
        ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1_box, title2, title3, colors_raw_bm3)

    fname = plot_dir.joinpath('overall_pm25_canmex_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Overall Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_times*n_obs_pm25_reg99_val, 6], np.nan)
    data1[:,0] = np.ravel(cams_ra_raw_pm25_con_val[:,inds_pm25_reg99_val])
    data1[:,1] = np.ravel(cams_fc_raw_pm25_con_val[:,inds_pm25_reg99_val])
    data1[:,2] = np.ravel(cmaq_fc_raw_pm25_con_val[:,inds_pm25_reg99_val])
    data1[:,3] = np.ravel(cams_fc_bm3_pm25_con_val[:,inds_pm25_reg99_val])
    data1[:,4] = np.ravel(cmaq_fc_bm3_pm25_con_val[:,inds_pm25_reg99_val])
    data1[:,5] = np.ravel(obs_pm25_con_val[:,inds_pm25_reg99_val])
    data2 = np.full([5], np.nan)
    data2[0] = cams_ra_raw_pm25_con_val_overall_rmse[11]
    data2[1] = cams_fc_raw_pm25_con_val_overall_rmse[11]
    data2[2] = cmaq_fc_raw_pm25_con_val_overall_rmse[11]
    data2[3] = cams_fc_bm3_pm25_con_val_overall_rmse[11]
    data2[4] = cmaq_fc_bm3_pm25_con_val_overall_rmse[11]
    data3 = np.full([5], np.nan)
    data3[0] = cams_ra_raw_pm25_con_val_overall_mbe[11]
    data3[1] = cams_fc_raw_pm25_con_val_overall_mbe[11]
    data3[2] = cmaq_fc_raw_pm25_con_val_overall_mbe[11]
    data3[3] = cams_fc_bm3_pm25_con_val_overall_mbe[11]
    data3[4] = cmaq_fc_bm3_pm25_con_val_overall_mbe[11]
    plot_stat_box_bar_3panel_horz(fname, data1, data2, data3, xlab1_raw_bm3, xlab2_raw_bm3, xlab3_raw_bm3,
        ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1_box, title2, title3, colors_raw_bm3,
        ymin_box=ymin_pm25_box, yscale_box=yscale_box_pm25)


    ## Second, make the diurnal plots
    ## Start making the plots for all stations, excluding the bias-corrected forecasts
    ## CONUS
    fname = plot_dir.joinpath('diurnal_o3_conus_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Diurnal-Average Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'All Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([24,4], np.nan)
    data1[:,0] = cams_ra_raw_o3_con_all_diurnal_mean[:,0]
    data1[:,1] = cams_fc_raw_o3_con_all_diurnal_mean[:,0]
    data1[:,2] = cmaq_fc_raw_o3_con_all_diurnal_mean[:,0]
    data1[:,3] = obs_o3_con_all_diurnal_mean[:,0]
    data2 = np.full([24,3], np.nan)
    data2[:,0] = cams_ra_raw_o3_con_all_diurnal_rmse[:,0]
    data2[:,1] = cams_fc_raw_o3_con_all_diurnal_rmse[:,0]
    data2[:,2] = cmaq_fc_raw_o3_con_all_diurnal_rmse[:,0]
    data3 = np.full([24,3], np.nan)
    data3[:,0] = cams_ra_raw_o3_con_all_diurnal_mbe[:,0]
    data3[:,1] = cams_fc_raw_o3_con_all_diurnal_mbe[:,0]
    data3[:,2] = cmaq_fc_raw_o3_con_all_diurnal_mbe[:,0]
    plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
        ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
        xlab_diurnal, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
        colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    fname = plot_dir.joinpath('diurnal_pm25_conus_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Diurnal-Average Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'All Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([24,4], np.nan)
    data1[:,0] = cams_ra_raw_pm25_con_all_diurnal_mean[:,0]
    data1[:,1] = cams_fc_raw_pm25_con_all_diurnal_mean[:,0]
    data1[:,2] = cmaq_fc_raw_pm25_con_all_diurnal_mean[:,0]
    data1[:,3] = obs_pm25_con_all_diurnal_mean[:,0]
    data2 = np.full([24,3], np.nan)
    data2[:,0] = cams_ra_raw_pm25_con_all_diurnal_rmse[:,0]
    data2[:,1] = cams_fc_raw_pm25_con_all_diurnal_rmse[:,0]
    data2[:,2] = cmaq_fc_raw_pm25_con_all_diurnal_rmse[:,0]
    data3 = np.full([24,3], np.nan)
    data3[:,0] = cams_ra_raw_pm25_con_all_diurnal_mbe[:,0]
    data3[:,1] = cams_fc_raw_pm25_con_all_diurnal_mbe[:,0]
    data3[:,2] = cmaq_fc_raw_pm25_con_all_diurnal_mbe[:,0]
    plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
        ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
        xlab_diurnal, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
        colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    ## Regions 1-10
    for rr in range(1,11):
        reg_str = 'reg'+str(rr).zfill(2)
        fname = plot_dir.joinpath('diurnal_o3_'+reg_str+'_3panel_'+date_range_file+'_all.'+plot_type)
        suptitle_prefix = 'Diurnal-Average Stats for '+mpl_o3+' '+em_dash+' '
        suptitle_suffix = 'All Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([24,4], np.nan)
        data1[:,0] = cams_ra_raw_o3_con_all_diurnal_mean[:,rr]
        data1[:,1] = cams_fc_raw_o3_con_all_diurnal_mean[:,rr]
        data1[:,2] = cmaq_fc_raw_o3_con_all_diurnal_mean[:,rr]
        data1[:,3] = obs_o3_con_all_diurnal_mean[:,rr]
        data2 = np.full([24,3], np.nan)
        data2[:,0] = cams_ra_raw_o3_con_all_diurnal_rmse[:,rr]
        data2[:,1] = cams_fc_raw_o3_con_all_diurnal_rmse[:,rr]
        data2[:,2] = cmaq_fc_raw_o3_con_all_diurnal_rmse[:,rr]
        data3 = np.full([24,3], np.nan)
        data3[:,0] = cams_ra_raw_o3_con_all_diurnal_mbe[:,rr]
        data3[:,1] = cams_fc_raw_o3_con_all_diurnal_mbe[:,rr]
        data3[:,2] = cmaq_fc_raw_o3_con_all_diurnal_mbe[:,rr]
        plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
            ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
            xlab_diurnal, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
            colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

        fname = plot_dir.joinpath('diurnal_pm25_'+reg_str+'_3panel_'+date_range_file+'_all.'+plot_type)
        suptitle_prefix = 'Diurnal-Average Stats for '+mpl_pm25+' '+em_dash+' '
        suptitle_suffix = 'All Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([24,4], np.nan)
        data1[:,0] = cams_ra_raw_pm25_con_all_diurnal_mean[:,rr]
        data1[:,1] = cams_fc_raw_pm25_con_all_diurnal_mean[:,rr]
        data1[:,2] = cmaq_fc_raw_pm25_con_all_diurnal_mean[:,rr]
        data1[:,3] = obs_pm25_con_all_diurnal_mean[:,rr]
        data2 = np.full([24,3], np.nan)
        data2[:,0] = cams_ra_raw_pm25_con_all_diurnal_rmse[:,rr]
        data2[:,1] = cams_fc_raw_pm25_con_all_diurnal_rmse[:,rr]
        data2[:,2] = cmaq_fc_raw_pm25_con_all_diurnal_rmse[:,rr]
        data3 = np.full([24,3], np.nan)
        data3[:,0] = cams_ra_raw_pm25_con_all_diurnal_mbe[:,rr]
        data3[:,1] = cams_fc_raw_pm25_con_all_diurnal_mbe[:,rr]
        data3[:,2] = cmaq_fc_raw_pm25_con_all_diurnal_mbe[:,rr]
        plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
            ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
            xlab_diurnal, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
            colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    ## Region 99 (outside U.S.)
    fname = plot_dir.joinpath('diurnal_o3_canmex_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Diurnal-Average Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'All Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([24,4], np.nan)
    data1[:,0] = cams_ra_raw_o3_con_all_diurnal_mean[:,11]
    data1[:,1] = cams_fc_raw_o3_con_all_diurnal_mean[:,11]
    data1[:,2] = cmaq_fc_raw_o3_con_all_diurnal_mean[:,11]
    data1[:,3] = obs_o3_con_all_diurnal_mean[:,11]
    data2 = np.full([24,3], np.nan)
    data2[:,0] = cams_ra_raw_o3_con_all_diurnal_rmse[:,11]
    data2[:,1] = cams_fc_raw_o3_con_all_diurnal_rmse[:,11]
    data2[:,2] = cmaq_fc_raw_o3_con_all_diurnal_rmse[:,11]
    data3 = np.full([24,3], np.nan)
    data3[:,0] = cams_ra_raw_o3_con_all_diurnal_mbe[:,11]
    data3[:,1] = cams_fc_raw_o3_con_all_diurnal_mbe[:,11]
    data3[:,2] = cmaq_fc_raw_o3_con_all_diurnal_mbe[:,11]
    plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
        ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
        xlab_diurnal, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
        colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    fname = plot_dir.joinpath('diurnal_pm25_canmex_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Diurnal-Average Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'All Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([24,4], np.nan)
    data1[:,0] = cams_ra_raw_pm25_con_all_diurnal_mean[:,11]
    data1[:,1] = cams_fc_raw_pm25_con_all_diurnal_mean[:,11]
    data1[:,2] = cmaq_fc_raw_pm25_con_all_diurnal_mean[:,11]
    data1[:,3] = obs_pm25_con_all_diurnal_mean[:,11]
    data2 = np.full([24,3], np.nan)
    data2[:,0] = cams_ra_raw_pm25_con_all_diurnal_rmse[:,11]
    data2[:,1] = cams_fc_raw_pm25_con_all_diurnal_rmse[:,11]
    data2[:,2] = cmaq_fc_raw_pm25_con_all_diurnal_rmse[:,11]
    data3 = np.full([24,3], np.nan)
    data3[:,0] = cams_ra_raw_pm25_con_all_diurnal_mbe[:,11]
    data3[:,1] = cams_fc_raw_pm25_con_all_diurnal_mbe[:,11]
    data3[:,2] = cmaq_fc_raw_pm25_con_all_diurnal_mbe[:,11]
    plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
        ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
        xlab_diurnal, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
        colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    ## Now make the plots for the validation stations, including the bias-corrected forecasts
    ## CONUS
    fname = plot_dir.joinpath('diurnal_o3_conus_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Diurnal-Average Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([24,6], np.nan)
    data1[:,0] = cams_ra_raw_o3_con_val_diurnal_mean[:,0]
    data1[:,1] = cams_fc_raw_o3_con_val_diurnal_mean[:,0]
    data1[:,2] = cmaq_fc_raw_o3_con_val_diurnal_mean[:,0]
    data1[:,3] = cams_fc_bm3_o3_con_val_diurnal_mean[:,0]
    data1[:,4] = cmaq_fc_bm3_o3_con_val_diurnal_mean[:,0]
    data1[:,5] = obs_o3_con_val_diurnal_mean[:,0]
    data2 = np.full([24,5], np.nan)
    data2[:,0] = cams_ra_raw_o3_con_val_diurnal_rmse[:,0]
    data2[:,1] = cams_fc_raw_o3_con_val_diurnal_rmse[:,0]
    data2[:,2] = cmaq_fc_raw_o3_con_val_diurnal_rmse[:,0]
    data2[:,3] = cams_fc_bm3_o3_con_val_diurnal_rmse[:,0]
    data2[:,4] = cmaq_fc_bm3_o3_con_val_diurnal_rmse[:,0]
    data3 = np.full([24,5], np.nan)
    data3[:,0] = cams_ra_raw_o3_con_val_diurnal_mbe[:,0]
    data3[:,1] = cams_fc_raw_o3_con_val_diurnal_mbe[:,0]
    data3[:,2] = cmaq_fc_raw_o3_con_val_diurnal_mbe[:,0]
    data3[:,3] = cams_fc_bm3_o3_con_val_diurnal_mbe[:,0]
    data3[:,4] = cmaq_fc_bm3_o3_con_val_diurnal_mbe[:,0]
    plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
        ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
        xlab_diurnal, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
        colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

    fname = plot_dir.joinpath('diurnal_pm25_conus_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Diurnal-Average Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([24,6], np.nan)
    data1[:,0] = cams_ra_raw_pm25_con_val_diurnal_mean[:,0]
    data1[:,1] = cams_fc_raw_pm25_con_val_diurnal_mean[:,0]
    data1[:,2] = cmaq_fc_raw_pm25_con_val_diurnal_mean[:,0]
    data1[:,3] = cams_fc_bm3_pm25_con_val_diurnal_mean[:,0]
    data1[:,4] = cmaq_fc_bm3_pm25_con_val_diurnal_mean[:,0]
    data1[:,5] = obs_pm25_con_val_diurnal_mean[:,0]
    data2 = np.full([24,5], np.nan)
    data2[:,0] = cams_ra_raw_pm25_con_val_diurnal_rmse[:,0]
    data2[:,1] = cams_fc_raw_pm25_con_val_diurnal_rmse[:,0]
    data2[:,2] = cmaq_fc_raw_pm25_con_val_diurnal_rmse[:,0]
    data2[:,3] = cams_fc_bm3_pm25_con_val_diurnal_rmse[:,0]
    data2[:,4] = cmaq_fc_bm3_pm25_con_val_diurnal_rmse[:,0]
    data3 = np.full([24,5], np.nan)
    data3[:,0] = cams_ra_raw_pm25_con_val_diurnal_mbe[:,0]
    data3[:,1] = cams_fc_raw_pm25_con_val_diurnal_mbe[:,0]
    data3[:,2] = cmaq_fc_raw_pm25_con_val_diurnal_mbe[:,0]
    data3[:,3] = cams_fc_bm3_pm25_con_val_diurnal_mbe[:,0]
    data3[:,4] = cmaq_fc_bm3_pm25_con_val_diurnal_mbe[:,0]
    plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
        ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
        xlab_diurnal, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
        colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

    ## Regions 1-10
    for rr in range(1,11):
        reg_str = 'reg'+str(rr).zfill(2)
        fname = plot_dir.joinpath('diurnal_o3_'+reg_str+'_3panel_'+date_range_file+'_val.'+plot_type)
        suptitle_prefix = 'Diurnal-Average Stats for '+mpl_o3+' '+em_dash+' '
        suptitle_suffix = 'Validation Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([24,6], np.nan)
        data1[:,0] = cams_ra_raw_o3_con_val_diurnal_mean[:,rr]
        data1[:,1] = cams_fc_raw_o3_con_val_diurnal_mean[:,rr]
        data1[:,2] = cmaq_fc_raw_o3_con_val_diurnal_mean[:,rr]
        data1[:,3] = cams_fc_bm3_o3_con_val_diurnal_mean[:,rr]
        data1[:,4] = cmaq_fc_bm3_o3_con_val_diurnal_mean[:,rr]
        data1[:,5] = obs_o3_con_val_diurnal_mean[:,rr]
        data2 = np.full([24,5], np.nan)
        data2[:,0] = cams_ra_raw_o3_con_val_diurnal_rmse[:,rr]
        data2[:,1] = cams_fc_raw_o3_con_val_diurnal_rmse[:,rr]
        data2[:,2] = cmaq_fc_raw_o3_con_val_diurnal_rmse[:,rr]
        data2[:,3] = cams_fc_bm3_o3_con_val_diurnal_rmse[:,rr]
        data2[:,4] = cmaq_fc_bm3_o3_con_val_diurnal_rmse[:,rr]
        data3 = np.full([24,5], np.nan)
        data3[:,0] = cams_ra_raw_o3_con_val_diurnal_mbe[:,rr]
        data3[:,1] = cams_fc_raw_o3_con_val_diurnal_mbe[:,rr]
        data3[:,2] = cmaq_fc_raw_o3_con_val_diurnal_mbe[:,rr]
        data3[:,3] = cams_fc_bm3_o3_con_val_diurnal_mbe[:,rr]
        data3[:,4] = cmaq_fc_bm3_o3_con_val_diurnal_mbe[:,rr]
        plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
            ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
            xlab_diurnal, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
            colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

        fname = plot_dir.joinpath('diurnal_pm25_'+reg_str+'_3panel_'+date_range_file+'_val.'+plot_type)
        suptitle_prefix = 'Diurnal-Average Stats for '+mpl_pm25+' '+em_dash+' '
        suptitle_suffix = 'Validation Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([24,6], np.nan)
        data1[:,0] = cams_ra_raw_pm25_con_val_diurnal_mean[:,rr]
        data1[:,1] = cams_fc_raw_pm25_con_val_diurnal_mean[:,rr]
        data1[:,2] = cmaq_fc_raw_pm25_con_val_diurnal_mean[:,rr]
        data1[:,3] = cams_fc_bm3_pm25_con_val_diurnal_mean[:,rr]
        data1[:,4] = cmaq_fc_bm3_pm25_con_val_diurnal_mean[:,rr]
        data1[:,5] = obs_pm25_con_val_diurnal_mean[:,rr]
        data2 = np.full([24,5], np.nan)
        data2[:,0] = cams_ra_raw_pm25_con_val_diurnal_rmse[:,rr]
        data2[:,1] = cams_fc_raw_pm25_con_val_diurnal_rmse[:,rr]
        data2[:,2] = cmaq_fc_raw_pm25_con_val_diurnal_rmse[:,rr]
        data2[:,3] = cams_fc_bm3_pm25_con_val_diurnal_rmse[:,rr]
        data2[:,4] = cmaq_fc_bm3_pm25_con_val_diurnal_rmse[:,rr]
        data3 = np.full([24,5], np.nan)
        data3[:,0] = cams_ra_raw_pm25_con_val_diurnal_mbe[:,rr]
        data3[:,1] = cams_fc_raw_pm25_con_val_diurnal_mbe[:,rr]
        data3[:,2] = cmaq_fc_raw_pm25_con_val_diurnal_mbe[:,rr]
        data3[:,3] = cams_fc_bm3_pm25_con_val_diurnal_mbe[:,rr]
        data3[:,4] = cmaq_fc_bm3_pm25_con_val_diurnal_mbe[:,rr]
        plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
            ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
            xlab_diurnal, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
            colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

    ## Region 99 (outside U.S.)
    fname = plot_dir.joinpath('diurnal_o3_canmex_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Diurnal-Average Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([24,6], np.nan)
    data1[:,0] = cams_ra_raw_o3_con_val_diurnal_mean[:,11]
    data1[:,1] = cams_fc_raw_o3_con_val_diurnal_mean[:,11]
    data1[:,2] = cmaq_fc_raw_o3_con_val_diurnal_mean[:,11]
    data1[:,3] = cams_fc_bm3_o3_con_val_diurnal_mean[:,11]
    data1[:,4] = cmaq_fc_bm3_o3_con_val_diurnal_mean[:,11]
    data1[:,5] = obs_o3_con_val_diurnal_mean[:,11]
    data2 = np.full([24,5], np.nan)
    data2[:,0] = cams_ra_raw_o3_con_val_diurnal_rmse[:,11]
    data2[:,1] = cams_fc_raw_o3_con_val_diurnal_rmse[:,11]
    data2[:,2] = cmaq_fc_raw_o3_con_val_diurnal_rmse[:,11]
    data2[:,3] = cams_fc_bm3_o3_con_val_diurnal_rmse[:,11]
    data2[:,4] = cmaq_fc_bm3_o3_con_val_diurnal_rmse[:,11]
    data3 = np.full([24,5], np.nan)
    data3[:,0] = cams_ra_raw_o3_con_val_diurnal_mbe[:,11]
    data3[:,1] = cams_fc_raw_o3_con_val_diurnal_mbe[:,11]
    data3[:,2] = cmaq_fc_raw_o3_con_val_diurnal_mbe[:,11]
    data3[:,3] = cams_fc_bm3_o3_con_val_diurnal_mbe[:,11]
    data3[:,4] = cmaq_fc_bm3_o3_con_val_diurnal_mbe[:,11]
    plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
        ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
        xlab_diurnal, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
        colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

    fname = plot_dir.joinpath('diurnal_pm25_canmex_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Diurnal-Average Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([24,6], np.nan)
    data1[:,0] = cams_ra_raw_pm25_con_val_diurnal_mean[:,11]
    data1[:,1] = cams_fc_raw_pm25_con_val_diurnal_mean[:,11]
    data1[:,2] = cmaq_fc_raw_pm25_con_val_diurnal_mean[:,11]
    data1[:,3] = cams_fc_bm3_pm25_con_val_diurnal_mean[:,11]
    data1[:,4] = cmaq_fc_bm3_pm25_con_val_diurnal_mean[:,11]
    data1[:,5] = obs_pm25_con_val_diurnal_mean[:,11]
    data2 = np.full([24,5], np.nan)
    data2[:,0] = cams_ra_raw_pm25_con_val_diurnal_rmse[:,11]
    data2[:,1] = cams_fc_raw_pm25_con_val_diurnal_rmse[:,11]
    data2[:,2] = cmaq_fc_raw_pm25_con_val_diurnal_rmse[:,11]
    data2[:,3] = cams_fc_bm3_pm25_con_val_diurnal_rmse[:,11]
    data2[:,4] = cmaq_fc_bm3_pm25_con_val_diurnal_rmse[:,11]
    data3 = np.full([24,5], np.nan)
    data3[:,0] = cams_ra_raw_pm25_con_val_diurnal_mbe[:,11]
    data3[:,1] = cams_fc_raw_pm25_con_val_diurnal_mbe[:,11]
    data3[:,2] = cmaq_fc_raw_pm25_con_val_diurnal_mbe[:,11]
    data3[:,3] = cams_fc_bm3_pm25_con_val_diurnal_mbe[:,11]
    data3[:,4] = cmaq_fc_bm3_pm25_con_val_diurnal_mbe[:,11]
    plot_stat_xy_3panel_horz(fname, data1, data2, data3, xarr_diurnal,
        ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
        xlab_diurnal, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
        colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

    ## Third, make the daily-average time series plots
    ## Start by making the plots for all stations, excluding the bias-corrected forecasts
    ## CONUS
    fname = plot_dir.joinpath('daily_o3_conus_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Daily-Average Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'All Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_days,4], np.nan)
    data1[:,0] = cams_ra_raw_o3_con_all_daily_mean[:,0]
    data1[:,1] = cams_fc_raw_o3_con_all_daily_mean[:,0]
    data1[:,2] = cmaq_fc_raw_o3_con_all_daily_mean[:,0]
    data1[:,3] = obs_o3_con_all_daily_mean[:,0]
    data2 = np.full([n_valid_days,3], np.nan)
    data2[:,0] = cams_ra_raw_o3_con_all_daily_rmse[:,0]
    data2[:,1] = cams_fc_raw_o3_con_all_daily_rmse[:,0]
    data2[:,2] = cmaq_fc_raw_o3_con_all_daily_rmse[:,0]
    data3 = np.full([n_valid_days,3], np.nan)
    data3[:,0] = cams_ra_raw_o3_con_all_daily_mbe[:,0]
    data3[:,1] = cams_fc_raw_o3_con_all_daily_mbe[:,0]
    data3[:,2] = cmaq_fc_raw_o3_con_all_daily_mbe[:,0]
    plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
        ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
        xlab_daily, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
        colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    fname = plot_dir.joinpath('daily_pm25_conus_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Daily-Average Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'All Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_days,4], np.nan)
    data1[:,0] = cams_ra_raw_pm25_con_all_daily_mean[:,0]
    data1[:,1] = cams_fc_raw_pm25_con_all_daily_mean[:,0]
    data1[:,2] = cmaq_fc_raw_pm25_con_all_daily_mean[:,0]
    data1[:,3] = obs_pm25_con_all_daily_mean[:,0]
    data2 = np.full([n_valid_days,3], np.nan)
    data2[:,0] = cams_ra_raw_pm25_con_all_daily_rmse[:,0]
    data2[:,1] = cams_fc_raw_pm25_con_all_daily_rmse[:,0]
    data2[:,2] = cmaq_fc_raw_pm25_con_all_daily_rmse[:,0]
    data3 = np.full([n_valid_days,3], np.nan)
    data3[:,0] = cams_ra_raw_pm25_con_all_daily_mbe[:,0]
    data3[:,1] = cams_fc_raw_pm25_con_all_daily_mbe[:,0]
    data3[:,2] = cmaq_fc_raw_pm25_con_all_daily_mbe[:,0]
    plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
        ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
        xlab_daily, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
        colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    ## Regions 1-10
    for rr in range(1,11):
        reg_str = 'reg'+str(rr).zfill(2)
        fname = plot_dir.joinpath('daily_o3_'+reg_str+'_3panel_'+date_range_file+'_all.'+plot_type)
        suptitle_prefix = 'Daily-Average Stats for '+mpl_o3+' '+em_dash+' '
        suptitle_suffix = 'All Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([n_valid_days,4], np.nan)
        data1[:,0] = cams_ra_raw_o3_con_all_daily_mean[:,rr]
        data1[:,1] = cams_fc_raw_o3_con_all_daily_mean[:,rr]
        data1[:,2] = cmaq_fc_raw_o3_con_all_daily_mean[:,rr]
        data1[:,3] = obs_o3_con_all_daily_mean[:,rr]
        data2 = np.full([n_valid_days,3], np.nan)
        data2[:,0] = cams_ra_raw_o3_con_all_daily_rmse[:,rr]
        data2[:,1] = cams_fc_raw_o3_con_all_daily_rmse[:,rr]
        data2[:,2] = cmaq_fc_raw_o3_con_all_daily_rmse[:,rr]
        data3 = np.full([n_valid_days,3], np.nan)
        data3[:,0] = cams_ra_raw_o3_con_all_daily_mbe[:,rr]
        data3[:,1] = cams_fc_raw_o3_con_all_daily_mbe[:,rr]
        data3[:,2] = cmaq_fc_raw_o3_con_all_daily_mbe[:,rr]
        plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
            ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
            xlab_daily, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
            colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

        fname = plot_dir.joinpath('daily_pm25_'+reg_str+'_3panel_'+date_range_file+'_all.'+plot_type)
        suptitle_prefix = 'Daily-Average Stats for '+mpl_pm25+' '+em_dash+' '
        suptitle_suffix = 'All Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([n_valid_days,4], np.nan)
        data1[:,0] = cams_ra_raw_pm25_con_all_daily_mean[:,rr]
        data1[:,1] = cams_fc_raw_pm25_con_all_daily_mean[:,rr]
        data1[:,2] = cmaq_fc_raw_pm25_con_all_daily_mean[:,rr]
        data1[:,3] = obs_pm25_con_all_daily_mean[:,rr]
        data2 = np.full([n_valid_days,3], np.nan)
        data2[:,0] = cams_ra_raw_pm25_con_all_daily_rmse[:,rr]
        data2[:,1] = cams_fc_raw_pm25_con_all_daily_rmse[:,rr]
        data2[:,2] = cmaq_fc_raw_pm25_con_all_daily_rmse[:,rr]
        data3 = np.full([n_valid_days,3], np.nan)
        data3[:,0] = cams_ra_raw_pm25_con_all_daily_mbe[:,rr]
        data3[:,1] = cams_fc_raw_pm25_con_all_daily_mbe[:,rr]
        data3[:,2] = cmaq_fc_raw_pm25_con_all_daily_mbe[:,rr]
        plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
            ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
            xlab_daily, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
            colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    ## Region 99 (outside U.S.)
    fname = plot_dir.joinpath('daily_o3_canmex_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Daily-Average Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'All Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_days,4], np.nan)
    data1[:,0] = cams_ra_raw_o3_con_all_daily_mean[:,11]
    data1[:,1] = cams_fc_raw_o3_con_all_daily_mean[:,11]
    data1[:,2] = cmaq_fc_raw_o3_con_all_daily_mean[:,11]
    data1[:,3] = obs_o3_con_all_daily_mean[:,11]
    data2 = np.full([n_valid_days,3], np.nan)
    data2[:,0] = cams_ra_raw_o3_con_all_daily_rmse[:,11]
    data2[:,1] = cams_fc_raw_o3_con_all_daily_rmse[:,11]
    data2[:,2] = cmaq_fc_raw_o3_con_all_daily_rmse[:,11]
    data3 = np.full([n_valid_days,3], np.nan)
    data3[:,0] = cams_ra_raw_o3_con_all_daily_mbe[:,11]
    data3[:,1] = cams_fc_raw_o3_con_all_daily_mbe[:,11]
    data3[:,2] = cmaq_fc_raw_o3_con_all_daily_mbe[:,11]
    plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
        ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
        xlab_daily, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
        colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    fname = plot_dir.joinpath('daily_pm25_canmex_3panel_'+date_range_file+'_all.'+plot_type)
    suptitle_prefix = 'Daily-Average Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'All Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_days,4], np.nan)
    data1[:,0] = cams_ra_raw_pm25_con_all_daily_mean[:,11]
    data1[:,1] = cams_fc_raw_pm25_con_all_daily_mean[:,11]
    data1[:,2] = cmaq_fc_raw_pm25_con_all_daily_mean[:,11]
    data1[:,3] = obs_pm25_con_all_daily_mean[:,11]
    data2 = np.full([n_valid_days,3], np.nan)
    data2[:,0] = cams_ra_raw_pm25_con_all_daily_rmse[:,11]
    data2[:,1] = cams_fc_raw_pm25_con_all_daily_rmse[:,11]
    data2[:,2] = cmaq_fc_raw_pm25_con_all_daily_rmse[:,11]
    data3 = np.full([n_valid_days,3], np.nan)
    data3[:,0] = cams_ra_raw_pm25_con_all_daily_mbe[:,11]
    data3[:,1] = cams_fc_raw_pm25_con_all_daily_mbe[:,11]
    data3[:,2] = cmaq_fc_raw_pm25_con_all_daily_mbe[:,11]
    plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
        ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
        xlab_daily, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
        colors_raw, dashes_raw, markers_raw, lglab_raw_obs)

    ## Now make the plots for validation stations, including the bias-corrected forecasts
    ## CONUS
    fname = plot_dir.joinpath('daily_o3_conus_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Daily-Average Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_days,6], np.nan)
    data1[:,0] = cams_ra_raw_o3_con_val_daily_mean[:,0]
    data1[:,1] = cams_fc_raw_o3_con_val_daily_mean[:,0]
    data1[:,2] = cmaq_fc_raw_o3_con_val_daily_mean[:,0]
    data1[:,3] = cams_fc_bm3_o3_con_val_daily_mean[:,0]
    data1[:,4] = cmaq_fc_bm3_o3_con_val_daily_mean[:,0]
    data1[:,5] = obs_o3_con_val_daily_mean[:,0]
    data2 = np.full([n_valid_days,5], np.nan)
    data2[:,0] = cams_ra_raw_o3_con_val_daily_rmse[:,0]
    data2[:,1] = cams_fc_raw_o3_con_val_daily_rmse[:,0]
    data2[:,2] = cmaq_fc_raw_o3_con_val_daily_rmse[:,0]
    data2[:,3] = cams_fc_bm3_o3_con_val_daily_rmse[:,0]
    data2[:,4] = cmaq_fc_bm3_o3_con_val_daily_rmse[:,0]
    data3 = np.full([n_valid_days,5], np.nan)
    data3[:,0] = cams_ra_raw_o3_con_val_daily_mbe[:,0]
    data3[:,1] = cams_fc_raw_o3_con_val_daily_mbe[:,0]
    data3[:,2] = cmaq_fc_raw_o3_con_val_daily_mbe[:,0]
    data3[:,3] = cams_fc_bm3_o3_con_val_daily_mbe[:,0]
    data3[:,4] = cmaq_fc_bm3_o3_con_val_daily_mbe[:,0]
    plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
        ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
        xlab_daily, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
        colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

    fname = plot_dir.joinpath('daily_pm25_conus_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Daily-Average Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in CONUS (All Regions)'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_days,6], np.nan)
    data1[:,0] = cams_ra_raw_pm25_con_val_daily_mean[:,0]
    data1[:,1] = cams_fc_raw_pm25_con_val_daily_mean[:,0]
    data1[:,2] = cmaq_fc_raw_pm25_con_val_daily_mean[:,0]
    data1[:,3] = cams_fc_bm3_pm25_con_val_daily_mean[:,0]
    data1[:,4] = cmaq_fc_bm3_pm25_con_val_daily_mean[:,0]
    data1[:,5] = obs_pm25_con_val_daily_mean[:,0]
    data2 = np.full([n_valid_days,5], np.nan)
    data2[:,0] = cams_ra_raw_pm25_con_val_daily_rmse[:,0]
    data2[:,1] = cams_fc_raw_pm25_con_val_daily_rmse[:,0]
    data2[:,2] = cmaq_fc_raw_pm25_con_val_daily_rmse[:,0]
    data2[:,3] = cams_fc_bm3_pm25_con_val_daily_rmse[:,0]
    data2[:,4] = cmaq_fc_bm3_pm25_con_val_daily_rmse[:,0]
    data3 = np.full([n_valid_days,5], np.nan)
    data3[:,0] = cams_ra_raw_pm25_con_val_daily_mbe[:,0]
    data3[:,1] = cams_fc_raw_pm25_con_val_daily_mbe[:,0]
    data3[:,2] = cmaq_fc_raw_pm25_con_val_daily_mbe[:,0]
    data3[:,3] = cams_fc_bm3_pm25_con_val_daily_mbe[:,0]
    data3[:,4] = cmaq_fc_bm3_pm25_con_val_daily_mbe[:,0]
    plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
        ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
        xlab_daily, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
        colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

    ## Regions 1-10
    for rr in range(1,11):
        reg_str = 'reg'+str(rr).zfill(2)
        fname = plot_dir.joinpath('daily_o3_'+reg_str+'_3panel_'+date_range_file+'_val.'+plot_type)
        suptitle_prefix = 'Daily-Average Stats for '+mpl_o3+' '+em_dash+' '
        suptitle_suffix = 'Validation Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([n_valid_days,6], np.nan)
        data1[:,0] = cams_ra_raw_o3_con_val_daily_mean[:,rr]
        data1[:,1] = cams_fc_raw_o3_con_val_daily_mean[:,rr]
        data1[:,2] = cmaq_fc_raw_o3_con_val_daily_mean[:,rr]
        data1[:,3] = cams_fc_bm3_o3_con_val_daily_mean[:,rr]
        data1[:,4] = cmaq_fc_bm3_o3_con_val_daily_mean[:,rr]
        data1[:,5] = obs_o3_con_val_daily_mean[:,rr]
        data2 = np.full([n_valid_days,5], np.nan)
        data2[:,0] = cams_ra_raw_o3_con_val_daily_rmse[:,rr]
        data2[:,1] = cams_fc_raw_o3_con_val_daily_rmse[:,rr]
        data2[:,2] = cmaq_fc_raw_o3_con_val_daily_rmse[:,rr]
        data2[:,3] = cams_fc_bm3_o3_con_val_daily_rmse[:,rr]
        data2[:,4] = cmaq_fc_bm3_o3_con_val_daily_rmse[:,rr]
        data3 = np.full([n_valid_days,5], np.nan)
        data3[:,0] = cams_ra_raw_o3_con_val_daily_mbe[:,rr]
        data3[:,1] = cams_fc_raw_o3_con_val_daily_mbe[:,rr]
        data3[:,2] = cmaq_fc_raw_o3_con_val_daily_mbe[:,rr]
        data3[:,3] = cams_fc_bm3_o3_con_val_daily_mbe[:,rr]
        data3[:,4] = cmaq_fc_bm3_o3_con_val_daily_mbe[:,rr]
        plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
            ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
            xlab_daily, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
            colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

        fname = plot_dir.joinpath('daily_pm25_'+reg_str+'_3panel_'+date_range_file+'_val.'+plot_type)
        suptitle_prefix = 'Daily-Average Stats for '+mpl_pm25+' '+em_dash+' '
        suptitle_suffix = 'Validation Stations in EPA Region '+str(rr)
        suptitle = suptitle_prefix+suptitle_suffix
        data1 = np.full([n_valid_days,6], np.nan)
        data1[:,0] = cams_ra_raw_pm25_con_val_daily_mean[:,rr]
        data1[:,1] = cams_fc_raw_pm25_con_val_daily_mean[:,rr]
        data1[:,2] = cmaq_fc_raw_pm25_con_val_daily_mean[:,rr]
        data1[:,3] = cams_fc_bm3_pm25_con_val_daily_mean[:,rr]
        data1[:,4] = cmaq_fc_bm3_pm25_con_val_daily_mean[:,rr]
        data1[:,5] = obs_pm25_con_val_daily_mean[:,rr]
        data2 = np.full([n_valid_days,5], np.nan)
        data2[:,0] = cams_ra_raw_pm25_con_val_daily_rmse[:,rr]
        data2[:,1] = cams_fc_raw_pm25_con_val_daily_rmse[:,rr]
        data2[:,2] = cmaq_fc_raw_pm25_con_val_daily_rmse[:,rr]
        data2[:,3] = cams_fc_bm3_pm25_con_val_daily_rmse[:,rr]
        data2[:,4] = cmaq_fc_bm3_pm25_con_val_daily_rmse[:,rr]
        data3 = np.full([n_valid_days,5], np.nan)
        data3[:,0] = cams_ra_raw_pm25_con_val_daily_mbe[:,rr]
        data3[:,1] = cams_fc_raw_pm25_con_val_daily_mbe[:,rr]
        data3[:,2] = cmaq_fc_raw_pm25_con_val_daily_mbe[:,rr]
        data3[:,3] = cams_fc_bm3_pm25_con_val_daily_mbe[:,rr]
        data3[:,4] = cmaq_fc_bm3_pm25_con_val_daily_mbe[:,rr]
        plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
            ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
            xlab_daily, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
            colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

    ## Region 99 (outside U.S.)
    fname = plot_dir.joinpath('daily_o3_canmex_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Daily-Average Stats for '+mpl_o3+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_days,6], np.nan)
    data1[:,0] = cams_ra_raw_o3_con_val_daily_mean[:,11]
    data1[:,1] = cams_fc_raw_o3_con_val_daily_mean[:,11]
    data1[:,2] = cmaq_fc_raw_o3_con_val_daily_mean[:,11]
    data1[:,3] = cams_fc_bm3_o3_con_val_daily_mean[:,11]
    data1[:,4] = cmaq_fc_bm3_o3_con_val_daily_mean[:,11]
    data1[:,5] = obs_o3_con_val_daily_mean[:,11]
    data2 = np.full([n_valid_days,5], np.nan)
    data2[:,0] = cams_ra_raw_o3_con_val_daily_rmse[:,11]
    data2[:,1] = cams_fc_raw_o3_con_val_daily_rmse[:,11]
    data2[:,2] = cmaq_fc_raw_o3_con_val_daily_rmse[:,11]
    data2[:,3] = cams_fc_bm3_o3_con_val_daily_rmse[:,11]
    data2[:,4] = cmaq_fc_bm3_o3_con_val_daily_rmse[:,11]
    data3 = np.full([n_valid_days,5], np.nan)
    data3[:,0] = cams_ra_raw_o3_con_val_daily_mbe[:,11]
    data3[:,1] = cams_fc_raw_o3_con_val_daily_mbe[:,11]
    data3[:,2] = cmaq_fc_raw_o3_con_val_daily_mbe[:,11]
    data3[:,3] = cams_fc_bm3_o3_con_val_daily_mbe[:,11]
    data3[:,4] = cmaq_fc_bm3_o3_con_val_daily_mbe[:,11]
    plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
        ymin_o3_mean, ymax_o3_mean, ymin_o3_rmse, ymax_o3_rmse, ymin_o3_mbe, ymax_o3_mbe,
        xlab_daily, ylab_o3_mean, ylab_o3_rmse, ylab_o3_mbe, suptitle, title1, title2, title3,
        colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)

    fname = plot_dir.joinpath('daily_pm25_canmex_3panel_'+date_range_file+'_val.'+plot_type)
    suptitle_prefix = 'Daily-Average Stats for '+mpl_pm25+' '+em_dash+' '
    suptitle_suffix = 'Validation Stations in Canada & Mexico'
    suptitle = suptitle_prefix+suptitle_suffix
    data1 = np.full([n_valid_days,6], np.nan)
    data1[:,0] = cams_ra_raw_pm25_con_val_daily_mean[:,11]
    data1[:,1] = cams_fc_raw_pm25_con_val_daily_mean[:,11]
    data1[:,2] = cmaq_fc_raw_pm25_con_val_daily_mean[:,11]
    data1[:,3] = cams_fc_bm3_pm25_con_val_daily_mean[:,11]
    data1[:,4] = cmaq_fc_bm3_pm25_con_val_daily_mean[:,11]
    data1[:,5] = obs_pm25_con_val_daily_mean[:,11]
    data2 = np.full([n_valid_days,5], np.nan)
    data2[:,0] = cams_ra_raw_pm25_con_val_daily_rmse[:,11]
    data2[:,1] = cams_fc_raw_pm25_con_val_daily_rmse[:,11]
    data2[:,2] = cmaq_fc_raw_pm25_con_val_daily_rmse[:,11]
    data2[:,3] = cams_fc_bm3_pm25_con_val_daily_rmse[:,11]
    data2[:,4] = cmaq_fc_bm3_pm25_con_val_daily_rmse[:,11]
    data3 = np.full([n_valid_days,5], np.nan)
    data3[:,0] = cams_ra_raw_pm25_con_val_daily_mbe[:,11]
    data3[:,1] = cams_fc_raw_pm25_con_val_daily_mbe[:,11]
    data3[:,2] = cmaq_fc_raw_pm25_con_val_daily_mbe[:,11]
    data3[:,3] = cams_fc_bm3_pm25_con_val_daily_mbe[:,11]
    data3[:,4] = cmaq_fc_bm3_pm25_con_val_daily_mbe[:,11]
    plot_stat_vs_time_3panel_vert(fname, data1, data2, data3, dt_dates_valid,
        ymin_pm25_mean, ymax_pm25_mean, ymin_pm25_rmse, ymax_pm25_rmse, ymin_pm25_mbe, ymax_pm25_mbe,
        xlab_daily, ylab_pm25_mean, ylab_pm25_rmse, ylab_pm25_mbe, suptitle, title1, title2, title3,
        colors_raw_bm3, dashes_raw_bm3, markers_raw_bm3, lglab_raw_bm3_obs)


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
