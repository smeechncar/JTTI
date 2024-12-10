'''
get_cams_eac4_pm25_o3.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 1 Nov 2021

This script downloads CAMS Global Reanalysis (EAC4) variables (PM2.5 and surface ozone) using the ECMWF ADS API.
https://ads.atmosphere.copernicus.eu/api-how-to
https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-reanalysis-eac4
'''

import cdsapi
import pathlib
import sys
import os
import datetime as dt
import pandas as pd

beg_date = '2022-01-01'
end_date = '2022-12-31'
hours = ['00', '03', '06', '09', '12', '15', '18', '21']
#variables = ['ozone','particulate_matter_2.5um']
#variables = ['particulate_matter_2.5um']
#variables = ['ozone']

fmt_date = '%Y-%m-%d'
beg_dt = pd.to_datetime(beg_date, format=fmt_date)
end_dt = pd.to_datetime(end_date, format=fmt_date)
all_dt = pd.date_range(start=beg_dt, end=end_dt, freq='D')
n_dates = len(all_dt)
n_hours = len(hours)

#out_dir_parent = pathlib.Path('/','glade','p','ral','wsap','jaredlee','NOAA_CMAQ_AnEn','cams','eac4')
#out_dir_parent = pathlib.Path('/','glade','derecho','scratch','jaredlee','NOAA_CMAQ_AnEn','cams','eac4')
out_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','eac4')

## Loop over dates
for dd in range(n_dates):
	this_dt = all_dt[dd]
	this_date = this_dt.strftime(fmt_date)
	this_yr = this_dt.strftime('%Y')
	this_mo = this_dt.strftime('%m')

	out_dir = out_dir_parent.joinpath(this_yr, this_mo)
	pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
	print('Processing date '+this_date)

	## Loop over hours
	for hh in range(n_hours):
		this_hr = hours[hh]
		out_path_pm25 = out_dir.joinpath('pm2.5_'+this_date+'_'+this_hr+'00.nc')
		out_path_o3 = out_dir.joinpath('o3sfc_'+this_date+'_'+this_hr+'00.nc')

		c = cdsapi.Client()
		c.retrieve(
			'cams-global-reanalysis-eac4',
			{
				'date': this_date,
				'format': 'netcdf',
				'variable': 'particulate_matter_2.5um',
				'model_level': '60',
				'time': this_hr+':00',
				'area': [55, -135, 20, -60],
			},
			str(out_path_pm25))

		c.retrieve(
			'cams-global-reanalysis-eac4',
			{
				'date': this_date,
				'format': 'netcdf',
				'variable': 'ozone',
				'model_level': '60',
				'time': this_hr+':00',
				'area': [55, -135, 20, -60]
			},
			str(out_path_o3))
