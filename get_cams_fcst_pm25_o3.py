'''
get_cams_fcst_pm25_o3.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 12 Nov 2021

This script downloads CAMS Global Atmospheric Composition Forecasts variables (PM2.5 and surface ozone) using the ECMWF ADS API.
https://ads.atmosphere.copernicus.eu/api-how-to
https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-atmospheric-composition-forecasts
'''

import cdsapi
import pathlib
import sys
import os
import datetime as dt
import pandas as pd

beg_date = '2022-12-01'
end_date = '2022-12-31'
cycles = ['00', '12']
#cycles = ['12']
get_pm25 = True
get_o3 = True
#lead_hrs_pm25 = ['04', '05', '06', '07', '08', '09', '10', '11']
lead_hrs_pm25 = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11']
lead_hrs_o3 = ['00', '03', '06', '09']
#lead_hrs_o3 = ['09']
#variables = ['ozone','particulate_matter_2.5um']
#variables = ['particulate_matter_2.5um']
#variables = ['ozone']

fmt_date = '%Y-%m-%d'
beg_dt = pd.to_datetime(beg_date, format=fmt_date)
end_dt = pd.to_datetime(end_date, format=fmt_date)
all_dt = pd.date_range(start=beg_dt, end=end_dt, freq='D')
n_dates = len(all_dt)
n_cycles = len(cycles)
n_leads_pm25 = len(lead_hrs_pm25)
n_leads_o3 = len(lead_hrs_o3)

out_dir_parent = pathlib.Path('/','glade','campaign','ral','nsap','JTTI','cams','fcst')

## Loop over dates
for dd in range(n_dates):
	this_dt = all_dt[dd]
	this_date = this_dt.strftime(fmt_date)
	this_yr = this_dt.strftime('%Y')
	this_mo = this_dt.strftime('%m')
	this_dy = this_dt.strftime('%d')

	## Loop over cycles
	for cc in range(n_cycles):
		this_hr = cycles[cc]
		init_dt = pd.to_datetime(this_date+'_'+this_hr, format=fmt_date+'_%H')

		out_dir = out_dir_parent.joinpath(this_yr, this_mo, this_yr+this_mo+this_dy+'_'+this_hr)
		pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
		print('Processing date '+this_date)

		## Loop over lead times
		if get_pm25:
			for ff in range(n_leads_pm25):
				this_lead = lead_hrs_pm25[ff]
				this_lead_num = str(int(this_lead))
				lead_dt = init_dt + dt.timedelta(hours=int(this_lead))
				lead_dt_str = lead_dt.strftime(fmt_date+'_%H00')
#				out_path_pm25 = out_dir.joinpath('pm2.5_'+this_date+'_'+this_hr+'00_f+'+this_lead+'.nc')
				out_path_pm25 = out_dir.joinpath('pm2.5_'+lead_dt_str+'.nc')
				c = cdsapi.Client()
				c.retrieve(
					'cams-global-atmospheric-composition-forecasts',
					{
						'date': this_date,
						'type': 'forecast',
						'format': 'netcdf',
						'variable': 'particulate_matter_2.5um',
						'time': this_hr+':00',
						'area': [55, -135, 20, -60],
						'leadtime_hour': this_lead_num,
					},
					str(out_path_pm25))

		if get_o3:
			for ff in range(n_leads_o3):
				this_lead = lead_hrs_o3[ff]
				this_lead_num = str(int(this_lead))
				lead_dt = init_dt + dt.timedelta(hours=int(this_lead))
				lead_dt_str = lead_dt.strftime(fmt_date+'_%H00')
#				out_path_o3 = out_dir.joinpath('o3sfc_'+this_date+'_'+this_hr+'00_f+'+this_lead+'.nc')
				out_path_o3 = out_dir.joinpath('o3sfc_'+lead_dt_str+'.nc')
				c = cdsapi.Client()
				c.retrieve(
					'cams-global-atmospheric-composition-forecasts',
					{
						'date': this_date,
						'type': 'forecast',
						'format': 'netcdf',
						'variable': 'ozone',
						'model_level': '137',
						'time': this_hr+':00',
						'area': [55, -135, 20, -60],
						'leadtime_hour': this_lead_num,
					},
					str(out_path_o3))
