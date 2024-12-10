"""
rename_vars.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 23 Aug 2024

This script changes the name of a variable in a series of files.
"""

import pathlib
import subprocess
import pandas as pd
import datetime as dt

data_dir = pathlib.Path('/', 'glade', 'campaign', 'ral', 'nsap', 'JTTI', 'noaa_bc_cmaq_06z')
old_var = 'pm25_anen_cmaq_bm1'
new_var = 'pm25_anen_cmaq_bm3'

beg_date = pd.to_datetime('20210101', format='%Y%m%d')
end_date = pd.to_datetime('20211231', format='%Y%m%d')
dates = pd.date_range(start=beg_date, end=end_date, freq='1D')
n_dates = len(dates)

# Loop over dates
for dd in range(n_dates):
    this_date = dates[dd]
    this_date_str = this_date.strftime('%Y%m%d')
    fname = data_dir.joinpath('bc_cmaq_06z_pm2.5_o3_06z_'+this_date_str + '.nc')
    print('Renaming variable in ' + str(fname))
    subprocess.run(['ncrename', '-h', '-O', '-v', old_var+','+new_var, str(fname)])
