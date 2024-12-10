'''
randomly_withhold_obs.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 6 Jun 2022

This script demonstrates reading in a file that contains AirNOW observations.
There is one row per station, and columns with metadata plus 24 hourly observations.
For each hour/column, the script identifies the non-missing observations, randomly flags a fraction
of them for withholding, then writes out a file with the mask for repeatability and validation.
The output file has the same format, but an M/U/W for missing/used/withheld.
'''

import sys
import pathlib
import datetime as dt
import numpy as np
import pandas as pd

withhold_fraction = 0.333

data_dir = pathlib.Path('/','glade','scratch','jaredlee','NOAA_CMAQ_AnEn','obs')

file_in = data_dir.joinpath('airnow_20200915_pm25.csv')
file_out = data_dir.joinpath('airnow_20200915_pm25_mask.csv')

## Read in the file to a pandas dataframe
df_in = pd.read_csv(file_in, header=None)
#print(df_in)

## Make a copy of the dataframe to modify later
df_out = df_in.copy()

## Col 0: Station ID
## Col 1: CMAQ i-grid point
## Col 2: CMAQ j-grid point
## Col 3: Station longitude (deg E)
## Col 4: Station latitude (deg N)
## Cols 5–28: Hourly data values
sid = df_in[0].to_numpy()

## Loop through hours
n_hours = 24
for hh in range(n_hours):
	## Extract data values
	vals = df_in[hh+5].to_numpy()
	n_vals = len(vals)

	## Get indices of missing and valid values
	inds_missing = np.where(vals == -999.0)[0]
	n_missing = len(inds_missing)
	inds_valid = np.where(vals != -999.0)[0]
	n_valid = len(inds_valid)

	## How many of the valid obs should be withheld for validation?
	n_withheld = int(n_valid * withhold_fraction)

	## Randomly select without replacement the indices to be withheld. Use all other indices.
	inds_withheld = np.sort(np.random.choice(range(n_valid), n_withheld, replace=False))
	inds_used = np.delete(range(n_valid), inds_withheld)
	inds_valid_withheld = inds_valid[inds_withheld]
	inds_valid_used     = inds_valid[inds_used]
	n_used = len(inds_valid_used)

	## Define a three-value masked array: 'U' = used, 'W' = withheld, 'M' = missing
	mask_col = np.full(n_vals, 'U')
	mask_col[inds_valid_withheld] = 'W'
	mask_col[inds_missing] = 'M'

	## Now replace the column in the output dataframe with the three-value mask
	df_out[hh+5] = mask_col

## Now write out the modified data array to a csv file
print('Writing '+str(file_out))
df_out.to_csv(file_out)


