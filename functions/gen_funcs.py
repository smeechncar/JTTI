'''
gen_funcs.py

Written by: Jared A. Lee
            jaredlee@ucar.edu

This file contains several general-use functions.

calc_corrcoef
	-- Added by JAL on 19 Sep 2019
	-- This function calculates the Pearson-r correlation coefficient between two time series.
	-- Inputs:
		- fcst: a 1-dimensional numpy array of forecast/predicted values
		- obs: a 1-dimensional numpy array of observations
	-- Output:
		- corrcoef: a scalar value of the Pearson r

calc_fge
	-- Added by JAL on 18 Feb 2020
	-- This function calculates the fractional gross error (FGE) from forecast and observation arrays.
		The FGE ranges between 0 (best) and 2 (worst).
	-- Reference: Remy et al. (2019), https://doi.org/10.5194/gmd-12-4627-2019, p. 4648.
	-- Inputs:
		- fcst: an N-dimensional numpy array of forecast/predicted values
		- obs: an N-dimensional numpy array of observations
	-- Output:
		- fge: a numpy array containing the FGE

calc_mae
	-- Added by JAL on 3 Jun 2019
	-- This function calculates the mean absolute error (MAE) from forecast and observation arrays.
	-- Inputs:
		- fcst: an N-dimensional numpy array of forecast/predicted values
		- obs: an N-dimensional numpy array of observations
	-- Output:
		- mae: a numpy array containing the MAE

calc_mbe / calc_me
	-- Added by JAL on 3 Jun 2019 (calc_me)
	-- Added by JAL on 14 Apr 2020 (calc_mbe)
	-- This function calculates the mean (bias) error (ME or MBE) from forecast and observation arrays.
	-- Inputs:
		- fcst: an N-dimensional numpy array of forecast/predicted values
		- obs: an N-dimensional numpy array of observations
	-- Output:
		- me (mbe): a numpy array containing the ME (MBE)

calc_mape
	-- (Re-)Added by JAL on 17 Apr 2020 (after having accidentally deleted it previously)
	-- This function calculates the mean absolute percentage error (MAPE) from forecast and observation arrays.
	-- Inputs:
		- fcst: an N-dimensional numpy array of forecast/predicted values
		- obs: an N-dimensional numpy array of observations
	-- Output:
		- mape: a numpy array containing the MAPE

calc_mpe
	-- (Re-)Added by JAL on 17 Apr 2020 (after having accidentally deleted it previously)
	-- This function calculates the mean percentage error (MPE) from forecast and observation arrays.
	-- Inputs:
		- fcst: an N-dimensional numpy array of forecast/predicted values
		- obs: an N-dimensional numpy array of observations
	-- Output:
		- mpe: a numpy array containing the MPE

calc_mnmb
	-- Added by JAL on 18 Feb 2020
	-- This function calculates the modified normalized mean bias (MNMB) from forecast and observation arrays.
		The MNMB ranges between -2 and 2, with 0 being best.
	-- Reference: Remy et al. (2019), https://doi.org/10.5194/gmd-12-4627-2019, p. 4648.
	-- Inputs:
		- fcst: an N-dimensional numpy array of forecast/predicted values
		- obs: an N-dimensional numpy array of observations
	-- Output:
		- mnmb: a numpy array containing the MNMB

calc_rmse
	-- Added by JAL on 3 Jun 2019
	-- This function calculates the root mean squared error (RMSE) from forecast and observation arrays.
	-- Inputs:
		- fcst: an N-dimensional numpy array of forecast/predicted values
		- obs: an N-dimensional numpy array of observations
	-- Output:
		- me: a numpy array containing the RMSE

make_ngl_xy_plot
	-- Added by JAL on 4 Jun 2019
	-- This function generates an xy plot using PyNGL.
	-- Inputs:
		- plot_type: string of file format/suffix
		- plot_file: string or Path object of file name (minus file extension)
		- x: numeric array of values on the x axis
		- y: numeric array of values on the y axis
		- res: Ngl.Resources object controlling the main plot
		- lg_txt: array of strings for legend labels
		- lg_llx: float, NDC x coordinate of lower left corner of legend
		- lg_lly: float, NDC y coordinate of lower left corner of legend
		- lgres: Ngl.Resources object controlling the legend

set_ngl_sliceres
	-- Added by JAL on 1 Aug 2019
	-- This function sets some Ngl plot resources for vertical slice plots (e.g., from WRF output)
	-- Inputs: none
	-- Output: Ngl.Resources object

find_grid_inds(grid_lat, grid_lon, lat_vals, lon_vals):
	-- Added 10 Mar 2023
	-- This function finds the closest j,i index in the grid lat/lon of a provided list/array of lat/lon values
	-- Inputs:
		- grid_lat: 2D latitude array
		- grid_lon: 2D longitude array
		- lat_vals: 1D latitude array
		- lon_vals: 1D longitude array
	-- Output:
		- j, i: 1D arrays of indices

find_nearest_index
	-- Added by JAL on 8 Apr 2019
	-- This function finds the index of the nearest value in a 1-D array to a requested value.
	-- Inputs:
		- array: a 1-D numpy array
		- value: a numeric value for which the nearest-value index is desired
	-- Output:
		- idx: the index value of the array that is closest to the requested number

find_nearest_value
	-- Added by JAL on 8 Apr 2019
	-- This function finds the nearest value in a 1-D array to a requested value.
	-- Inputs:
		- array: a 1-D numpy array
		- value: a numeric value for which the nearest value is desired
	-- Output:
		- array[idx]: the value in the array that is closest to the requested number

Data_Array
	-- Added by JAL on 10 Apr 2019
	-- This class of functions allows for a new _FillValue attribute to be added to an array.
	-- Usage: x = gen_funcs.Data_Array(input_array, _FillValue={numeric value})

'''

import numpy as np
import warnings
#import Ngl
import sys

def quad_poly_func(x, a, b, c):
	return a*x**2 + b*x + c

def calc_aod550(aeronet_aod):
	from scipy.optimize import curve_fit

	## AERONET doesn't report AOD550 directly. Derive it from surrounding AOD values.
	## Follow Gueymard and Yang (2020, https://doi.org/10.1016/j.atmosenv.2019.117216)
	## Use: log(tau_lambda) = a0 + a1*log(lambda) + a2*(log(lambda))^2
	
	## Alternatively, we could leverage the AE (but we don't do that here):
	## tau_550 = tau_440 * (550/440)^(-alpha_440-675) or
	## tau_550 = tau_500 * (550/500)^(-alpha_440-675) or
	## tau_550 = tau_675 * (550/675)^(-alpha_440-675)

	print('Interpolating AERONET AOD data to AOD550')
	n_aeronet_times = aeronet_aod.shape[0]
	n_aeronet_sites = aeronet_aod.shape[1]
	n_aeronet_bands = aeronet_aod.shape[2]

	if n_aeronet_bands != 8:
		print('ERROR! Function calc_aod550 assumes that aeronet_aod has 8 wavelength bands:')
		print('       1640, 1020, 870, 675, 500, 440, 380, and 340 nm')
		print('       This array only has '+str(n_aeronet_bands)+' bands. Exiting!')
		sys.exit()

	aeronet_aod_550  = np.full([n_aeronet_times, n_aeronet_sites], np.nan)
	aeronet_aod_380  = aeronet_aod[:,:,6]
	aeronet_aod_440  = aeronet_aod[:,:,5]
	aeronet_aod_500  = aeronet_aod[:,:,4]
	aeronet_aod_675  = aeronet_aod[:,:,3]
	aeronet_aod_870  = aeronet_aod[:,:,2]
	aeronet_aod_1020 = aeronet_aod[:,:,1]

	## Gueymard and Yang (2020) method, using four AOD measurements to interpolate a quadratic polynomial
	for tt in range(n_aeronet_times):
		for ss in range(n_aeronet_sites):
			aod380  = aeronet_aod_380[tt,ss]
			aod440  = aeronet_aod_440[tt,ss]
			aod500  = aeronet_aod_500[tt,ss]
			aod675  = aeronet_aod_675[tt,ss]
			aod870  = aeronet_aod_870[tt,ss]
			aod1020 = aeronet_aod_1020[tt,ss]

			## Gueymard and Yang (2020) throw out dubious AOD values <=0.00 or >5.00
			## Negative values cause np.log(aod) to be inf or nan below, so they need to be removed
			if aod380 <= 0.0 or aod380 > 5.0:
				aod380 = np.nan
			if aod440 <= 0.0 or aod440 > 5.0:
				aod440 = np.nan
			if aod500 <= 0.0 or aod500 > 5.0:
				aod500 = np.nan
			if aod675 <= 0.0 or aod675 > 5.0:
				aod675 = np.nan
			if aod870 <= 0.0 or aod870 > 5.0:
				aod870 = np.nan
			if aod1020 <= 0.0 or aod1020 > 5.0:
				aod1020 = np.nan

			is380   = ~np.isnan(aod380)
			is440   = ~np.isnan(aod440)
			is500   = ~np.isnan(aod500)
			is675   = ~np.isnan(aod675)
			is870   = ~np.isnan(aod870)
			is1020  = ~np.isnan(aod1020)
			isaod   = np.array([is440, is500, is675, is870, is380, is1020])
			wl_all  = np.array([440, 500, 675, 870, 380, 1020])
			aod_all = np.array([aod440, aod500, aod675, aod870, aod380, aod1020])

			## How many good values do we have? Need at least 4. If fewer, then cycle the loop
			## These tests and conditions are explicitly mentioned in Gueymard and Yang (2020)
			num_aod = np.count_nonzero(isaod)
			if num_aod < 4:
				continue
			else:
				## Take the indices of the first four good AOD values
				inds = np.where(isaod)[0]
				wl   = wl_all[inds[0:4]]
				aod  = aod_all[inds[0:4]]
				## Exception: If 380 and 1020 are the last two elements, don't use it, and cycle the loop
				if wl[2] == 380 and wl[3] == 1020:
					continue
				ln_wl = np.log(wl)
				ln550 = np.log(550)
				lnaod = np.log(aod)

				## Find the quadratic polynomial coefficients using scipy.curve_fit
				## See: https://www.datatechnotes.com/2020/09/curve-fitting-with-curve-fit-function-in-python.html
				params, covs = curve_fit(quad_poly_func, ln_wl, lnaod)
				a2, a1, a0 = params[0], params[1], params[2]

				## Now calculate AOD550
				aeronet_aod_550[tt,ss] = np.exp(a0 + a1*ln550 + a2*(ln550**2))

	return aeronet_aod_550


def calc_corrcoef(fcst, obs):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)
		warnings.simplefilter('ignore', category=UserWarning)
		if fcst.ndim != 1:
			print('ERROR! Input fcst in function calc_corrcoef has more than 1 dimension. Unexpected behavior may result. Exiting!')
			sys.exit()
		if obs.ndim != 1:
			print('ERROR! Input obs in function calc_corrcoef has more than 1 dimension. Unexpected behavior may result. Exiting!')
			sys.exit()
		if all(np.isnan(fcst)) or all(np.isnan(obs)):
			## Do this to avoid warning messages if either of the input arrays is all missing
			corrcoef = np.nan
		else:
			if any(np.isnan(fcst)) or any(np.isnan(obs)):
				## If any missing data is in the arrays, then we need to create masked arrays
				fcst_ma  = np.ma.masked_invalid(fcst)
				obs_ma   = np.ma.masked_invalid(obs)
				corrcoef_matrix   = np.ma.corrcoef(fcst_ma, obs_ma, allow_masked=True)
			else:
				corrcoef_matrix = np.corrcoef(fcst, obs)
			corrcoef = corrcoef_matrix[0,1]
	return corrcoef

def calc_fge(fcst, obs, axis=None):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)
		fge = 2.0 * np.nanmean( np.absolute( (fcst-obs) / (fcst+obs) ), axis=axis )
	return fge

def calc_mae(fcst, obs, axis=None):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)
		mae = np.nanmean(np.absolute(fcst-obs), axis=axis)
	return mae

def calc_mbe(fcst, obs, axis=None):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)
		mbe = np.nanmean(fcst-obs, axis=axis)
	return mbe

def calc_me(fcst, obs, axis=None):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)
		me = np.nanmean(fcst-obs, axis=axis)
	return me

def calc_mape(fcst, obs, axis=None):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)
		mape = np.nanmean(np.abs((fcst-obs) / obs), axis=axis) * 100
	return mape

def calc_mpe(fcst, obs, axis=None):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)
		mpe = np.nanmean( (fcst-obs) / obs , axis=axis) * 100
	return mpe

def calc_mnmb(fcst, obs, axis=None):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)
		mnmb = 2.0 * np.nanmean( (fcst - obs) / (fcst + obs) , axis=axis)
	return mnmb

def calc_rmse(fcst, obs, axis=None):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)
		rmse = np.sqrt(np.nanmean((fcst-obs)**2, axis=axis))
	return rmse

def find_grid_inds(grid_lat, grid_lon, lat_vals, lon_vals):
	## Find closest j,i index in the grid lat/lon of a provided list/array of lat/lon values
	if len(lat_vals) != len(lon_vals):
		print('ERROR! gen_funcs.find_grid_inds: Length of lat_vals and lon_vals not the same. Exiting!')
		sys.exit()
	if grid_lat.shape != grid_lon.shape:
		print('ERROR! gen_funcs.find_grid_inds: Shape of grid_lat and grid_lon not the same. Exiting!')
		sys.exit()

	n_inds = len(lat_vals)
	j = np.full(n_inds, 0)
	i = np.full(n_inds, 0)
	for nn in range(n_inds):
		dist = abs(grid_lat-lat_vals[nn])+abs(grid_lon-lon_vals[nn])
		j[nn], i[nn] = np.unravel_index(dist.argmin(), dist.shape)
	return j, i

def find_nearest_index(array, value):
	array	= np.asarray(array)
	idx	= (np.abs(array - value)).argmin()
	return idx

def find_nearest_value(array, value):
	array	= np.asarray(array)
	idx	= (np.abs(array - value)).argmin()
	return array[idx]

def make_ngl_xy_plot(plot_type, plot_file, x, y, res, lg_txt, lg_llx, lg_lly, lgres):
	print('-- Creating plot: '+str(plot_file)+'.'+plot_type)
	wks	= Ngl.open_wks(plot_type, str(plot_file))
	plot	= Ngl.xy(wks, x, y, res)
	if res.nglDraw == False:
		Ngl.draw(plot)
	Ngl.legend_ndc(wks, len(lg_txt), lg_txt, lg_llx, lg_lly, lgres)
	if res.nglFrame == False:
		Ngl.frame(wks)
	Ngl.delete_wks(wks)

def set_ngl_sliceres():
	res = Ngl.Resources()
	res.nglDraw			= False
	res.nglFrame		= False
	res.nglMaximize	= True

	res.cnFillOn					= True
	res.cnFillPalette				= 'NCL_default'
	res.cnLinesOn					= False
	res.cnLineLabelsOn			= False
	res.cnInfoLabelOn				= False
	res.cnSpanFillPalette		= True
	res.cnFillMode					= 'RasterFill'
	res.cnRasterSmoothingOn		= False
	res.cnLevelSelectionMode	= 'ManualLevels'

	res.lbLabelBarOn			= True
	res.lbOrientation			= 'horizontal'
	res.lbTitleOn				= True
	res.lbTitlePosition		= 'Bottom'
	res.lbTitleFontHeightF	= 0.015
	res.lbLabelFontHeightF	= 0.015

	res.tiMainFontHeightF	= 0.020
	res.tiXAxisFontHeightF	= 0.015
	res.tiYAxisFontHeightF	= 0.015

	res.pmTickMarkDisplayMode	= 'Always'

	res.tmXBOn	= True
	res.tmXTOn	= False
	res.tmYLOn	= True
	res.tmYROn	= True
	res.tmXBLabelFontHeightF	= 0.008
	res.tmYLLabelFontHeightF	= 0.008

	return res


'''
def make_ngl_xy_plot_text(plot_type, plot_file, x, y, res, lg_txt, lg_llx, lg_lly, lgres, txt_strs, txt_x, txt_y, txres)
	print('-- Creating plot: '+str(plot_file)+'.'+plot_type)
	res.nglDraw		= False
	res.nglFrame	= False
	wks	= Ngl.open_wks(plot_type, str(plot_file))
	plot	= Ngl.xy(wks, x, y, res)
	Ngl.draw(plot)
	Ngl.legend_ndc(wks, len(lg_txt), lg_txt, lg_llx, lg_lly, lgres)
	txt_ids	= 
	for ii in range(len(txt_strs)):
		
	Ngl.frame(wks)
	Ngl.delete_wks(wks)
'''

class Data_Array(np.ndarray):
#  def __new__(subtype, shape, dtype=float, buffer=None, offset=0, strides=None, order=None, info=None):
#     ## Create the ndarray instance of our type, given the usual ndarray arguments.
#     ## This will call the standard ndarray constructor, but return an object of our type.
#     ## It also triggers a call to Data_Array.__array_finalize__
#     obj   = super(Data_Array, subtype).__new__(subtype, shape, dtype, buffer, offset, strides, order)
#     ## Set the new '_FillValue' attribute to the value passed
#     obj._FillValue = info
#     ## Finally, we must return the newly created object
#     return obj

   def __new__(cls, input_array, _FillValue=None):
      ## Input array is an already formed ndarray instance
      ## We first cast to be our class type
      obj   = np.asarray(input_array).view(cls)
      ## Add the new attribute to the created instance
      obj._FillValue = _FillValue
      ## Finally, we must return the newly created object
      return obj

   def __array_finalize__(self, obj):
      ## 'self' is a new object resulting from ndarray.__new__(Data_Array, ...), therefore it only has
      ## attributes that the ndarray.__new__ constructor gave it —- i.e., those of a standard ndarray.
      ## We could have gotten to the ndarray.__new__ call in 3 ways:
      ## 1. From an explicit constructor -- e.g., Data_Array():
      ##    obj is None
      ##    (we're in the middle of the Data_Array.__new__ constructor, and self._FillValue will be set when we return to Data_Array.__new__)
      if obj is None: return
      ## 2. From view casting -- e.g., arr.view(Data_Array):
      ##    obj is arr
      ##    (type(obj) can be Data_Array)
      ## 3. From new-from-template -- e.g., infoarr[:3]
      ##    type(obj) is Data_Array
      ##
      ## Note that it is here, rather than in the __new__ method, that we set the default value for _FillValue,
      ## because this method sees all creation of default objects -- with the Data_Array.__new__ constructor,
      ## but also with arr.view(Data_Array).
      self._FillValue   = getattr(obj, '_FillValue', None)
