# JTTI
Scripts used to create NOAA JTTI CMAQ / CAMS data



Scripts used to create JTTI data based on these steps outlined below:

#  (Ju-Hye) Reading Bufr files for Airnow Stations, generating csv files without 500 micg/m3 threshold (/glade/campaign/ral/nsap/JTTI/bcdata/airnow)
Input:
/glade/campaign/ral/nsap/JTTI/bcdata/airnow/YYYY/YYYYMMDD/b008/xx021(o3)
/glade/campaign/ral/nsap/JTTI/bcdata/airnow/YYYY/YYYYMMDD/b008/xx031(pm25)
Output: 
/glade/campaign/ral/nsap/JTTI/bcdata/airnow/csv_pm25_noQC
/glade/campaign/ral/nsap/JTTI/bcdata/airnow/csv_o3

# (Jared) Generate again netcdf filed with airnow observations splitting in training and verification (*use* and *val*)
python withhold_airnow_stations.py 20200801_00 -l 20211231_23
Inputs:
AirNow PM2.5 daily obs CSV files in: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/csv_pm25_noQC
AirNow O3 daily obs CSV files in: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/csv_o3
Outputs:
AirNow obs, training stations (sites vary each hour): /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_vary/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_use.nc
AirNow obs, validation stations (sites vary each hour): /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_vary/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_val.nc
python tally_airnow_valid_obs.py 20200801_00 -l 20211231_23
Inputs:
AirNow obs, training stations (sites vary each hour): /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_vary/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_use.nc
AirNow obs, validation stations (sites vary each hour): /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_vary/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_val.nc
Outputs:
AirNow obs, all stations for full period: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/airnow_pm2.5_o3_valid_uptime_20200801_00-20211231_23_all.nc
AirNow obs, training stations for full period: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/airnow_pm2.5_o3_valid_uptime_20200801_00-20211231_23_use.nc
AirNow obs, validation stations for full period: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/airnow_pm2.5_o3_valid_uptime_20200801_00-20211231_23_val.nc
pm2.5 stations: all = 961, use = 574, val = 387
o3 stations: all = 1192, use = 517, val = 675
python split_airnow_static_use_val.py 20200801_00 -l 20211231_23
Inputs:
AirNow obs, all stations over full period: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/airnow_pm2.5_o3_valid_uptime_20200801_00-20211231_23_all.nc
AirNow obs, training stations over full period: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/airnow_pm2.5_o3_valid_uptime_20200801_00-20211231_23_use.nc
AirNow obs, validation stations over full period: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/airnow_pm2.5_o3_valid_uptime_20200801_00-20211231_23_val.nc
Outputs:
AirNow obs, training stations (sites static every hour) /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_use.nc
AirNow obs, validation stations (sites static every hour): /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_val.nc

# (Jared) Generate files with CAMS raw (first 12 hours) and obs, CMAQ raw (first 6 hours) and obs at AirNow sites. This is the analysis time series.
python extract_cams_at_airnow_use_val.py 20200801_00 -l 20211231_23 -f
Inputs:
CAMS O3 on native 0.4° lat/lon grid time series, every 3 h: /glade/campaign/ral/nsap/JTTI/cams/fcst/YYYY/MM/YYYYMMDD_HH/o3sfc_YYYY-MM-DD_HH00.nc
CAMS PM2.5 on native 0.4° lat/lon grid time series, every 1 h: /glade/campaign/ral/nsap/JTTI/cams/fcst/YYYY/MM/YYYYMMDD_HH/pm2.5_YYYY-MM-DD_HH00.nc
AirNow obs, training stations: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_use.nc
AirNow obs, validation stations: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_val.nc
Outputs:
CAMS analysis time series and AirNow obs, training stations: /glade/campaign/ral/nsap/JTTI/cams/fcst/YYYY/MM/YYYYMMDD_HH/cams_airnow_pm2.5_o3_static_YYYYMMDD_HH00_use.nc
CAMS analysis time series and AirNow obs, validation stations: /glade/campaign/ral/nsap/JTTI/cams/fcst/YYYY/MM/YYYYMMDD_HH/cams_airnow_pm2.5_o3_static_YYYYMMDD_HH00_val.nc
python extract_cmaq_at_airnow_use_val.py 20200801_00 -l 20211231_23
Inputs:
CMAQ 12-km native grid, O3 & PM2.5: /glade/campaign/ral/nsap/JTTI/bcdata/cmaq_hourly/cmaq_airnow_pm2.5_o3_YYYYMMDDHH.nc
AirNow obs, training sites: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_use.nc
AirNow obs, validation sites: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_val.nc
Outputs:
CMAQ analysis time series and AirNow obs, training stations: /glade/campaign/ral/nsap/JTTI/bcdata/cmaq_hourly/sites_static/YYYY/MM/cmaq_airnow_pm2.5_o3_static_YYYYMMDD_HH00_use.nc
CMAQ analysis time series and AirNow obs, validation stations: /glade/campaign/ral/nsap/JTTI/bcdata/cmaq_hourly/sites_static/YYYY/MM/cmaq_airnow_pm2.5_o3_static_YYYYMMDD_HH00_val.nc

# (Ju-Hye) Generate BM3 for CAMS and CMAQ for PM2.5
Input: 
/glade/campaign/ral/nsap/JTTI/cams/fcst/YYYY/MM/YYYYMMDD_HH/cams_airnow_pm2.5_o3_static_YYYYMMDD_HH00_use.nc
/glade/campaign/ral/nsap/JTTI/bcdata/cmaq_hourly/sites_static/YYYY/MM/cmaq_airnow_pm2.5_o3_static_YYYYMMDD_HH00_use.nc
cams_o3_pm25_regrid_cmaq_YYYYMMDD_HH00.nc
/glade/campaign/ral/nsap/JTTI/bcdata/cmaq_hourly/cmaq_airnow_pm2.5_o3_YYYYMMDDHH.nc
Output:
/glade/campaign/ral/nsap/JTTI/merge/BM3/static/cams_include0/cams_pm25_regrid_cmaq_YYYYMMDD_HH00_BM3_static_noQC.nc
/glade/campaign/ral/nsap/JTTI/merge/BM3/static/cmaq_include0_nn/cmaq_pm25_YYYYMMDD_HH00_BM3_static_noQC.nc

# (Scott) Evaluating BM models and analysis (CMAQ raw and CAMS raw)
Input:
CAMS/AirNOW: /glade/campaign/ral/nsap/JTTI/cams/fcst/YYYY/MM/YYYYMMDD_HH/cams_airnow_pm2.5_o3_static_YYYYMMDD_HH00_val.nc
CMAQ: /glade/campaign/ral/nsap/JTTI/bcdata/grid/*
BM CAMS: /glade/campaign/ral/nsap/JTTI/merge/BM3/static/cams_include0/cams_pm25_regrid_cmaq_YYYYMMDD_HH00_BM3_static_noQC.nc
BM CMAQ: /glade/campaign/ral/nsap/JTTI/merge/BM3/static/cmaq_include0_nn/cmaq_pm25_YYYYMMDD_HH00_BM3_static_noQC.nc

# (Jared) Preparing CMAQ forecast, AirNOW, BM_CAMS, BM_CMAQ collocated data for 0–72 hours lead time at the training and verification stations
python extract_cmaq_raw_bm_cams_bm_at_airnow_use_val_leadtimeloop.py 20200801_06 -l 20211228_06 -f 72
Inputs:
AirNow obs, training stations: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_use.nc
AirNow obs, validation stations: /glade/campaign/ral/nsap/JTTI/bcdata/airnow/split/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_val.nc
CAMS BM3 PM2.5: /glade/campaign/ral/nsap/JTTI/merge/BM3/static/cams_include0/cams_pm25_regrid_cmaq_YYYYMMDD_HH00_BM3_static_noQC.nc
CMAQ BM3 PM2.5:
/glade/campaign/ral/nsap/JTTI/merge/BM3/static/cmaq_include0_nn/cmaq_pm25_YYYYMMDD_HH00_BM3_static_noQC.nc
CMAQ BM3 O3: /glade/campaign/ral/nsap/JTTI/merge/BM3/static/cams_include0/cams_o3_regrid_cmaq_YYYYMMDD_HH00_BM3_static_include0.nc
CAMS BM3 O3: /glade/campaign/ral/nsap/JTTI/merge/BM3/static/cmaq_include0_nn/cmaq_o3_YYYYMMDD_HH00_BM3_static_inc0_nn.nc
CMAQ raw forecasts with 9 variables for AnEn predictors: /glade/campaign/ral/nsap/JTTI/06z_cmaq_predictors/cmaq_9variables_YYYYMMDD.06z.nc
Outputs:
/glade/campaign/ral/nsap/JTTI/cmaq/raw/airnow/sites_static/YYYY/MM/cmaq_raw_cmaq_bm_cams_bm_airnow_pm2.5_o3_YYYYMMDD_0600_use.nc
/glade/campaign/ral/nsap/JTTI/cmaq/raw/airnow/sites_static/YYYY/MM/cmaq_raw_cmaq_bm_cams_bm_airnow_pm2.5_o3_YYYYMMDD_0600_val.nc

# (Ju-Hye, task 10 in G space) Prepare input  to the AnEn CMAQ file already available with the predictors. We need one additional file with the obs which are BM CAMS at the training stations and another additional file for the obs which is CMAQ BM3 for both O3 and PM25.
Input:
s
Output:
/glade/campaign/ral/nsap/JTTI/anen_input/InputAnEn.O3_rawOBS_fcstCMAQ.517st_use.72hr.20200801-20211231.All.06z.add_camsBM3_cmaqBM3.nc
/glade/campaign/ral/nsap/JTTI/anen_input/InputAnEn.PM25_rawOBS_fcstCMAQ.574st_use.72hr.20200801-20211231.All.06z.add_camsBM3_cmaqBM3.nc
(Ju-Hye, task 11 in Gspace) input  to the AnEn CMAQ file with the raw CMAQ predictors at the verification stations. we need one additional file with the obs which are BM CAMS at the verification stations and another additional file for the obs which is CMAQ BM3 for both O3 and PM25
Input:


Output:
/glade/campaign/ral/nsap/JTTI/anen_input/InputAnEn.O3_rawOBS_fcstCMAQ.675st_val.72hr.20200801-20211231.All.06z.add_camsBM3_cmaqBM3.nc
/glade/campaign/ral/nsap/JTTI/anen_input/InputAnEn.PM25_rawOBS_fcstCMAQ.387st_val.72hr.20200801-20211231.All.06z.add_camsBM3_cmaqBM3.nc
(Stefano) Running the AnEn over the training stations for PM2.5 and Ozone. 
Input
/glade/campaign/ral/nsap/JTTI/anen_input/InputAnEn.O3_rawOBS_fcstCMAQ.675st_val.72hr.20200801-20211231.All.06z.add_camsBM3_cmaqBM3.nc
/glade/campaign/ral/nsap/JTTI/anen_input/InputAnEn.PM25_rawOBS_fcstCMAQ.387st_val.72hr.20200801-20211231.All.06z.add_camsBM3_cmaqBM3.nc
Output

/glade/u/home/alessand/anen_output_o3_fromobs_use.nc 
/glade/u/home/alessand/anen_output_pm25_fromobs_use.nc

# (Stefano) Compute the difference between the mean of the AnEn forecast and CMAQ raw forecast (which is input to the AnEn)
Input
/glade/u/home/alessand/anen_output_o3_fromobs_use.nc 
/glade/u/home/alessand/anen_output_pm25_fromobs_use.nc
Output
/glade/u/home/alessand/anen_output_mean-cmaq_PM25_conus_use.nc
/glade/u/home/alessand/anen_output_mean_PM25_conus_use.nc
/glade/u/home/alessand/anen_output_mean_O3_conus_use.nc
/glade/u/home/alessand/anen_output_mean-cmaq_O3_conus_use.nc

 (Ju-Hye) Spread the difference from step 9 over CONUS
Input : AnEn outputs from 9) and 10)
Output : /glade/scratch/jkim/CMAQ/bm3_cmaq_anen_pm25_final/cmaq_pm25_AnEn_bm3_${date}_06z_LT${time}00.nc

# (Ju-Hye) For each forecast lead time and run, compute the field resulting from adding step 10) output to the original CMAQ field (PM25 and Ozone).
Input: Outputs of 11)
Output
/glade/campaign/ral/nsap/JTTI/bc_cmaq_06z_using_anen/bc_cmaq_06z_pm2.5_o3_06z_YYYYMMDD.nc
(Jared) Extract the values of Bias corrected CMAQ from step 12 at the verification station over CONUS (this script also automatically does the same over just the California subdomain). This is the currently operational NOAA approach.
python extract_anen_cmaq_bm_at_airnow_use_val_leadtime.py 20210101_06 -l 20211228_06
Inputs:
/glade/campaign/ral/nsap/JTTI/bcdata/airnow/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_use.nc
/glade/campaign/ral/nsap/JTTI/bcdata/airnow/sites_static/20200801_00-20211231_23/YYYY/MM/airnow_pm2.5_o3_YYYYMMDD_HH00_val.nc
/glade/campaign/ral/nsap/JTTI/bc_cmaq_06z_using_anen/bc_cmaq_06z_pm2.5_o3_06z_YYYYMMDD.nc
Outputs:
/glade/campaign/ral/nsap/JTTI/cmaq/anen_bc/airnow/sites_static/.YYYY/MM/conus_anen_bc_cmaq_airnow_pm2.5_o3_YYYYMMDD_HH00_use.nc
/glade/campaign/ral/nsap/JTTI/cmaq/anen_bc/airnow/sites_static/.YYYY/MM/conus_anen_bc_cmaq_airnow_pm2.5_o3_YYYYMMDD_HH00_val.nc
/glade/campaign/ral/nsap/JTTI/cmaq/anen_bc/airnow/sites_static/.YYYY/MM/calif_anen_bc_cmaq_airnow_pm2.5_o3_YYYYMMDD_HH00_use.nc
/glade/campaign/ral/nsap/JTTI/cmaq/anen_bc/airnow/sites_static/.YYYY/MM/calif_anen_bc_cmaq_airnow_pm2.5_o3_YYYYMMDD_HH00_val.nc
(Stefano) Running the AnEn over the verification stations for PM2.5 and Ozone over CONUS (both using BM CAMS and BM CMAQ analysis)
/glade/p/ral/nsap/alessand/CMAQ_JTTI/New_obs_Conus/Val/O3/Anen_output_o3_fromcamsbm3_val.nc
/glade/p/ral/nsap/alessand/CMAQ_JTTI/New_obs_Conus/Val/O3/Anen_output_o3_fromcmaqbm3_val.nc
/glade/p/ral/nsap/alessand/CMAQ_JTTI/New_obs_Conus/Val/PM25/Anen_output_pm25_fromcmaqbm3_val.nc
/glade/p/ral/nsap/alessand/CMAQ_JTTI/New_obs_Conus/Val/PM25/anen_output_pm25_fromcamsbm3_val.nc

# (Stefano) Final verification against obs over CONUS using AnEn over verification stations with BM  CAMS and BM CMAQ analysis and output from step 12 (which is the current NOAA operational approach)

