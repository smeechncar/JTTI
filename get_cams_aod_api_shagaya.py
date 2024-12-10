from ecmwfapi import ECMWFDataServer
import sys
import datetime as dt

# 2018-09-22_00 read error (tried twice)
# 2018-10-06_12 error (try again later)
# 2018-10-07_12 error (try again later)
# 2018-10-08_12 error (try again later)
beg_date = '2019-12-15'
end_date = '2019-12-31'
cycles = ['00', '12']
#cycles = ['12']

target_dir = '/d1/shagaya_data/cams-aerosol_forecast/raw_archive'
param_str = 'aod550'	# add more parameters (variables) with a slash delimiter in the string

server = ECMWFDataServer()

## Parse dates
beg_yr = int(beg_date.split('-')[0])
beg_mo = int(beg_date.split('-')[1])
beg_dy = int(beg_date.split('-')[2])

end_yr = int(end_date.split('-')[0])
end_mo = int(end_date.split('-')[1])
end_dy = int(end_date.split('-')[2])

beg_dt = dt.datetime(beg_yr, beg_mo, beg_dy, 0, 0, 0)
end_dt = dt.datetime(end_yr, end_mo, end_dy, 0, 0, 0)

this_dt = beg_dt

## Loop through each date
while this_dt <= end_dt:
#	print(this_dt)
	this_yr = this_dt.strftime('%Y')
	this_mo = this_dt.strftime('%m')
	this_dy = this_dt.strftime('%d')

	## Loop through each cycle for each date
	for hh in cycles:
		this_hr = hh
		target_fname = 'cams_'+this_yr+this_mo+this_dy+'_'+this_hr+'.nc'

		## Download the data
		## Get one full forecast cycle per request to manage the number of data requests.
		## The real-time files come as one valid time per file (every 1 h), while the archive files will be every 3-h output per cycle per file.
		## The real-time files also don't have a scale factor or offset for aod550, unlike the archived files.
		## The real-time files also have aod550 as a float, whereas in the archive aod550 is a short.
		## Thus, the downloaded archived files will need to be subsequently post-processed to look like the real-time archive.
		server.retrieve({
			'dataset' : 'cams_nrealtime',
			'stream'  : 'oper',
			'date'    : this_yr+'-'+this_mo+'-'+this_dy,
			'time'    : this_hr,
			'levtype' : 'sfc',
			'param'   : param_str,
			'step'    : '0/3/6/9/12/15/18/21/24/27/30/33/36/39/42/45/48/51/54/57/60/63/66/69/72/75/78/81/84/87/90/93/96/99/102/105/108/111/114/117/120',
			'expver'  : '0001',
			'type'    : 'fc',
			'class'   : 'mc',
			'format'  : 'netcdf',
			'grid'    : '0.4/0.4',
			'target'  : target_dir+'/'+target_fname
		})

	this_dt = this_dt + dt.timedelta(days=1)

'''
server.retrieve({
    'dataset' : 'cams_nrealtime',
    'stream'  : 'oper',
    'date'    : '2019-11-01',
    'time'    : '00',
    'levtype' : 'sfc',
    'param'   : 'aod550',
    'step'    : '0/3/6/9/12/15/18/21/24/27/30/33/36/39/42/45/48/51/54/57/60/63/66/69/72/75/78/81/84/87/90/93/96/99/102/105/108/111/114/117/120',
    'expver'  : '0001',
    'type'    : 'fc',
    'class'   : 'mc',
    'format'  : 'netcdf',
    'grid'    : '0.4/0.4',
    'target'  : '/d1/shagaya_data/cams-aerosol_forecast/raw_archive/cams_20191101_00.nc'
})
'''
