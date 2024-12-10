"""
plot_map_epa_regions.py

Created by: Jared A. Lee (jaredlee@ucar.edu)
Created on: 16 Aug 2024

This script generates a map of CONUS, with states color-coded by what EPA Region they're in.

Derived from: https://scitools.org.uk/cartopy/docs/v0.15/examples/hurricane_katrina.html
"""
from typing import Dict, Any
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

plot_type = 'eps'
plot_dir = pathlib.Path('/', 'glade', 'campaign', 'ral', 'nsap', 'JTTI', 'plots')
plot_file = 'map_epa_regions.' + plot_type
fname = plot_dir.joinpath(plot_file)

mpl.rcParams['figure.figsize'] = (6, 4)
mpl.rcParams['savefig.bbox'] = 'tight'

data_crs = ccrs.LambertConformal()
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection=data_crs)
ax.set_extent([-120, -72, 22, 50], ccrs.Geodetic())

shapename = 'admin_1_states_provinces_lakes'
states_shp = shpreader.natural_earth(resolution='10m', category='cultural', name=shapename)

epa_regions = {
    'Maine': 1,
    'New Hampshire': 1,
    'Vermont': 1,
    'Massachusetts': 1,
    'Rhode Island': 1,
    'Connecticut': 1,
    'New York': 2,
    'New Jersey': 2,
    'Pennsylvania': 3,
    'Maryland': 3,
    'Delaware': 3,
    'West Virginia': 3,
    'Virginia': 3,
    'Kentucky': 4,
    'Tennessee': 4,
    'North Carolina': 4,
    'South Carolina': 4,
    'Mississippi': 4,
    'Alabama': 4,
    'Georgia': 4,
    'Florida': 4,
    'Minnesota': 5,
    'Wisconsin': 5,
    'Michigan': 5,
    'Illinois': 5,
    'Indiana': 5,
    'Ohio': 5,
    'New Mexico': 6,
    'Texas': 6,
    'Oklahoma': 6,
    'Arkansas': 6,
    'Louisiana': 6,
    'Nebraska': 7,
    'Iowa': 7,
    'Kansas': 7,
    'Missouri': 7,
    'North Dakota': 8,
    'South Dakota': 8,
    'Montana': 8,
    'Wyoming': 8,
    'Utah': 8,
    'Colorado': 8,
    'California': 9,
    'Nevada': 9,
    'Arizona': 9,
    'Washington': 10,
    'Oregon': 10,
    'Idaho': 10,
    }

ax.set_title('EPA Regions', fontsize=16)

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
        'tab:brown', 'tab:olive', 'tab:gray', 'tab:pink', 'tab:cyan', 'black']

for astate in shpreader.Reader(states_shp).records():

    edgecolor = 'black'

    try:
        # use the name of this state to get the EPA region
        state_reg = epa_regions[astate.attributes['name']]
    except:
        state_reg = 0

    # simple scheme to assign color to each state
    if state_reg == 1:
        facecolor = colors[0]
    elif state_reg == 2:
        facecolor = colors[1]
    elif state_reg == 3:
        facecolor = colors[2]
    elif state_reg == 4:
        facecolor = colors[3]
    elif state_reg == 5:
        facecolor = colors[4]
    elif state_reg == 6:
        facecolor = colors[5]
    elif state_reg == 7:
        facecolor = colors[6]
    elif state_reg == 8:
        facecolor = colors[7]
    elif state_reg == 9:
        facecolor = colors[8]
    elif state_reg == 10:
        facecolor = colors[9]
    else:
        facecolor = 'white'

    # `astate.geometry` is the polygon to plot
    ax.add_geometries([astate.geometry], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)

    ax.text(-70.0, 45.0, '1', size=16, weight='heavy', transform=ccrs.PlateCarree())
    ax.text(-76.0, 42.5, '2', size=16, weight='heavy', transform=ccrs.PlateCarree())
    ax.text(-79.25, 37.25, '3', size=16, weight='heavy', transform=ccrs.PlateCarree())
    ax.text(-84.75, 33.0, '4', size=16, weight='heavy', transform=ccrs.PlateCarree())
    ax.text(-90.5, 43.5, '5', size=16, weight='heavy', transform=ccrs.PlateCarree())
    ax.text(-99.5, 31.0, '6', size=16, weight='heavy', transform=ccrs.PlateCarree())
    ax.text(-99.0, 40.5, '7', size=16, weight='heavy', transform=ccrs.PlateCarree())
    ax.text(-108.0, 42.0, '8', size=16, weight='heavy', transform=ccrs.PlateCarree())
    ax.text(-117.0, 38.0, '9', size=16, weight='heavy', transform=ccrs.PlateCarree())
    ax.text(-121.0, 43.5, '10', size=16, weight='heavy', transform=ccrs.PlateCarree())


# plt.show()
print('Plotting '+str(fname))
plt.savefig(fname)
plt.close()
fname.chmod(0o666)