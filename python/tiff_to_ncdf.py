import sys
import glob
import numpy as np
import scipy as sp
import pandas as pd
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print("Urban area or pollutant specified improperly--using defaults...")
    # set defaults
    urban_area = "AtlantaGA"
    pollutant = "NO2"
else:
    urban_area = sys.argv[1]
    pollutant = sys.argv[2]

def time_index_from_filenames(filenames):
    '''generate DatetimeIndex from filenames'''
    return pd.DatetimeIndex([pd.Timestamp(f[-12:-4]) for f in filenames])

files = glob.glob(f"./UAGeoTiffs/{urban_area}/{urban_area}{pollutant}data/**/*.tif", recursive=True)

ordered_files = sorted(files)

# reading one file
test = rxr.open_rasterio(ordered_files[0], default_name=f"{pollutant} concentration")
print(test)
print(test.shape)
test.plot(vmin=0, vmax=30)
plt.savefig("./Tests/test_plot.png")

time = xr.Variable('time', time_index_from_filenames(ordered_files))
chunks = {'x': test.shape[2], 'y': test.shape[1], 'band': test.shape[0]}
concat_data = xr.concat([rxr.open_rasterio(f, chunks=chunks, default_name=f"{pollutant} concentration") for f in ordered_files], dim=time)
concat_data.to_netcdf(f"./NetCDFs/{urban_area}{pollutant}.nc")