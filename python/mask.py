import sys
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("Urban area specified improperly--using default...")
    # set defaults
    urban_area = "AtlantaGA"
else:
    urban_area = sys.argv[1]

pollutant = "O3"
input_ncdf = f"./NetCDFs/{urban_area}{pollutant}.nc"
dataset = xr.open_dataset(input_ncdf)
map = dataset.sel(time='2000-07-01')
map = map[f"{pollutant} concentration"].squeeze(drop=True)
no_nans = map.fillna(0)
bin_mask = no_nans.where(no_nans == 0, 1)
mask = bin_mask.rename("Urban Area").drop_vars("time")
print(mask.shape)
mask.to_netcdf(f"./Masks/{urban_area}mask.nc")