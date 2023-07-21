import sys
import numpy as np
import scipy as sp
import pandas as pd
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt

if len(sys.argv) != 4:
    print("Urban area or pollutants specified improperly--using defaults...")
    # set defaults
    urban_area = "AtlantaGA"
    pollutant1 = "O3"
    pollutant2 = "NO2"
else:
    urban_area = sys.argv[1]
    pollutant1 = sys.argv[2]
    pollutant2 = sys.argv[3]

num_ncdf = f"./NetCDFs/{urban_area}{pollutant1}.nc"
num_dataset = xr.open_dataset(num_ncdf)

denom_ncdf = f"./NetCDFs/{urban_area}{pollutant2}.nc"
denom_dataset = xr.open_dataset(denom_ncdf)

ratios = (num_dataset[f'{pollutant1} concentration'].to_numpy() / denom_dataset[f'{pollutant2} concentration'].to_numpy()).reshape(tuple(num_dataset.dims.values()))
time = num_dataset['time'].values
band = num_dataset['band'].values
y = num_dataset['y'].values
x = num_dataset['x'].values
ratio_dr = xr.DataArray(ratios, coords={'time': time, 'band': band, 'y': y, 'x': x}, 
                        dims=['time', 'band', 'y', 'x'])

num_dataset[f"{pollutant1}{pollutant2}ratio concentration"] = ratio_dr
ratio_dataset = num_dataset.drop(f"{pollutant1} concentration")
ratio_dataset.to_netcdf(f"./NetCDFs/{urban_area}{pollutant1}{pollutant2}ratio.nc")