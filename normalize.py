import sys
import numpy as np
import scipy as sp
import pandas as pd
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print("Urban area or pollutants specified improperly--using defaults...")
    # set defaults
    urban_area = "AtlantaGA"
    pollutant = "O3"
else:
    urban_area = sys.argv[1]
    pollutant = sys.argv[2]

data_ncdf = f"./NetCDFs/{urban_area}{pollutant}.nc"
pol_dataset = xr.open_dataset(data_ncdf)

norm = (pol_dataset[f'{pollutant} concentration'].to_numpy() / np.nanmax(pol_dataset[f'{pollutant} concentration'].to_numpy())).reshape(tuple(pol_dataset.dims.values()))
time = pol_dataset['time'].values
band = pol_dataset['band'].values
y = pol_dataset['y'].values
x = pol_dataset['x'].values
norm_dr = xr.DataArray(norm, coords={'time': time, 'band': band, 'y': y, 'x': x}, 
                        dims=['time', 'band', 'y', 'x'])

pol_dataset[f"{pollutant}norm concentration"] = norm_dr
norm_dataset = pol_dataset.drop(f"{pollutant} concentration")
norm_dataset.to_netcdf(f"./NetCDFs/{urban_area}{pollutant}norm.nc")