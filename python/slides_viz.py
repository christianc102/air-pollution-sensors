import sys
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from statsmodels.tsa.stattools import adfuller
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.tsa.arima.model import ARIMA
from sklearn.preprocessing import StandardScaler

sns.set_theme()
sns.set(font_scale=1.5)
sns.color_palette("deep")

units_dict = {"NO2": " (ppb)", "O3": " (ppb)", "PM": u" (\u03bcg m$^{-3}$)", "O3NO2ratio": "", "O3PMratio": " (ppb/u\u03bcg m$^{-3}$)"}

if len(sys.argv) != 2:
    print("Urban area or pollutant specified improperly--using defaults...")
    # set defaults
    urban_area = "AtlantaGA"
else:
    urban_area = sys.argv[1]

pollutants = ["NO2", "PM", "O3"]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8,10), sharex=True)
axes = [ax1, ax2, ax3]

for i in range(len(pollutants)):
    input_ncdf = f"./NetCDFs/{urban_area}{pollutants[i]}.nc"
    dataset = xr.open_dataset(input_ncdf)
    # .sel(time=slice('2000-01-01', '2000-02-01'))
    times = dataset['time'].dt.strftime("%m/%d/%Y")

    mean_data_t_series = dataset.mean(dim=["x", "y"])
    mean_data_t_series_df = mean_data_t_series.to_dataframe()[[f"{pollutants[i]} concentration"]].droplevel(1)
    mean_data_t_series_df = mean_data_t_series_df.rename(columns={f"{pollutants[i]} concentration": f"{pollutants[i]}{units_dict[pollutants[i]]}"})
    window_length = [50, 200, 730]
    mean_data_t_series_df[f'{window_length[0]}-day moving average'] = mean_data_t_series_df[f"{pollutants[i]}{units_dict[pollutants[i]]}"].rolling(window=window_length[0]).mean()
    mean_data_t_series_df[f'{window_length[1]}-day moving average'] = mean_data_t_series_df[f"{pollutants[i]}{units_dict[pollutants[i]]}"].rolling(window=window_length[1]).mean()
    mean_data_t_series_df[f'{window_length[2]//365}-yr moving average'] = mean_data_t_series_df[f"{pollutants[i]}{units_dict[pollutants[i]]}"].rolling(window=window_length[2]).mean()
    sns.lineplot(x=mean_data_t_series_df.index, y=f"{pollutants[i]}{units_dict[pollutants[i]]}", ax=axes[i], data=mean_data_t_series_df, alpha=0.5, linewidth=0.5, label="Data", legend=False)
    sns.lineplot(x=mean_data_t_series_df.index, y=f"{window_length[0]}-day moving average", ax=axes[i], data=mean_data_t_series_df, label=f"{window_length[0]}-day moving average", legend=False)
    sns.lineplot(x=mean_data_t_series_df.index, y=f"{window_length[1]}-day moving average", ax=axes[i], data=mean_data_t_series_df, label=f"{window_length[1]}-day moving average", legend=False)
    sns.lineplot(x=mean_data_t_series_df.index, y=f"{window_length[2]//365}-yr moving average", ax=axes[i], data=mean_data_t_series_df, label=f"{window_length[2]//365}-year moving average", legend=False)
    plt.xticks(fontsize=12)
    plt.tight_layout()
    plt.savefig(f"./Visualizations/{urban_area}mean_time_series.png")