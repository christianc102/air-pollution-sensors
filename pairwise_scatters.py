import sys
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt

sns.set_theme()
sns.color_palette("deep")

units_dict = {"NO2": " (ppb)", "O3": " (ppb)", "PM": u" (\u03bcg m$^{-3}$)"}

if len(sys.argv) != 2:
    print("Urban area or pollutant specified improperly--using defaults...")
    # set defaults
    urban_area = "AtlantaGA"
else:
    urban_area = sys.argv[1]

def load_mean_t_series_df(pollutant):
    input_ncdf = f"./NetCDFs/{urban_area}{pollutant}.nc"
    dataset = xr.open_dataset(input_ncdf)
    times = dataset['time'].dt.strftime("%m/%d/%Y")
    mean_data_t_series = dataset.mean(dim=["x", "y"])
    mean_data_t_series_df = mean_data_t_series.to_dataframe()[[f"{pollutant} concentration"]].droplevel(1)
    mean_data_t_series_df = mean_data_t_series_df.rename(columns={f"{pollutant} concentration": f"{pollutant}{units_dict[pollutant]}"})
    return mean_data_t_series_df

def get_season(month):
    if 3 <= month <= 5:
        return "MAM"
    elif 6 <= month <= 8:
        return "JJA"
    elif 9 <= month <= 11:
        return "SON"
    else:
        return "DJF"
    

no2_df = load_mean_t_series_df("NO2")
o3_df = load_mean_t_series_df("O3")
pm_df = load_mean_t_series_df("PM")
full_df = pd.concat([no2_df, o3_df, pm_df], axis=1)

full_df['Season'] = full_df.index.month.map(get_season)

fig = plt.figure()
sns.scatterplot(full_df, x=f"PM{units_dict['PM']}", y=f"O3{units_dict['O3']}", hue="Season")
plt.xticks(fontsize=12)
plt.xlim(0, 112)
plt.title(f"Seasonal O3 vs. PM - {urban_area}")
plt.savefig(f"./Visualizations/{urban_area}pairwise_pm_o3.png")

fig = plt.figure()
sns.scatterplot(full_df, x=f"NO2{units_dict['NO2']}", y=f"O3{units_dict['O3']}", hue="Season")
plt.xticks(fontsize=12)
plt.xlim(0, 112)
plt.title(f"Seasonal O3 vs. NO2 - {urban_area}")
plt.savefig(f"./Visualizations/{urban_area}pairwise_no2_o3.png")

fig = plt.figure()
sns.scatterplot(full_df, x=f"NO2{units_dict['NO2']}", y=f"PM{units_dict['PM']}", hue="Season")
plt.xticks(fontsize=12)
plt.title(f"Seasonal PM vs. NO2 - {urban_area}")
plt.savefig(f"./Visualizations/{urban_area}pairwise_no2_pm.png")