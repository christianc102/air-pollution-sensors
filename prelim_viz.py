import sys
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt
import matplotlib.animation as animation

if len(sys.argv) != 3:
    print("Urban area or pollutant specified improperly--using defaults...")
    # set defaults
    urban_area = "AtlantaGA"
    pollutant = "NO2"
else:
    urban_area = sys.argv[1]
    pollutant = sys.argv[2]

input_ncdf = f"./NetCDFs/{urban_area}{pollutant}.nc"
dataset = xr.open_dataset(input_ncdf)
# .sel(time=slice('2000-01-01', '2000-02-01'))

mean_data_t_series = dataset.mean(dim=["x", "y"])
fig = plt.figure()
mean_data_t_series_df = mean_data_t_series.to_dataframe()[[f"{pollutant} concentration"]].droplevel(1)
sns.lineplot(data=mean_data_t_series_df)
plt.xticks(rotation=90)
plt.savefig(f"./Visualizations/{urban_area}{pollutant}mean_time_series.png")

bounds_dict = {"NO2": 53, "O3": 70, "PM": 35}

pollution_data = dataset[f"{pollutant} concentration"].squeeze('band')
fig = plt.figure()

def animate(frame):
    fig.clear()
    ax = plt.subplot()
    current_data = pollution_data[frame, :, :]
    img = ax.imshow(current_data, origin='lower', vmin=0, vmax=bounds_dict[pollutant])
    ax.set_title(f'{pollutant} Frame {frame}')
    plt.colorbar(img, ax=ax)

num_frames = pollution_data.shape[0]
anim = animation.FuncAnimation(fig, animate, frames=num_frames, interval=250)
writer = animation.FFMpegWriter(fps=8)
anim.save(f"./Visualizations/{urban_area}{pollutant}.mp4", writer=writer)