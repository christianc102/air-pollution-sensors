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
sns.color_palette("deep")

units_dict = {"NO2": " (ppb)", "O3": " (ppb)", "PM": u" (\u03bcg m$^{-3}$)", "O3NO2ratio": "", "O3PMratio": u" (ppb/\u03bcg m$^{-3}$)", "PMNO2ratio": u" (\u03bcg m$^{-3}$/ppb)", "O3normNO2normratio": ""}

if len(sys.argv) != 3:
    print("Urban area or pollutant specified improperly--using defaults...")
    # set defaults
    urban_area = "AtlantaGA"
    pollutant = "NO2"
else:
    urban_area = sys.argv[1]
    pollutant = sys.argv[2]

input_ncdf = f"./NetCDFs/{urban_area}{pollutant}.nc"
dataset = xr.open_dataset(input_ncdf).sel(time=slice('2008-11-01', '2008-11-30'))
times = dataset['time'].dt.strftime("%m/%d/%Y")
print(dataset[f'{pollutant} concentration'])

# mean_data_t_series = dataset.mean(dim=["x", "y"])
# fig = plt.figure()
# mean_data_t_series_df = mean_data_t_series.to_dataframe()[[f"{pollutant} concentration"]].droplevel(1)
# mean_data_t_series_df = mean_data_t_series_df.rename(columns={f"{pollutant} concentration": f"{pollutant}{units_dict[pollutant]}"})
# print(mean_data_t_series_df.info())
# # scaler = StandardScaler()
# # mean_concs = mean_data_t_series_df[f"{pollutant}{units_dict[pollutant]}"].replace(np.inf, 0)
# # standardized_mean_t_series = scaler.fit_transform(mean_concs.to_numpy().reshape(-1,1))
# # standardized_mean_t_series_df = pd.DataFrame(standardized_mean_t_series, index=mean_data_t_series_df.index, columns=[f"{pollutant}{units_dict[pollutant]}"])
# window_length = [50, 200, 730]
# mean_data_t_series_df[f'{window_length[0]}-day moving average'] = mean_data_t_series_df[f"{pollutant}{units_dict[pollutant]}"].rolling(window=window_length[0]).mean()
# mean_data_t_series_df[f'{window_length[1]}-day moving average'] = mean_data_t_series_df[f"{pollutant}{units_dict[pollutant]}"].rolling(window=window_length[1]).mean()
# mean_data_t_series_df[f'{window_length[2]//365}-yr moving average'] = mean_data_t_series_df[f"{pollutant}{units_dict[pollutant]}"].rolling(window=window_length[2]).mean()
# sns.lineplot(x=mean_data_t_series_df.index, y=f"{pollutant}{units_dict[pollutant]}", data=mean_data_t_series_df, alpha=0.5, linewidth=0.5, label="Data", legend=False)
# sns.lineplot(x=mean_data_t_series_df.index, y=f"{window_length[0]}-day moving average", data=mean_data_t_series_df, label=f"{window_length[0]}-day moving average", legend=False)
# sns.lineplot(x=mean_data_t_series_df.index, y=f"{window_length[1]}-day moving average", data=mean_data_t_series_df, label=f"{window_length[1]}-day moving average", legend=False)
# sns.lineplot(x=mean_data_t_series_df.index, y=f"{window_length[2]//365}-yr moving average", data=mean_data_t_series_df, label=f"{window_length[2]//365}-year moving average", legend=False)
# plt.xticks(fontsize=12)
# plt.savefig(f"./Visualizations/{urban_area}{pollutant}mean_time_series.png")

# fig = plt.figure()
# standardized_mean_t_series_df[f'{window_length[0]}-day moving average'] = standardized_mean_t_series_df[f"{pollutant}{units_dict[pollutant]}"].rolling(window=window_length[0]).mean()
# standardized_mean_t_series_df[f'{window_length[1]}-day moving average'] = standardized_mean_t_series_df[f"{pollutant}{units_dict[pollutant]}"].rolling(window=window_length[1]).mean()
# standardized_mean_t_series_df[f'{window_length[2]//365}-yr moving average'] = standardized_mean_t_series_df[f"{pollutant}{units_dict[pollutant]}"].rolling(window=window_length[2]).mean()
# print(standardized_mean_t_series_df.info())
# sns.lineplot(x=standardized_mean_t_series_df.index, y=f"{pollutant}{units_dict[pollutant]}", data=standardized_mean_t_series_df, alpha=0.5, linewidth=0.5, label="Data")
# sns.lineplot(x=standardized_mean_t_series_df.index, y=f"{window_length[0]}-day moving average", data=standardized_mean_t_series_df, label=f"{window_length[0]}-day moving average")
# sns.lineplot(x=standardized_mean_t_series_df.index, y=f"{window_length[1]}-day moving average", data=standardized_mean_t_series_df, label=f"{window_length[1]}-day moving average")
# sns.lineplot(x=standardized_mean_t_series_df.index, y=f"{window_length[2]//365}-yr moving average", data=standardized_mean_t_series_df, label=f"{window_length[2]//365}-year moving average")
# plt.xticks(fontsize=12)
# plt.savefig(f"./Visualizations/{urban_area}{pollutant}std_mean_time_series.png")

# # find order of d in ARIMA
# result = adfuller(mean_data_t_series_df[f"{pollutant} concentration"].to_numpy())
# print(f'ADF Statistic: {result[0]}')
# print(f'p-value: {result[1]}')

# # find order of p in ARIMA
# fig, axes = plt.subplots(1, 2, sharex=True)
# axes[0].plot(np.ediff1d(mean_data_t_series_df[f"{pollutant} concentration"].to_numpy()))
# axes[0].set_title('1st Differencing')
# axes[1].set(ylim=(0,5), xlim=(0,40))
# plot_pacf(np.ediff1d(mean_data_t_series_df[f"{pollutant} concentration"].to_numpy()), ax=axes[1])
# plt.savefig("./Tests/PACF.png")

# # find order of q in ARIMA
# fig, axes = plt.subplots(1, 2, sharex=True)
# axes[0].plot(np.ediff1d(mean_data_t_series_df[f"{pollutant} concentration"].to_numpy()))
# axes[0].set_title('1st Differencing')
# axes[1].set(ylim=(0,5), xlim=(0,40))
# plot_acf(np.ediff1d(mean_data_t_series_df[f"{pollutant} concentration"].to_numpy()), ax=axes[1])
# plt.savefig("./Tests/ACF.png")

# # build ARIMA model and plot
# model = ARIMA(mean_data_t_series_df[f"{pollutant} concentration"], order=(0,0,1))
# model_fit = model.fit()
# print(model_fit.summary())

# residuals = pd.DataFrame(model_fit.resid)
# fig, ax = plt.subplots(1,2)
# residuals.plot(title="Residuals", ax=ax[0])
# residuals.plot(kind='kde', title='Density', ax=ax[1])
# plt.savefig("./Tests/Residuals.png")

# mean_data_t_series_df["ARIMA"] = model_fit.predict()
# fig = plt.figure()
# sns.lineplot(data=mean_data_t_series_df)
# plt.xticks(rotation=90)
# plt.savefig(f"./Visualizations/{urban_area}{pollutant}ARIMAmean_time_series.png")

bounds_dict = {"NO2": 53, "O3": 70, "PM": 35, "O3NO2ratio": 4, "O3PMratio": 8, "PMNO2ratio": 2, "O3normNO2normratio": 4}

pollution_data = dataset[f"{pollutant} concentration"].squeeze('band')
fig = plt.figure()

def animate(frame):
    fig.clear()
    plt.set_cmap('inferno_r')
    ax = plt.subplot()
    current_data = pollution_data[frame, :, :]
    img = ax.imshow(current_data, vmin=0, vmax=bounds_dict[pollutant])
    current_time = times[frame].values
    ax.set_title(f'{pollutant}{units_dict[pollutant]} - {current_time}')
    plt.colorbar(img, ax=ax)

num_frames = pollution_data.shape[0]
anim = animation.FuncAnimation(fig, animate, frames=num_frames, interval=250)
writer = animation.FFMpegWriter(fps=8)
anim.save(f"./Visualizations/{urban_area}{pollutant}.mp4", writer=writer)