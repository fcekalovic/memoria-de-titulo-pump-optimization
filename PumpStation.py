import numpy as np

# Time periods and hours_delta time
days = 2
hours_delta = 1
periods_per_hour = round(hours_delta ** -1)
periods = 24 * days * periods_per_hour

# Number of pumps, flow rate and power rating
number_of_pumps = 7
pump_power = 1.343    # MW
pump_flow_rate = 700  # m3/hour

# Reservoir size, reservoir bounds, reservoir initial level
reservoir_size = 60000                          # m3
reservoir_lb = 0.25 * reservoir_size            # m3
reservoir_ub = 0.75 * reservoir_size            # m3
reservoir_initial_level = 0.5 * reservoir_size
reservoir_maintain_level = True

# Water demand
water_demand_horizon = np.array([*np.repeat(2800, 15), *np.repeat(2100, 9)])
water_demand_horizon = np.tile(water_demand_horizon, days)                       # Repeat the water demand profile times "days"
water_demand_horizon = water_demand_horizon.repeat(periods_per_hour, axis=0)      # Oversample by periods per hour
water_demand_horizon = water_demand_horizon * hours_delta                        # Downscale by hours per period
water_demand_cumulative = np.cumsum(water_demand_horizon)

# Solar energy
pv_scale = 5
hourly_solar_profile = np.array([0.      , 0.      , 0.      , 0.      , 0.      , 0.010513,
       0.171737, 0.421572, 0.617406, 0.755007, 0.835454, 0.858055,
       0.828236, 0.744472, 0.600042, 0.407912, 0.16919 , 0.010574,
       0.      , 0.      , 0.      , 0.      , 0.      , 0.      ])

hourly_solar_profile = np.tile(hourly_solar_profile, days)                       # Repeat the solar profile times "days"
hourly_solar_profile = hourly_solar_profile.repeat(periods_per_hour, axis=0)      # Oversample by periods per hour
hourly_solar_profile = hourly_solar_profile * hours_delta                        # Downscale by hours per period
hourly_solar_profile = hourly_solar_profile * pv_scale

# Battery Bank
# https://www.tesla.com/megapack

battery_number = 4
# 2-hour standard megapack
battery_roundtrip_eff = 0.87
battery_size = 2.529
battery_max_power = 1.264
battery_power_step_size = 0.0843 

# battery operational conditions
battery_min_level = 0
battery_max_level = 1
battery_initial_level = 0.5
battery_mantain_level = True

# RTP tariff
tariff_RTP_SARIMA = [125.65, 118.23, 124.06, 117.17, 122.01, 121.09, 122.46, 123.36, 80.09, 48.16, 38.02, 37.58, 30.85, 19.37, 18.44, 18.43, 18.72, 20.42, 58.74, 118.74, 119.27, 118.86, 118.26, 118.03, 118.8, 112.65, 117.97, 112.13, 116.47, 115.75, 117.01, 117.84, 80.13, 52.31, 43.48, 43.11, 37.24, 27.24, 26.43, 26.43, 26.68, 28.16, 61.59, 113.92, 114.38, 114.02, 113.5, 113.3]
tariff_RTP = [130.01,122.52,112.71,100.91,117.45,121.37,138.58,148.88,66.67,56.05,52.87,59.51,52.16,52.16,33.90,0.00,0.00,0.00,56.70,138.83,138.83,138.83,138.83,131.12,126.49,116.92,110.09,110.09,110.09,110.09,110.09,108.25,69.21,17.04,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,96.07,141.00,141.00,141.00,141.00,124.40]

tariff_RTP_SARIMA = np.array(tariff_RTP_SARIMA[0:days*24]).repeat(periods_per_hour, axis=0)
tariff_RTP = np.array(tariff_RTP[0:days*24]).repeat(periods_per_hour, axis=0)

prices = tariff_RTP_SARIMA

buying_prices = np.array(prices) * (1 + 0.025)
selling_prices = np.array(prices) * (1 - 0.025) * -1

# MD tariff
md_price = 12.5657           # USD/kW
md_price = md_price * 1000   # USD/MW
md_price_control_horizon = md_price * (days / 30)