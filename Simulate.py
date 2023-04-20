from PumpStation import *
from OptimizationProblem import PumpStationOptimization
import matplotlib.pyplot as plt
import numpy as np

f = PumpStationOptimization()
f.simulate_maximum_demand = True
f.simulate_solar = True
f.simulate_solar_exchange = True
f.simulate_battery = True

f.build_problem()
f.build_constraints()
f.solve()

print("Objective function", f.res.fun)
print("Energy purchased", np.dot(f.solution['energy_purchased'], buying_prices))
print("Maximum demand", f.solution['maximum_demand']*md_price_control_horizon)
print("Energy solar used", sum(f.solution['energy_solar']))
print("Energy battery discharging", sum(f.solution['energy_battery_discharging'])*battery_power_step_size)

# # plot results

fig, axes = plt.subplots(nrows=7, ncols=1, figsize=(5,8))

for ax in axes:
    ax.grid(True)
    ax.set_xlim(0,periods)
    ax.set_xticks(np.arange(0, periods+1, 1))

    if ax != axes[-1]:
        ax.get_xaxis().set_ticklabels([])    

def step_plot(ax, data, style='-'):
    ax.step(range(0,len(data)+1), [*data, data[-1]], where='post', linestyle=style)

step_plot(axes[0], tariff_RTP_SARIMA)
step_plot(axes[0], tariff_RTP)
axes[0].set_ylabel('RTP USD/MWh')

axes[1].set_yticks([reservoir_lb, reservoir_initial_level, reservoir_ub])
axes[1].plot([0, periods], [reservoir_lb, reservoir_lb], '--r', alpha=0.5)
axes[1].plot([0, periods], [reservoir_initial_level, reservoir_initial_level], '--k', alpha=0.5)
axes[1].plot([0, periods], [reservoir_ub, reservoir_ub], '--r', alpha=0.5)
axes[1].plot(f.solution['reservoir_level'])
axes[1].set_ylabel('Reservoir level m3')

step_plot(axes[2], f.solution['pump_status'])
axes[2].set_yticks([1,3,5,7])
axes[2].set_ylim(0,7+0.5)
axes[2].set_ylabel('Active pumps')

step_plot(axes[3], np.array(f.solution['energy_purchased']))
axes[3].set_ylabel('$E_p$ MWh')

step_plot(axes[4], f.solution['energy_solar'])
axes[4].set_ylabel('$E_s$ MWh')

axes[5].set_ylabel('$E_d$ MWh')
step_plot(axes[5], (f.solution['energy_battery_discharging'] - f.solution['energy_battery_charging'])*battery_power_step_size)

axes[6].set_ylabel('SoC MWh')
axes[6].plot(f.solution['battery_level'])

fig.tight_layout()
plt.show()
