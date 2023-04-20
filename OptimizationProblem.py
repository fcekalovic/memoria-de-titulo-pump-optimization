from PumpStation import *
import copy
import matplotlib.pyplot as plt
from scipy.optimize import milp

class PumpStationOptimization:

    """
    Class for optimizing the operation of a pump station with various simulation options, time-based tariff, maximum-demand tariff, solar power,
    grid energy exchange, battery.     
    
    The `build_problem()` method sets up the problem by initializing null decisions, setting bounds on decision
    variables, and defining costs for optimization.
    The `build_constraints()` method sets up all the operationals constraints calling it's own defined methods.
    The `solve()` method uses mixed-integer linear programming (MILP) to find the optimal solution for pump operation and energy management.
    The calculated solution is stored in the `solution` attribute of the object.
    """

    simulate_solar = False
    simulate_solar_exchange = False
    simulate_battery = False
    simulate_maximum_demand = False

    def solve(self):
        # Call the solver
        res = milp(
            c=self.costs,
            constraints=self.constraints,
            integrality=self.integrality,
            bounds=self.bounds,
        )

        if not res.success:
            print(res.message)
            return
    
        self.res = res
        self.solution = {}

        if self.simulate_maximum_demand:
            maximum_demand = res.x[-1]
            self.solution['maximum_demand'] = maximum_demand
            out_vector = res.x[0:-1]
        else:
            out_vector = res.x

        simulation_modes = {
                (False, False, False): 2,
                (True, False, False): 3,
                (True, True, False): 4,
                (False, False, True): 6,
                (True, False, True): 7,
                (True, True, True): 8,
            }

        num_splits = simulation_modes[(self.simulate_solar, self.simulate_solar_exchange, self.simulate_battery)]
        out_splits = np.array_split(out_vector, num_splits)

        self.solution['pump_status'] = out_splits[0]
        self.solution['energy_purchased'] = out_splits[1]

        if self.simulate_solar:
            self.solution['energy_solar'] = out_splits[2]
        if self.simulate_solar_exchange:
            self.solution['energy_solar_exchange'] = out_splits[3]

        if self.simulate_battery:
            battery_initial_split = 2
            if self.simulate_solar:
                battery_initial_split += 1
            if self.simulate_solar_exchange:
                battery_initial_split += 1
            self.solution['battery_charging'] = out_splits[battery_initial_split]
            self.solution['battery_discharging'] = out_splits[battery_initial_split+1]
            self.solution['energy_battery_charging'] = out_splits[battery_initial_split+2]
            self.solution['energy_battery_discharging'] = out_splits[battery_initial_split+3]   

        self.solution['reservoir_level'] = self.calculate_reservoir_levels()

        if self.simulate_battery:
            self.solution['battery_level'] = self.calculate_battery_levels()

    def calculate_battery_levels(self):
        total_charge_bat_level = self.solution['energy_battery_discharging'] - self.solution['energy_battery_charging']

        battery_status = [battery_initial_level * battery_number * battery_size]
        battery_level = battery_initial_level * battery_number * battery_size

        for battery_state_level in total_charge_bat_level:
            battery_level -= battery_state_level * battery_power_step_size * hours_delta
            battery_status.append(battery_level)

        return battery_status

    def calculate_reservoir_levels(self):
        pump_status = self.solution['pump_status']

        reservoir_status = [reservoir_initial_level]
        reservoir_level = reservoir_initial_level

        for time, pump_number in enumerate(pump_status):
            reservoir_level += pump_number * pump_flow_rate * hours_delta
            reservoir_level -= water_demand_horizon[time]        
            reservoir_status.append(reservoir_level)

        return reservoir_status

    def build_problem(self):
        def build_null_decisions(self):
            self.null_dict = {}

            self.null_dict["pump_status"] = np.zeros(periods)
            self.null_dict["energy_purchased"] = np.zeros(periods)

            if self.simulate_solar:
                self.null_dict["energy_solar"] = np.zeros(periods)
                if self.simulate_solar_exchange:
                    self.null_dict["energy_solar_sold"] = np.zeros(periods)

            if self.simulate_battery:
                self.null_dict["battery_charging"] = np.zeros(periods)
                self.null_dict["battery_discharging"] = np.zeros(periods)
                self.null_dict["energy_battery_charging"] = np.zeros(periods)
                self.null_dict["energy_battery_discharging"] = np.zeros(periods)

            if self.simulate_maximum_demand:
                self.null_dict["maximum_demand"] = [0]

        def build_costs(self):
            self.costs_dict = copy.deepcopy(self.null_dict)

            self.costs_dict["energy_purchased"] = buying_prices

            if self.simulate_solar and self.simulate_solar_exchange:
                self.costs_dict["energy_solar_sold"] = selling_prices

            if self.simulate_maximum_demand:
                self.costs_dict["maximum_demand"] = [md_price_control_horizon]

            self.costs = np.concatenate(list(self.costs_dict.values()))

        def build_bounds(self):
            lower_bounds, upper_bounds = self.bounds_dict = [
                copy.deepcopy(self.null_dict),
                copy.deepcopy(self.null_dict),
            ]

            upper_bounds["pump_status"] = np.ones(periods) * number_of_pumps
            upper_bounds["energy_purchased"] = np.ones(periods) * np.inf

            if self.simulate_solar:
                upper_bounds["energy_solar"] = hourly_solar_profile
                if self.simulate_solar_exchange:
                    upper_bounds["energy_solar_sold"] = hourly_solar_profile

            if self.simulate_battery:
                upper_bounds["battery_charging"] = np.ones(periods)
                upper_bounds["battery_discharging"] = np.ones(periods)
                upper_bounds["energy_battery_charging"] = (
                    np.ones(periods) * battery_number*round(battery_max_power/battery_power_step_size)
                )
                upper_bounds["energy_battery_discharging"] = (
                    np.ones(periods) * battery_number*round(battery_max_power/battery_power_step_size)
                )

            if self.simulate_maximum_demand:
                upper_bounds["maximum_demand"] = [np.inf]

            self.bounds = [
                np.concatenate(list(self.bounds_dict[0].values())),
                np.concatenate(list(self.bounds_dict[1].values())),
            ]

        def build_integrality(self):
            self.integrality_dict = copy.deepcopy(self.null_dict)

            self.integrality_dict["pump_status"] = np.ones(periods)

            if self.simulate_battery:
                self.integrality_dict["battery_charging"] = np.ones(periods)
                self.integrality_dict["battery_discharging"] = np.ones(periods)
                self.integrality_dict["energy_battery_charging"] = np.ones(periods)
                self.integrality_dict["energy_battery_discharging"] = np.ones(periods)

            self.integrality = np.concatenate(list(self.integrality_dict.values()))
            self.integrality = [int(value) for value in self.integrality]

        build_null_decisions(self)
        build_costs(self)
        build_bounds(self)
        build_integrality(self)

    def build_constraints(self):

        def energy_balance_constraint(self):
            A, LB, LU = [], [], []
            for time in range(1, periods + 1):
                A0 = copy.deepcopy(self.null_dict)

                A0["pump_status"] = np.identity(periods)[time - 1] * pump_power * hours_delta * -1
                A0["energy_purchased"] = np.identity(periods)[time - 1]

                if self.simulate_solar:
                    A0["energy_solar"] = np.identity(periods)[time - 1]

                if self.simulate_battery:
                    A0["energy_battery_discharging"] = (np.identity(periods)[time - 1] * battery_power_step_size * hours_delta)
                    A0["energy_battery_charging"] = (np.identity(periods)[time - 1] * -1 / battery_roundtrip_eff * battery_power_step_size * hours_delta)

                lb0 = 0
                lu0 = 0

                A.append(np.concatenate(list(A0.values())))
                LB.append(lb0)
                LU.append(lu0)

            self.constraints.append([A, LB, LU])

        def reservoir_level_constraint(self):
            A, LB, LU = [], [], []
            for time in range(1, periods + 1):
                A0 = copy.deepcopy(self.null_dict)

                A0["pump_status"] = [
                    *np.repeat(pump_flow_rate, time) * hours_delta,
                    *np.zeros(periods - time),
                ]

                lb0 = (
                    reservoir_lb
                    + water_demand_cumulative[time - 1]
                    - reservoir_initial_level
                )
                lu0 = (
                    reservoir_ub
                    + water_demand_cumulative[time - 1]
                    - reservoir_initial_level
                )

                A.append(np.concatenate(list(A0.values())))
                LB.append(lb0)
                LU.append(lu0)
            
            if reservoir_initial_level:
                LB[-1] = water_demand_cumulative[periods - 1]

            self.constraints.append([A, LB, LU])

        def maximum_demand_constraint(self):
            A, LB, LU = [], [], []
            for time in range(0, periods):
                A0 = copy.deepcopy(self.null_dict)

                A0["energy_purchased"] = np.identity(periods)[time] / hours_delta
                A0["maximum_demand"] = [-1]

                lb0 = np.inf * -1
                lu0 = 0

                A.append(np.concatenate(list(A0.values())))
                LB.append(lb0)
                LU.append(lu0)

            self.constraints.append([A, LB, LU])

        def solar_exchange_balance(self):
            A, LB, LU = [], [], []
            for time in range(0, periods):
                A0 = copy.deepcopy(self.null_dict)

                A0["energy_solar"] = np.identity(periods)[time]
                A0["energy_solar_sold"] = np.identity(periods)[time]

                lb0 = 0
                lu0 = hourly_solar_profile[time]

                A.append(np.concatenate(list(A0.values())))
                LB.append(lb0)
                LU.append(lu0)

            self.constraints.append([A, LB, LU])

        def battery_charge_discharge_constraint(self):
            A, LB, LU = [], [], []
            for time in range(0, periods):
                A0 = copy.deepcopy(self.null_dict)

                A0["battery_charging"] = np.identity(periods)[time]
                A0["battery_discharging"] = np.identity(periods)[time]

                lb0 = 0
                lu0 = 1

                A.append(np.concatenate(list(A0.values())))
                LB.append(lb0)
                LU.append(lu0)

            self.constraints.append([A, LB, LU])

        def battery_coupling(self):
            for field in ["charging", "discharging"]:

                A, LB, LU = [], [], []
                for time in range(0, periods):
                    A0 = copy.deepcopy(self.null_dict)

                    A0[f"battery_{field}"] = np.identity(periods)[time] * battery_number * battery_max_power * hours_delta
                    A0[f"energy_battery_{field}"] = np.identity(periods)[time] * -1 * battery_power_step_size * hours_delta

                    lb0 = 0
                    lu0 = np.inf

                    A.append(np.concatenate(list(A0.values())))
                    LB.append(lb0)
                    LU.append(lu0)

                self.constraints.append([A, LB, LU])

        def battery_state_of_charge(self):
            A, LB, LU = [], [], []
            for time in range(1, periods + 1):
                A0 = copy.deepcopy(self.null_dict)

                A0["energy_battery_discharging"] = [
                    *np.ones(time),
                    *np.zeros(periods - time),
                ]

                A0["energy_battery_charging"] = [
                    *np.ones(time) * -1,
                    *np.zeros(periods - time),
                ]

                lb0 = (battery_min_level - battery_initial_level)*battery_number*battery_size/battery_power_step_size/hours_delta
                lu0 = (battery_max_level - battery_initial_level)*battery_number*battery_size/battery_power_step_size/hours_delta

                A.append(np.concatenate(list(A0.values())))
                LB.append(lb0)
                LU.append(lu0)

            if battery_mantain_level:
                LB[-1] = 0
                LU[-1] = 0

            self.constraints.append([A, LB, LU])

        self.constraints = []
        energy_balance_constraint(self)
        reservoir_level_constraint(self)

        if self.simulate_maximum_demand:
            maximum_demand_constraint(self)

        if self.simulate_solar_exchange:
            solar_exchange_balance(self)

        if self.simulate_battery:
            battery_charge_discharge_constraint(self)
            battery_coupling(self)
            battery_state_of_charge(self)

