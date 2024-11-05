from circuit_reader import *
from circuit_calcs import *

# SETTINGS

# runner mode plots various variables in one circuit
# looper mode plots change of results of one circuit due to change of params
mode = "runner"  # "looper"  # 
template = "series"  # "kent"  # 
param_file = None  # "params/kent_params_240708.txt"  # 
variation = None  # "no_b2"  # 
results_to_watch = ["v(101)"]  # , "i(l3)"]  # "i(l0)", "i(l1)", "i(l2)"]  # phase, leff, etc.

# looper settings
if mode == "looper":
    param_to_change = "idc_mag"  # "idc_mag"
    param_list = np.linspace(0.8194e-6, 0.81942e-6, 101)  # 0e-6, 5e-6, 101)
    results_list = {}
    for result in results_to_watch:
        results_list[result] = []

# Settings end, begin simulating results

# make CircuitData object
cd = CircuitData()
cd.simulation_cycle(template, param_file, variation)  # initialize

if mode == "runner":
    # simulate circuit
    cd.simulation_cycle(template, param_file, variation)
    time_array = cd.data["time"]
    for result in results_to_watch:
        plt.plot(time_array, cd.data[result], label=f"{result}")
        plt.xlabel("Time (s)")
        plt.tight_layout()
        plt.grid()
        plt.legend()
        plt.show()
    print(np.average(cd.data["v(101)"]))
    
elif mode == "looper":
    t1 = time.time()
    for param in param_list:
        print(param)
        cd.change_param(param_to_change, param)
        cd.simulation_cycle(template, param_file, variation)
        time_array = cd.data["time"]
        for result in results_to_watch:
            result_value = cd.data[result].to_numpy()[-1]  # final value only for now
            results_list[result].append(result_value)
    t2 = time.time()
    print(f"Simulation loop took {(t2 - t1)/60} minutes.")
    
    plt.plot(param_list, results_list[results_to_watch[0]], "x", label=f"{results_to_watch[0]}")
    plt.xlabel(f"{param_to_change}")
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.show()

'''
df = {"ibias": param_list, "phase": results_list[results_to_watch[0]], "i_l1": results_list[results_to_watch[1]], "i_l2": results_list[results_to_watch[2]], "i_l3": results_list[results_to_watch[3]]}
df = pd.DataFrame(df)
df.to_csv("Equiv_Weak.csv", index=False)
'''