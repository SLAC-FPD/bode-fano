from circuit_reader import *
from circuit_calcs import *
from scipy.optimize import fsolve

# SETTINGS

# runner mode plots various variables in one circuit
# looper mode plots change of results of one circuit due to change of params
mode = "looper"  # "looper"  # 
template = "kent"
param_file = "params/kent_params_250219.txt"  # None  # "params/parallel_params_240709.txt"  # 
variation = "no_b2_v"  # "no_current_source"  # "add_idc"  # "biased_jj"  # 
# results_to_watch = ["v(101)"]  # , "i(l0)", "i(l2)"]  # , "v(102)"]  # phase, leff, etc.

lj = calc_lj(1.6e-7, 0)
func = lambda phi: lg/lj*np.sin(np.pi - phi) + phi + np.pi

withcap = True
withcap_ = ""
if withcap: withcap_ = "_withcap"

# looper settings
if mode == "looper":
    lg = 3.5e-9
    expected_phase = np.abs(fsolve(func, np.pi/2)[0])
    lj_ = calc_lj(1.6e-7, expected_phase)
    param_to_change = "k1_mag"  # "phi1_mag"
    param_list = np.linspace(0.001, 0.999, 999)
    expected_list = {"jj_phase": np.ones(len(param_list))*expected_phase, "ltot": [], "jj_current": []}
    results_list = {"jj_phase": [], "current_amplitude": [], "ltot": [], "phase_diff": [], "jj_current": []}

# Settings end, begin simulating results

# make CircuitData object
cd = CircuitData()
cd.simulation_cycle(template, param_file, variation)  # needed to initialize

if mode == "runner":
    print("Not used")
    
elif mode == "looper":  # fixed vals: l0_mag, i1_freq, ics1_mag, cpic_mag
    t1 = time.time()
    for param in param_list:
        print(param)
        cd.change_param(param_to_change, param)
        l2_val = lg/(1 - param**2)
        cd.change_param("l1_mag", l2_val)
        if withcap: cd.change_param("cpic_mag", 0.7e-9)
        cd.simulation_cycle(template, param_file, variation)
        time_array = cd.data["time"]
        jj_phase = cd.data["v(101)"].to_numpy()[19800:]  # final value only for now
        results_list["jj_phase"].append(np.average(jj_phase))
        time_ = time_array.to_numpy()[19800:]
        current = cd.data["i(l0)"].to_numpy()[19800:]
        voltage = cd.data["v(1)-v(0)"].to_numpy()[19800:]
        jj_current = cd.data["i(l2)"].to_numpy()[19800:]
        results_list["jj_current"].append(np.average(jj_current))
        expected_jj_current = 1.6e-7 * np.sin(expected_phase)
        expected_list["jj_current"].append(expected_jj_current)
        omega = 1e12
        current_amp = (np.max(current) - np.min(current))/2
        voltage_amp = (np.max(voltage) - np.min(voltage))/2
        results_list["current_amplitude"].append(current_amp)
        cur_phase = time_[np.where(current == np.max(current))[0][0]]
        vol_phase = time_[np.where(voltage == np.max(voltage))[0][0]]
        phase_diff = (cur_phase - vol_phase)*omega
        # print(phase_diff)
        if phase_diff < -np.pi: phase_diff += 2*np.pi
        elif phase_diff > np.pi: phase_diff += -2*np.pi
        results_list["phase_diff"].append(phase_diff)
        expected_ltot = 1e-9*(1 - param**2)*(lg + lj_)/(lg + (1-param**2)*lj_)
        expected_list["ltot"].append(expected_ltot)
        ltot = voltage_amp / current_amp / omega
        if phase_diff < 0: ltot = -ltot
        results_list["ltot"].append(ltot)

    t2 = time.time()
    print(f"Simulation loop took {(t2 - t1)/60} minutes.")

for result in results_list:
    try: plt.plot(param_list, expected_list[result], "x", label=f"Expected {result}")
    except KeyError: pass
    plt.plot(param_list, results_list[result], ".", label=f"WRSPICE {result}")
    if result == "ltot":
        max_val = 1.2*np.max(results_list["ltot"])
        min_val = 1.2*np.min(results_list["ltot"])
        if min_val > 0: min_val = np.min(results_list["ltot"]) - 0.2*np.max(results_list["ltot"])
        plt.ylim(min_val, max_val)
    plt.xlabel("Coupling constant")
    plt.ylabel(f"{result}")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"lg_{lg}_{result}{withcap_}.png")
    plt.show()
    plt.cla()
    