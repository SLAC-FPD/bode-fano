from circuit_reader import *
from circuit_calcs import *

# SETTINGS

# runner mode plots various variables in one circuit
# looper mode plots change of results of one circuit due to change of params
mode = "runner"  # "looper"  # "runner"  # 
template = "ind_cancel"  # "impedance_match"  # "kent"  # "floquet"  # "kent_equiv"  # "parallel"  # 
param_file = "params/ind_cancel_params_250606.txt"  # None  #  # "params/kent_params_250303.txt"  # 
variation = None  # "couple_out"  # "biased_jj"  # "v_bias_cancel"  # "no_source"  # "no_b2"  # 
results_to_watch = ["v(101)"]  # , "i(la)", "i(lb)", "i(lout)"]  # , "v(102)"]  # , "i(l2)"]  # , "i(l1)", "i(l2)"]  # , "i(l0)"]  #   # phase, leff, etc.
save_results = False

# looper settings
if mode == "looper":
    param_to_change = "idc_mag"  # "lb_mag"  # "phi1_mag"  # "l1_mag"  # 
    param_list = np.linspace(1.0339169242309647e-05 * 0.998, 1.0339169242309647e-05 * 1.002, 41)
    # param_list = np.linspace(6.1e-10, 7e-10, 10)
    results_list = {param_to_change: param_list}
    results_list = {"ltot": [], "current_amp": [], "phase_diff": []}
    # all_accels_list = {"phase": [], "accel": [], "accel_mod": [], "init_val": []}
    for result in results_to_watch:
        results_list[result] = []

# Settings end, begin simulating results

# make CircuitData object
cd = CircuitData()

if mode == "runner":
    # simulate circuit
    cd.simulation_cycle(template, param_file, variation)  # needed to initialize
    # cd.change_param("phi1_mag", -4.8)
    # cd.change_param("idc_mag", 4e-6)  # avg to 4.5e-6
    cd.change_param("tran_step", 0.5e-11)
    cd.change_param("tran_stop", 1e-5)
    cd.change_param("tran_start", 0e-5)
    # cd.change_param("i1_max", 2*np.pi*cd.params["i1_freq"]*cd.params["tran_stop"])
    cd.simulation_cycle(template, param_file, variation)
    time_array = cd.data["time"]
    for result in results_to_watch:
        plt.plot(time_array, cd.data[result], label=f"{result}")
        plt.xlabel("Time (s)")
        plt.tight_layout()
        plt.grid()
        plt.legend()
        plt.show()
    # print(np.average(cd.data["v(101)"]))
    
elif mode == "looper":
    count = 1
    t1 = time.time()
    for param in param_list:
        print(param)
        if param == param_list[0]:
            cd.simulation_cycle(template, param_file, variation)  # needed to initialize
            # cd.change_param("i1_mag", 0)
            cd.change_param("tran_step", 0.5e-10)
            cd.change_param("tran_stop", 0.75e-6)  # 2e-10)
            cd.change_param("tran_start", 0e-10)
            # cd.change_param("phi1_mag", 0)
        cd.change_param(param_to_change, param)
        cd.simulation_cycle(template, param_file, variation)
        time_array = cd.data["time"]
        
        if template == "ind_cancel" and variation == "couple_out": current = cd.data["i(lout)"].to_numpy()[1000:]
        else: current = cd.data["i(lb)"].to_numpy()[1000:]
        voltage = cd.data["v(1)-v(0)"].to_numpy()[1000:]
        # jj_current = cd.data["i(l3)"].to_numpy()[19800:]
        # results_list["jj_current"].append(np.average(jj_current))
        # expected_jj_current = 1.6e-7 * np.sin(expected_phase)
        # expected_list["jj_current"].append(expected_jj_current)
        omega = cd.params["i1_freq"] * 2 * np.pi  # 1e12 for 0225 1e10 for 0303
        current_amp = (np.max(current) - np.min(current))/2
        voltage_amp = (np.max(voltage) - np.min(voltage))/2
        results_list["current_amp"].append(current_amp)
        cur_phase = time_array[np.where(current == np.max(current))[0][0]]
        vol_phase = time_array[np.where(voltage == np.max(voltage))[0][0]]
        phase_diff = (cur_phase - vol_phase)*omega
        # print(phase_diff)
        if phase_diff < -np.pi: phase_diff += 2*np.pi
        elif phase_diff > np.pi: phase_diff += -2*np.pi
        results_list["phase_diff"].append(phase_diff)
        # expected_ltot = 1e-9*(1 - param**2)*(lg + lj_)/(lg + (1-param**2)*lj_)
        # expected_list["ltot"].append(expected_ltot)
        ltot = voltage_amp / current_amp / omega
        # if phase_diff < 0: ltot = -ltot
        print(f"LTOT: {ltot} H")
        results_list["ltot"].append(ltot)
        
        for result in results_to_watch:
            result_value_ = cd.data[result].to_numpy()  # final value only for now
            result_value = np.average(result_value_[9900:])  # CHANGE THIS
            results_list[result].append(result_value)
            '''
            plt.plot(time_array, cd.data[result], label=f"{result}, {param_to_change}={param}")
            plt.tight_layout()
            plt.grid()
            plt.legend()
            plt.savefig(f"phi1_graphs/{count}_{param_to_change}.png")
            plt.clf()
            # setup
            delta_t = (time_array[1] - time_array[0])
            # cap = cd.data["@b1[cap]"][0]
            # cond = cd.data["@b1[g0]"][0]
            
            phase_x = cd.data[result].to_numpy()  # "position"
            phase_v = (phase_x[:-1] - phase_x[1:])/delta_t  # "velocity"
            phase_for_v = phase_x[:-1]
            phase_a = (phase_v[:-1] - phase_v[1:])/delta_t  # "acceleration"
            phase_for_a = phase_x[1:-1]
            phasev_for_a = (phase_v[:-1] + phase_v[1:])/2  # "average velocity"
            phase_a_mod = phase_a*cap - phasev_for_a*cond
            a, b = np.polyfit(phase_for_a, phase_a*cap, 1)
            x = phase_for_a
            fit = a*x + b
            # plt.plot(phase_for_a, phase_a*cap, "x", label=f"Acceleration, Initial {param_to_change}={param}")
            # plt.plot(x, fit, "r--", label=f"Fit, Initial {param_to_change}={param}")
            plt.plot(phase_for_a, phase_a*cap - fit, "x", label=f"Acceleration (Fit Removed), Initial {param_to_change}={param}")
            plt.plot(phase_for_a, phasev_for_a * cond, "x", label=f"Velocity, Initial {param_to_change}={param}")
            plt.plot(phase_for_a, phase_a_mod - fit, "x", label=f"Mod?, Initial {param_to_change}={param}")
            plt.tight_layout()
            plt.grid()
            plt.legend()
            plt.savefig(f"{param_to_change}_graphs/{count}_{param_to_change}_accel.png")
            plt.clf()
            phase_a_mod = phase_a - phasev_for_a*cond/cap  # normalize to phase_a levels
            for phase_ in phase_for_a:
                all_accels_list["phase"].append(phase_)
            for accel_ in phase_a:
                all_accels_list["accel"].append(accel_)
            for accel_mod_ in phase_a_mod:
                all_accels_list["accel_mod"].append(accel_mod_)
                all_accels_list["init_val"].append(param)
            '''
        count += 1
    t2 = time.time()
    print(f"Simulation loop took {(t2 - t1)/60} minutes.")
    
    for result in results_to_watch:
        phases = np.array(results_list[result]) % (2*np.pi)
        plt.plot(param_list, phases, "x", label=f"{result}")

    plt.xlabel(f"{param_to_change}")
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f"{result}_{param_to_change}.png")
    # plt.clf()
    plt.show()
    
    # '''
    plt.plot(param_list, results_list["ltot"], "x")  # , label=f"{result}")
    plt.xlabel(f"DC Current Bias (A)")
    plt.ylabel("Total Inductance Magnitude (H)")
    plt.grid()
    # plt.legend()
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(f"{result}_ltot_amps.png")
    # plt.clf()
    plt.show()
    # '''
    
    # '''
    plt.plot(param_list, results_list["phase_diff"], ".")  # , label=f"{result}")
    plt.xlabel(f"DC Current Bias (A)")
    plt.ylabel("Current-Voltage Phase difference (rad)")
    plt.grid()
    # plt.legend()
    plt.yscale("linear")
    plt.tight_layout()
    plt.savefig(f"{result}_ltot_phase_diff.png")
    # plt.clf()
    plt.show()
    # '''
    
    # '''
    plt.plot(param_list, results_list["current_amp"], ".")  # , label=f"{result}")
    plt.xlabel(f"DC Current Bias (A)")
    plt.ylabel("Current Amplitude in readout (A)")
    plt.grid()
    # plt.legend()
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(f"{result}_ltot_current_amp.png")
    # plt.clf()
    plt.show()
    # '''
    
    '''
    plt.plot(all_accels_list["phase"], all_accels_list["accel"], "x", label="Acceleration")
    plt.plot(all_accels_list["phase"], all_accels_list["accel_mod"], "x", label="Acceleration Mod")
    plt.xlabel("Phase")
    plt.ylabel("Acceleration")
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f"{result}_accel_{param_to_change}.png")
    plt.show()
    
    # ad-hoc
    l0, l1, k = 5e-10, 5e-10, 0.5
    # phases = results_list[result]
    phases = np.linspace(0, 2*np.pi, 721)
    xval = phases  # 
    ljs = calc_lj(1.1e-6, phases)
    expected_l = l0*(1 - k**2*l1/(l1+ljs))
    plt.plot(xval, expected_l, "x")
    plt.xlabel("Phase (rad)")
    plt.ylabel("Expected LTOT (H)")
    base_l = np.ones(len(phases)) * l0*(1 - k**2)
    plt.plot(xval, base_l, "--")
    plt.plot
    plt.xlim(1.8, 4.6)
    plt.ylim(-0.3e-9, 1.18e-9)
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.show()
    '''


# current = cd.data["i(lout)"].to_numpy()[1000:]
voltage = cd.data["v(1)-v(0)"].to_numpy()[1000:]
# voltage = cd.data["v(0)-v(6)"].to_numpy()[5000:]
if variation == "couple_out":
    current = cd.data["i(lout)"].to_numpy()[1000:]
    voltage_l0 = cd.data["v(0)-v(7)"].to_numpy()[1000:]
    voltage_l1 = cd.data["v(7)-v(2)"].to_numpy()[1000:]
else:
    current = cd.data["i(la)"].to_numpy()[1000:]
    voltage_l0 = cd.data["v(2)-v(1)"].to_numpy()[1000:]
    voltage_l1 = cd.data["v(0)-v(2)"].to_numpy()[1000:]
time_array2 = time_array[1000:]
fig, ax1 = plt.subplots()

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Voltage (V)')
ax1.plot(time_array2, voltage_l0, color="lightblue", label="Inductance to cancel")
ax1.plot(time_array2, voltage_l1, color="cyan", label="Inductance cancelling device")
ax1.plot(time_array2, voltage, color="blue", label="Coupled-in voltage")
ax1.tick_params(axis='y', labelcolor="blue")

ax2 = ax1.twinx()
ax2.set_ylabel('Current (A)')
ax2.plot(time_array2, current, color="orange")
ax2.tick_params(axis='y', labelcolor="orange")

ax1.legend()

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.grid()
plt.savefig(f"current_and_voltage_250314.png")
plt.show()

cur_amp = (np.max(current) - np.min(current)) / 2
vol_amp = (np.max(voltage) - np.min(voltage)) / 2
print(f"\nCURRENT AMPLITUDE: {cur_amp}, VOLTAGE AMPLITUDE: {vol_amp}")
omega = cd.params["i1_freq"] * 2 * np.pi
ltot = vol_amp / (cur_amp * omega)

period = (2*np.pi) / omega
cur_phase = time_array[np.where(current==np.max(current[50:]))[0][0]]
vol_phase = time_array[np.where(voltage==np.max(voltage[50:]))[0][0]]
phase_diff = (cur_phase - vol_phase) % period / period * 2 * np.pi
if phase_diff > np.pi: phase_diff = phase_diff - 2 * np.pi

print(f"EFFECTIVE INDUCTANCE MAGNITUDE: {ltot}, PHASE LAG: {phase_diff}\n")

# save collected data
'''
all_accels_list_pd = pd.DataFrame(all_accels_list)
all_accels_list_pd = all_accels_list_pd.sort_values("phase")
all_accels_list_pd.to_csv(f"{result}_accel_{param_to_change}.csv", header=True, index=False)
'''

if save_results:
    if mode == "runner": data_pd = pd.DataFrame(cd.data)
    elif mode == "looper": data_pd = pd.DataFrame(results_list)
    data_pd.to_csv(f"{template}_data_{mode}.csv", header=True, index=False)

# FIT zone
'''
plt.plot(param_list, results_list[result], "x", label=f"{result}")
x_vals = param_list * 1.e9
y_vals = results_list[result]
model = exp_inverse_model  # exp_alt_model  # 
model_fit = model.fit(y_vals, x=x_vals, a=1, b=2.057, c=0, d=0)
plt.plot(param_list, model_fit.best_fit, label="Fit")

plt.xlabel(f"{param_to_change}")
plt.tight_layout()
plt.grid()
plt.legend()
# plt.savefig(f"{result}_{param_to_change}.png")
plt.show()

print(model_fit.fit_report())
'''
