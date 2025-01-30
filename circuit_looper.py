from circuit_reader import *
from circuit_calcs import *

# SETTINGS

# runner mode plots various variables in one circuit
# looper mode plots change of results of one circuit due to change of params
mode = "looper"  # "looper"  # 
template = "kent"  # "kent_equiv"  # "parallel"  # 
param_file = "params/kent_params_250128.txt"  # None  # "params/parallel_params_240709.txt"  # 
variation = "no_b2"  # "no_current_source"  # "add_idc"  # "biased_jj"  # 
results_to_watch = ["v(101)"]  # , "i(l0)", "i(l2)"]  # , "v(102)"]  # phase, leff, etc.

# looper settings
if mode == "looper":
    param_to_change = "phi1_mag"  # "l1_mag"  # "idc_mag"  # 
    # param_list = np.linspace(0e-9, 2.0e-9, 21)
    param_list = np.linspace(0.9999999999*np.pi, 1.0000000001*np.pi, 21)  # np.linspace(0e-6, 7.5e-6, 151)  # 
    results_list = {}
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
    cd.change_param("tran_step", 1e-13)
    cd.change_param("tran_stop", 5e-10)
    cd.change_param("tran_start", 3e-10)
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
    count = 1
    t1 = time.time()
    for param in param_list:
        print(param)
        if param == param_list[0]:
            cd.simulation_cycle(template, param_file, variation)  # needed to initialize
            # cd.change_param("i1_mag", 0)
            # cd.change_param("tran_step", 1e-13)
            # cd.change_param("tran_stop", 2e-9)  # 2e-10)
            # cd.change_param("tran_start", 0e-10)
            # cd.change_param("phi1_mag", 0)
        cd.change_param(param_to_change, param)
        cd.simulation_cycle(template, param_file, variation)
        time_array = cd.data["time"]
        for result in results_to_watch:
            result_value = cd.data[result].to_numpy()[-1]  # final value only for now
            results_list[result].append(result_value)
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
            
            '''
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
        plt.plot(param_list, results_list[result], "x", label=f"{result}")

    plt.xlabel(f"{param_to_change}")
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f"{result}_{param_to_change}.png")
    # plt.clf()
    plt.show()
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
    
# save collected data
'''
all_accels_list_pd = pd.DataFrame(all_accels_list)
all_accels_list_pd = all_accels_list_pd.sort_values("phase")
all_accels_list_pd.to_csv(f"{result}_accel_{param_to_change}.csv", header=True, index=False)
'''

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
