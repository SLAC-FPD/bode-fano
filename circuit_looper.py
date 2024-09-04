from circuit_reader import *
from circuit_calcs import *

# SETTINGS

# runner mode plots various variables in one circuit
# looper mode plots change of results of one circuit due to change of params
mode = "looper"  # "looper"  # 
template = "kent"  # "parallel"  # 
param_file = "params/kent_params_240731.txt"  # "params/parallel_params_240709.txt"  # 
variation = "no_b2"  # "pulsed_no_b2"  # "add_idc"  # "biased_jj"  # 
results_to_watch = ["v(101)"]  # , "i(l0)", "i(l2)"]  # , "v(102)"]  # phase, leff, etc.

# looper settings
if mode == "looper":
    param_to_change = "idc_mag"  # "idc_mag"
    param_list = np.linspace(0e-6, 7.5e-6, 751)
    results_list = {}
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
            cd.change_param("tran_step", 1e-13)
            cd.change_param("tran_stop", 2e-10)  # 2e-10)
            cd.change_param("tran_start", 0e-10)
            cd.change_param("phi1_mag", 0)
        cd.change_param(param_to_change, param)
        cd.simulation_cycle(template, param_file, variation)
        time_array = cd.data["time"]
        for result in results_to_watch:
            result_value = cd.data[result].to_numpy()[-1]  # final value only for now
            results_list[result].append(result_value)
            
            plt.plot(time_array, cd.data[result], label=f"{result}, idc_mag={param}")
            plt.tight_layout()
            plt.grid()
            plt.legend()
            plt.savefig(f"kent_test_nocap/{count}_idc_mag.png")
            plt.clf()
            
        count += 1
    t2 = time.time()
    print(f"Simulation loop took {(t2 - t1)/60} minutes.")
    
    for result in results_to_watch:
        plt.plot(param_list, results_list[result], "x", label=f"{result}")
    plt.xlabel(f"{param_to_change}")
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.show()
    
    '''
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
