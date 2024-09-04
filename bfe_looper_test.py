from circuit_reader import *
from circuit_calcs import *

# SETTINGS

template = "kent"  # "parallel"  # 
param_file = "params/kent_params_240731.txt"  # "params/parallel_params_240709.txt"  # 
variation = "no_b2"  # "pulsed_no_b2"  # "add_idc"  # "biased_jj"  # 
results_to_watch = ["v(101)", "i(l0)", "i(l2)"]  # , "v(102)"]  # phase, leff, etc.

# Settings end, begin simulating results

# make CircuitData object
cd = CircuitData()

# average of dc bias amplitude
param_name1 = "i1_mag"
param_list1 = np.logspace(-10, -5, 51)

# amplitude of oscillations fed (if using idc_mag)
param_name2 = "idc_mag"  # "i1_freq"  # 
param_list2 = np.linspace(0.e-6, 7.5e-6, 76)  # np.linspace(2.5e11, 3e11, 51)  # 

phase_results_2d = []
amp_il2_results_2d = []

cd.simulation_cycle(template, param_file, variation)  # needed to initialize
cd.change_param("tran_step", 1e-13)
cd.change_param("tran_stop", 5e-10)
cd.change_param("tran_start", 3e-10)

# Fixed parameter (Frequency if using i1_freq)
fixed_param = "i1_freq"  # "idc_mag"  # 
cd.change_param(fixed_param, 2.5e11)  # 2.77e11)  # 4.15e-6)  # 

t1 = time.time()
for param1 in param_list1:
    phase_results_column = []
    amp_il2_results_column = []
    cd.change_param(param_name1, param1)
    for param2 in param_list2:
        print(param1, param2)
        cd.change_param(param_name2, param2)
        cd.simulation_cycle(template, param_file, variation)
        time_array = cd.data["time"]
        phase = cd.data[results_to_watch[0]]
        il0 = cd.data[results_to_watch[1]]
        il2 = cd.data[results_to_watch[2]]
        '''
        for result in results_to_watch:
            plt.plot(time_array, cd.data[result], label=f"{result}")
            plt.xlabel("Time (s)")
            plt.tight_layout()
            plt.grid()
           plt.legend()
           plt.show()
        '''
        phase_result = np.average(cd.data["v(101)"])
        phase_results_column.append(phase_result)
        amp_il2 = cd.data["i(l2)"] - np.average(cd.data["i(l2)"]) # centered at zero
        amp_il2_result = np.sqrt(2*np.average(amp_il2**2)) / param1  # rms avg, normalized (transfer-fn like)
        amp_il2_results_column.append(np.log(amp_il2_result))  # log scale
    phase_results_2d.append(phase_results_column)
    amp_il2_results_2d.append(amp_il2_results_column)
        
t2 = time.time()
print(f"Simulation loop took {(t2 - t1)/60} minutes.")
    
plt.imshow(phase_results_2d, extent=(param_list2[0], param_list2[-1], param_list1[-1], param_list1[0]), interpolation='none', aspect="auto", cmap="seismic")  # , norm=colors.LogNorm())
plt.xlabel(f"{param_name2}")
plt.ylabel(f"{param_name1}")
plt.yscale("log")

cbar = plt.colorbar()
# plt.clim(-5, 5)
# plt.clim(-5e-6, 5e-6)
cbar.set_label("Average Phase", rotation=90)

plt.tight_layout()
plt.savefig(f"bfe_phase.png")
plt.show()

plt.imshow(amp_il2_results_2d, extent=(param_list2[0], param_list2[-1], param_list1[-1], param_list1[0]), interpolation='none', aspect="auto")  # , cmap="seismic")  # , norm=colors.LogNorm())
plt.xlabel(f"{param_name2}")
plt.ylabel(f"{param_name1}")
plt.yscale("log")

cbar = plt.colorbar()
# plt.clim(-5, 5)
# plt.clim(-5e-6, 5e-6)
cbar.set_label("Relative loop circuit amplitude")

plt.tight_layout()
plt.savefig(f"bfe_amp.png")
plt.show()
