from circuit_reader import *
from circuit_calcs import *

# SETTINGS

template = "floquet"  # "snail"  # "kent"  # "parallel"  # 
param_file = None  # "params/kent_params_240905.txt"  # "params/parallel_params_240709.txt"  # 
variation = None  # "no_b2"  # "pulsed_no_b2"  # "add_idc"  # "biased_jj"  # 
results_to_watch = ["v(101)", "i(l1)"]  # , "v(102)"]  # phase, leff, etc.

# Settings end, begin simulating results

# make CircuitData object
cd = CircuitData()

# y-axis
param_name1 = "idc_mag"  # "cpic_mag" # "phi1_mag"  # "vm_mag"  # 
param_list1 = np.linspace(0, 1.0339169242309647e-05*3, 101)  # 0, 2*np.pi, 73)  # np.linspace(8.e-3, 100e-3, 185)  # np.logspace(-10, -5, 51)

# x-axis
param_name2 = "i1_mag"  # "i1_freq"  # 
param_list2 = np.linspace(0e-5, 5e-6, 251)  # np.linspace(2.5e6, 3e6, 51)  # 

result1_2d = []
result2_2d = []

cd.simulation_cycle(template, param_file, variation)  # needed to initialize
# cd.change_param("tran_step", 1e-13)
# cd.change_param("tran_stop", 5e-10)
# cd.change_param("tran_start", 3e-10)

# Fixed parameter (Frequency if using i1_freq)
# fixed_param = "i1_freq"  # "idc_mag"  # 
# cd.change_param(fixed_param, 1.e6)  #  4.15e-6)  # 

t1 = time.time()
for param1 in param_list1:
    result1_column = []
    result2_column = []
    cd.change_param(param_name1, param1)
    for param2 in param_list2:
        print(param1, param2)
        cd.change_param(param_name2, param2)
        cd.simulation_cycle(template, param_file, variation)
        time_array = cd.data["time"]

        '''
        for result in results_to_watch:
            plt.plot(time_array, cd.data[result], label=f"{result}")
            plt.xlabel("Time (s)")
            plt.tight_layout()
            plt.grid()
           plt.legend()
           plt.show()
        '''
        # set up results
        result1_ = cd.data[results_to_watch[0]]
        result1 = np.average(result1_[9500:])%(2*np.pi)
        result2_ = cd.data[results_to_watch[1]]
        result2 = np.max(result2_[9500:])  #  - np.min(result2_[9500:]))/2
        print(result2)
        
        result1_column.append(result1)
        result2_column.append(result2)
    result1_2d.append(result1_column)
    result2_2d.append(result2_column)
        
t2 = time.time()
print(f"Simulation loop took {(t2 - t1)/60} minutes.")
    
plt.imshow(result1_2d, extent=(param_list2[0], param_list2[-1], param_list1[-1], param_list1[0]), interpolation='none', aspect="auto", cmap="seismic")  # , norm=colors.LogNorm())
plt.xlabel(f"{param_name2}")
plt.ylabel(f"{param_name1}")
# plt.yscale("log")

cbar = plt.colorbar()
# plt.clim(-5, 5)
# plt.clim(-5e-6, 5e-6)
cbar.set_label("Average Phase", rotation=90)

plt.tight_layout()
plt.savefig(f"2d_plot_{param_name1}_{param_name2}_{results_to_watch[0]}.png")
plt.show()


plt.imshow(result2_2d, extent=(param_list2[0], param_list2[-1], param_list1[-1], param_list1[0]), interpolation='none', aspect="auto")  # , cmap="seismic")  # , norm=colors.LogNorm())
plt.xlabel(f"{param_name2}")
plt.ylabel(f"{param_name1}")
# plt.yscale("log")

cbar = plt.colorbar()
# plt.clim(-5, 5)
# plt.clim(-5e-6, 5e-6)
cbar.set_label("Loop circuit amplitude")

plt.tight_layout()
plt.savefig(f"2d_plot_{param_name1}_{param_name2}_{results_to_watch[1]}.png")
plt.show()

