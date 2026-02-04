# newly written for bfe template. use old_ic_looper for other templates

from circuit_reader import *
from circuit_calcs import *

cd = CircuitData()  # create instance.

# options
template = "bfe"
param_file = "params/bfe_params_251218.txt"  # None  # 
variation = "parallel_bfe_pulsed"  # "series_bfe" # _resonant"  # _no_bfe"  #  _bfe" # None  # "simple_lc"  # 
mode = ""
if variation is not None: mode += f"_{variation}"

loop_param = "v1_freq"  # "l2_mag"  # "idc_mag"  # "c1_offset"  # 
loop_list = [5.e7]  # np.linspace(3e7, 7e7, 41)  # np.linspace(3.4489e-10, 3.4491e-10, 21)  # np.linspace(1.0339169242309647e-05*0.99999, 1.0339169242309647e-05*1.00001, 21)  # 
loop_len = len(loop_list)

# If looping over multiple params, use a dedicated file. This overrides settings above.
loop_params_file = None  # "ic_results/loop_params.csv"  # 
if loop_params_file is not None:
    loop_df = pd.read_csv(loop_params_file)
    loop_param = list(loop_df.keys())
    loop_list = loop_df.to_numpy()
    loop_len = len(loop_df[loop_param[0]])

# initialize
cd.read_template(template, param_file, variation)
cd.change_param("filename", cd.params["filename"] + mode)

# time settings here
freq = cd.params["v1_freq"]
round_digits = np.floor(np.log10(freq)) + 1
step_time = 10**(-(round_digits+3)) * 10 # / 100  # * 10  # multiply or divide to this value
idx_ringing = 1000000  # 0  # 100000  # to remove beginning of simulation. 2e4 is default, 1e5 for paper
npts = 1000000 # 1e5 is default, 2e3 with step_time/100 for JJ ringing
# for resonant: step_time <=*10, npts >= 1e5
if "resonant" in variation: step_time *= 100
start_time = 0  # change if needed
stop_time = step_time * (npts + idx_ringing)
cd.change_param("tran_step", step_time)
cd.change_param("tran_stop", stop_time)  # 2e-10)
cd.change_param("tran_start", 0e-10)

# results found in these components
time_tag = "time"
iout_tag = "i(lout)"
iin_tag = "i(lin)"
iloop_tag = "i(l0)"
ijj_tag = "i(l2)"
vin_tag = "v(1)-v(0)"
vbase_tag = "v(2)-v(0)"
vout_tag = "v(8)-v(0)"
vcancel_tag = "v(3)-v(2)"
vjj_tag = "v(5)-v(4)"
phase_tag = "v(101)"

if variation == "simple_lc":
    iout_tag = "i(l0)"
    iin_tag = "i(l0)"
    # vcancel_tag = "v(3)-v(2)"
    vbase_tag = "v(2)-v(1)"
    vout_tag = "v(0)-v(3)"
    
if "parallel" in variation:
    iout_tag = "i(lout)"
    vcancel_tag = "v(2)-v(0)"

if "series" in variation:
    iout_tag = "i(l0)"
    iin_tag = "i(l0)"
    vcancel_tag = "v(3)-v(2)"
    vbase_tag = "v(2)-v(8)"

cur_loop_idx = 0
results_dict = {f"{loop_param}": [], "pin_sum": [], "pout_sum": []}  # let's add to the dictionary if needed
if variation is None or ("_bfe" in variation and "_no" not in variation): results_dict["phase_avg"] = []  # add phase to dictionary
t1 = time.time()

for loop_val in loop_list:
    results_dict[f"{loop_param}"].append(loop_val)
    if loop_param is not None:
        cur_loop_idx += 1
        print(f"\n({cur_loop_idx}/{loop_len}): {loop_param} = {loop_val}")
    # change params
    if type(loop_param) == list:  # put in multiple values to change at once
        for loop_idx in range(len(loop_param)):
            cd.change_param(loop_param[loop_idx], loop_val[loop_idx])
    else:
        cd.change_param(loop_param, loop_val)
    # run the simulation
    
    cd.simulation_cycle(template, None, variation)
    time_array = cd.data[time_tag].to_numpy()[idx_ringing:]
    iout_array = cd.data[iout_tag].to_numpy()[idx_ringing:]
    iin_array = cd.data[iin_tag].to_numpy()[idx_ringing:]
    iloop_array = cd.data[iloop_tag].to_numpy()[idx_ringing:]
    vin_array = cd.data[vin_tag].to_numpy()[idx_ringing:]
    vbase_array = cd.data[vbase_tag].to_numpy()[idx_ringing:]
    vout_array = cd.data[vout_tag].to_numpy()[idx_ringing:]
    vcancel_array = cd.data[vcancel_tag].to_numpy()[idx_ringing:]
    pin_array = iin_array * vin_array
    pout_array = iout_array * vout_array  # pout_array2 = vout_array ** 2 / cd.params["rout_mag"]
    p_idx_base = int(np.round(1 / freq / step_time))  # one frequency cycle
    p_idx = p_idx_base
    while p_idx_base + p_idx < len(time_array): p_idx = p_idx_base + p_idx  # assumes steady state
    # p_step = step_time*p_idx
    pin_sum = np.abs(np.sum(pin_array[-p_idx:]) / p_idx)
    pout_sum = np.abs(np.sum(pout_array[-p_idx:]) / p_idx)
    '''
    pin_freq, pin_power = do_fft(time_array, vin_array*iin_array)
    pout_freq, pout_power = do_fft(time_array, vout_array*iout_array)
    try: comp_idx = np.where(np.round(pin_freq) == 2*freq)[0][0]  # 0  # because power multiplies two sinusoids
    except IndexError: comp_idx = np.where(pin_power == np.max(pin_power))[0][0]
    start_idx = np.min([1, comp_idx-5])
    pfreq_fft = pin_freq[comp_idx]
    pin_sum = pin_power[comp_idx]
    pout_sum = pout_power[comp_idx]
    '''
    print(pin_sum, pout_sum, power2dB(pout_sum/pin_sum))
    results_dict["pin_sum"].append(pin_sum)
    results_dict["pout_sum"].append(pout_sum)
    
    if loop_len == 1:
        # Voltage plot
        # plt.plot(time_array, vin_array, label="Input voltage")
        plt.plot(time_array, vbase_array, label="Inductance to cancel")
        if "resonant" in variation or variation == "simple_lc": plt.plot(time_array, vcancel_array, label="Capacitor")
        elif variation is None or variation == "parallel_bfe": plt.plot(time_array, vcancel_array, label="Inductance canceled")
        # plt.plot(time_array, vout_array, label="Output voltage")
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage (V)")
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"ic_results/ic_voltage{mode}.png")
        plt.show()
        '''
        # Current plot
        plt.plot(time_array, iin_array - offset(iin_array), zorder=1, label="Input (centered)")
        plt.plot(time_array, iout_array - offset(iout_array), zorder=2, label="Output (centered)")
        plt.legend()
        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")
        plt.tight_layout()
        plt.grid()
        plt.savefig(f"ic_results/ic_current{mode}.png")
        plt.show()
        # Power plot
        plt.plot(time_array, pin_array, label="Input")
        plt.plot(time_array, pout_array, label="Output")
        plt.legend()
        plt.xlabel("Time (s)")
        plt.ylabel("Instantaneous Power (W)")
        plt.tight_layout()
        plt.grid()
        plt.savefig(f"ic_results/ic_power{mode}.png")
        plt.show()
        '''
        
    if variation is None or ("_bfe" in variation and "_no" not in variation):  # JJ related things are only added if they exist.
        ijj_array = cd.data[ijj_tag].to_numpy()[idx_ringing:]
        vjj_array = cd.data[vjj_tag].to_numpy()[idx_ringing:]
        phase_array = cd.data[phase_tag].to_numpy()[idx_ringing:]
        phase_avg = np.average(phase_array)
        results_dict["phase_avg"].append(phase_avg)
        print(f"Average phase: {phase_avg}")
        if loop_len <= 1:
            # JJ phase plot
            plt.plot(time_array, phase_array)
            plt.xlabel("Time (s)")
            plt.ylabel("Josephson junction phase")
            plt.grid()
            plt.tight_layout()
            plt.savefig(f"ic_results/ic_phase{mode}.png")
            plt.show()

t2 = time.time()
print(f"Full loop took {(t2 - t1)/60} minutes.")

# plot power transfer if looped among many sims
if loop_len > 1:
    pout_sum_list = np.array(results_dict["pout_sum"])
    pin_sum_list = np.array(results_dict["pin_sum"])
    if loop_params_file is not None: xparam = np.linspace(0, len(pout_sum_list)-1, len(pout_sum_list))
    else: xparam = loop_list
    plt.plot(xparam, power2dB(pout_sum_list/pin_sum_list), ".")
    if loop_param == "v1_freq": loop_xlabel = "Frequency (Hz)"
    elif loop_params_file is not None: loop_xlabel = "Configuration #"
    else: loop_xlabel = loop_param
    plt.xlabel(loop_xlabel)  # change if needed
    plt.ylabel("Power Transfer (dB)")
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"ic_results/ic_transfer{mode}.png")
    plt.show()
    #'''
    pout_sum_list = np.array(results_dict["pout_sum"])
    pin_sum_list = np.array(results_dict["pin_sum"])
    plt.plot(xparam, pin_sum_list, ".", label="Input Power")
    plt.plot(xparam, pout_sum_list, ".", label="Output Power")
    plt.xlabel(loop_xlabel)  # change if needed
    plt.ylabel("Average Power (W)")
    plt.yscale("log")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"ic_results/ic_power_sweep{mode}.png")
    plt.show()
    #'''
    if variation is None or variation == "parallel_bfe":
        phase_list = np.array(results_dict["phase_avg"])
        plt.plot(xparam, phase_list, ".", label="Input Power")
        plt.xlabel(loop_xlabel)  # change if needed
        plt.ylabel("Average Phase (rad)")
        plt.grid()
        plt.tight_layout()
        plt.savefig(f"ic_results/ic_phase_sweep{mode}.png")
        plt.show()
    # save data, for looped results
    for key in results_dict.keys():  # remove empty, make sure first value is not empty
        if len(results_dict[key]) == 0: results_dict[key] = np.zeros(keylen)
        else: keylen = len(results_dict[key])
    print("Saving Results.")
    data_pd = pd.DataFrame(results_dict)  # if loops
    data_pd.to_csv(f"ic_results/ic_results{mode}_data.csv", header=True, index=False)
    
# OLD STUFF
'''
    if plot_jj_info and draw_plots:
        fig, ax1 = plt.subplots()

        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Josephson Junction Voltage (V)', color="blue")
        ax1.plot(time_array, vjj, color="blue")
        ax1.tick_params(axis='y', labelcolor="blue")

        ax2 = ax1.twinx()
        ax2.set_ylabel('Josephson Junction Current (A)', color="orange")
        ax2.plot(time_array, ijj, color="orange")
        # ax2.plot(time_array, ics*np.sin(phase), "r--")
        ax2.tick_params(axis='y', labelcolor="orange")

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.grid()
        plt.savefig(f"ic_results/ic_jj_info.png")
        plt.show()
        plt.cla()

results_dict = {f"{loop_param_name}": [], "params_list": [], "phase_avg": [], "ijj_mag": [], "vjj_mag": [], 
            "iin_mag": [], "iloop_mag": [], "iout_mag": [], 
            "vinput_mag": [], "vin_mag": [], "vbase_mag": [], "vcancel_mag": [], "vout_mag": [], "voutput_mag": [],
            "leff_input": [], "energy_in": [], "energy_out": [], "energy_out_prev": [],
            "ifreq_fft": [], "iin_fft": [], "iout_fft": [], "iout_prev_fft": [], 
            "energy_in_fft": [], "energy_out_fft": [], "energy_out_prev_fft": [],
            "vfreq_fft": [], "vin_fft": [], "vout_fft": [], "pfreq_fft":[] , "pin_fft": [], "pout_fft": []}
no_bodefano = False  # True  # 
no_r_jj = False
no_c_jj = False
phase_bias = False
original_params = True  # False  # get params directly from params file
plot_jj_info = True
get_fft_data = True
draw_plots = False  # True  # 


'''
