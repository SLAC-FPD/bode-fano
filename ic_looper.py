from circuit_reader import *
from circuit_calcs import *

cd = CircuitData()

# to do list?
'''
- Send results to separate folder
- Make loopable
'''

# what do we want to do
template = "ind_cancel"  # "kent"  # 
param_file = "params/ind_cancel_params_250606.txt"  # "params/kent_params_250303.txt"  # None  # 
variation = "voltage_source_resonant"  # "couple_out_rout"  # _resonant"  # "voltage_source_output"  # "v_bias_cancel"  # None  # 
no_bodefano = False  # True  # 
no_r_jj = False
no_c_jj = False
phase_bias = False
original_params = False  # True  # get params directly from params file
plot_jj_info = True
get_fft_data = True
draw_plots = False  # True  # 
mode = ""

loop_param = "i1_freq"  # None  # "l1_mag"  # None  # triggers looper if not None
loop_list = np.linspace(1.e6, 1e9, 334)  # None  # np.linspace(3.333e-10, 3.335e-10, 21)  # None  # np.linspace(3.9e7, 4.0e7, 11)  # 

if "_resonant" in variation: no_bodefano = True  # for now, jj has convergence problem
if no_bodefano:
    plot_jj_info = False  # I don't need it
    mode += "_nobodefano"
if no_r_jj: mode += "_norjj"
if no_c_jj: mode += "_nocjj"
if phase_bias: mode += "_phasebias"
if loop_param is None: loop_list = [None]

# initialize
cd.read_template(template, param_file, variation)
# results_to_watch = ["v(101)"]
save_results = True

# set base parameters
ics = 1.09e-6  # 1.6e-7  # This is the actual jj junction critical current
ic = ics  # alias
icr = ics
lj0 = calc_lj(ics)
lj = lj0  # calc_lj(ics)
omega = 1.e7 * 2*np.pi  # 
freq = omega/(2*np.pi)
round_digits = np.floor(np.log10(freq)) + 1
idc = 1.0339169242309647e-05  # 0  # to current bias to pi
if idc == 0: mode += "_nojjbias"

step_time = 10**(-(round_digits+3))
idx_ringing = 20000  # 2000  # if needed to remove
npts = 100000  # 1000000 if need to go three digits accurate
stop_time = step_time * npts + idx_ringing * step_time

l0 = 3e-10  # 0  # 5e-10  # 
kin = np.sqrt(1./3.)  # used when coupling in and out
kout = np.sqrt(1./3.)  # kin  # these are free-ish parameters

linput = 0  # 0  # input inductance
loutput = 1e-2  # 2.25e-10  # 0  # output inductance in readout

k = 0.2  # between coupled circuits

###########################################################################
# THINGS THAT DON'T NEED MUCH CHANGE

if "couple_out" in variation:  
    lin = l0/2/(1 - kin**2)
    lout = l0/2/(1 - kout**2)
    cin = linput/lin
    cout = loutput/lout
    c = - 2  # + kout**2/(cout + 1)  # + kin**2  # - 2
    c_ = 1
    
    a1 = (c_ + 1)*(2 - kin**2/(cin + 1) - kout**2/(cout + 1))  # 4  # 
    a2 = c_*(2 - kin**2/(cin + 1) - kout**2/(cout + 1))  # 2  # 
    b = c_ - 1 + k**2  # 0.04  # 
    
    det_x0 = (k**2*a2 - a1 + b*c)**2 - 4*c*a1*(c_ - b)
    x0 = ((k**2*a2 - a1 + b*c) - np.sqrt(det_x0) ) / (2*c*a1)
    xj = (c*x0 + k**2 - 1)/(c*x0 - 1)
    
    l0 = 0
    lin_ = lin
    
else:
    lin = l0
    # linput = lin
    lout = l0
    x0 = 1.0
    xj = 1 - k**2/2
    lin_ = lin*(1 - kin**2*lout/(lout + loutput))
    
# print(f"Lin/L1 = {x0}, LJ/L2 = {xj}")
l1 = (l0 + lin_) / x0  # 3e-10  # 5e-10  # 
l2 = lj / xj  # lj * 1.0296298147012206  # lj/(1-k**2/2)  # 2.1178582977171555e-09  # 

# do a solver later with a lg, k, and ltot that you want
lg = l2 * (1 - k**2*l1/(l0 + l1))
ltot = l0 + l1 * (1 - k**2*l2/(l2-lj))
# print(f"LJ: {lj}, LG: {lg}, LTOT: {ltot}, L0: {l0}")

if original_params:  # does not change time scale, i1_freq/mag
    icr = cd.params["icrit_mag"]
    ics = cd.params["ics1_mag"]
    idc = cd.params["idc_mag"]
    l0 = cd.params["l0_mag"]
    l1 = cd.params["l1_mag"]
    l2 = cd.params["l2_mag"]
    k = cd.params["k1_mag"]
    kin = cd.params["k0_mag"]
    kout = cd.params["k0_mag"]
    lin = cd.params["la_mag"]
    lout = cd.params["lc_mag"]
    lj0 = calc_lj(ics)
    lj = lj0  # calc_lj(ics)
    
# change parameters from base
# cd.change_param("cpic_mag", 0)
# cd.change_param("rtype", 0)

cd.change_param("tran_step", step_time)
cd.change_param("tran_stop", stop_time)  # 2e-10)
cd.change_param("tran_start", 0e-10)
cd.change_param("i1_freq", freq)
# cd.change_param("i1_mag", 1e-12)
cd.change_param("icrit_mag", icr)
cd.change_param("ics1_mag", ics)
cd.change_param("idc_mag", idc)
cd.change_param("l0_mag", l0)
cd.change_param("l1_mag", 3.3339e-10)  # l1)
cd.change_param("l2_mag", l2)
cd.change_param("k1_mag", k)
cd.change_param("filename", f"ic_results/IndCancel_{variation}{mode}")

iout_tag = "i(lb)"
iin_tag = "i(la)"
iloop_tag = "i(lb)"
ijj_tag = "i(l2)"
vin_tag = "v(1)-v(0)"
vbase_tag = "v(2)-v(1)"
vout_tag = "v(2)-v(1)"
voutput_tag = "v(2)-v(1)"
vcancel_tag = "v(0)-v(2)"
vinput_tag = "v(0)-v(6)"
vjj_tag = "v(4)-v(3)"
phase_tag = "v(101)"

if "v_bias_cancel" in variation:
    iout_tag = "i(l1)"
    iin_tag = "i(l1)"
    iloop_tag = "i(l1)"
    vin_tag = "v(5)-v(1)"
    vbase_tag = "v(5)-v(1)"
    vcancel_tag = "v(0)-v(5)"
    vinput_tag = "v(1)-v(0)"
    phase_tag = "v(101)"
        
elif "couple_out" in variation:  # variation == "couple_out" or variation == "couple_out_load" or variation == "couple_out_open":
    cd.change_param("l0_mag", 0)
    cd.change_param("la_mag", lin)
    cd.change_param("lb_mag", lin)
    cd.change_param("k0_mag", kin)
    cd.change_param("lc_mag", lout)
    cd.change_param("ld_mag", lout)
    cd.change_param("k3_mag", kout)
    cd.change_param("lin_mag", linput)
    cd.change_param("lout_mag", loutput)
    vcancel_tag = "v(7)-v(2)"
    vout_tag = "v(0)-v(7)"
    voutput_tag = "v(0)-v(8)"
    iout_tag = "i(ld)"  # "i(lout)"  # 
    
elif "voltage_source" in variation:
    cd.change_param("l0_mag", l0)
    cd.change_param("la_mag", lin)
    cd.change_param("lb_mag", lin)
    cd.change_param("k0_mag", kin)
    cd.change_param("lout_mag", loutput)
    vcancel_tag = "v(6)-v(2)"
    vinput_tag = "v(1)-v(0)"
    vout_tag = "v(0)-v(6)"
    voutput_tag = "v(0)-v(7)"
    iin_tag = "i(la)"
    iloop_tag = "i(la)"
    iout_tag = "i(lout)"

if phase_bias:
    cd.change_param("phi1_mag", np.pi)
    cd.change_param("l3_mag", 0)
    cd.change_param("k2_mag", 0)
    cd.change_param("l4_mag", 0)
    cd.change_param("idc_mag", 0)

if no_bodefano:
    cd.change_param("l1_mag", 0)
    cd.change_param("k1_mag", 0)
    # cd.change_param("l2_mag", 0)

###########################################################################
results_dict = {f"{loop_param}": [], "phase_avg": [], "ijj_mag": [], "vjj_mag": [], 
            "iin_mag": [], "iloop_mag": [], "iout_mag": [], 
            "vinput_mag": [], "vin_mag": [], "vbase_mag": [], "vcancel_mag": [], "vout_mag": [], "voutput_mag": [],
            "leff_input": [], "energy_in": [], "energy_out": [], "energy_out_prev": [],
            "ifreq_fft": [], "iin_fft": [], "iout_fft": [], "iout_prev_fft": [], 
            "energy_in_fft": [], "energy_out_fft": [], "energy_out_prev_fft": []}
cur_loop_idx = 0
t1 = time.time()

for loop_val in loop_list:
    loop_len = len(loop_list)
    if loop_param is not None:
        cur_loop_idx += 1
        print(f"\n({cur_loop_idx}/{loop_len}): {loop_param} = {loop_val}")
        cd.change_param(loop_param, loop_val)
    if loop_param == "i1_freq":  # in this case, we need to change the time limits as well
        omega = loop_val * 2* np.pi
        freq = loop_val
        round_digits = np.floor(np.log10(freq)) + 1
        idc = 1.0339169242309647e-05  # 0  # to current bias to pi
        step_time = 10**(-(round_digits+3))
        stop_time = step_time * npts + idx_ringing * step_time
        cd.change_param("tran_step", step_time)
        cd.change_param("tran_stop", stop_time)  # 2e-10)
    cd.simulation_cycle(template, None, variation)

    time_array = cd.data["time"].to_numpy()[idx_ringing:]
    iin = cd.data[iin_tag].to_numpy()[idx_ringing:]
    iout = cd.data[iout_tag].to_numpy()[idx_ringing:]
    iloop = cd.data[iloop_tag].to_numpy()[idx_ringing:]
    vin = cd.data[vin_tag].to_numpy()[idx_ringing:]
    vinput = cd.data[vinput_tag].to_numpy()[idx_ringing:]
    vbase = cd.data[vbase_tag].to_numpy()[idx_ringing:]
    vcancel = cd.data[vcancel_tag].to_numpy()[idx_ringing:]
    vout = cd.data[vout_tag].to_numpy()[idx_ringing:]
    voutput = cd.data[voutput_tag].to_numpy()[idx_ringing:]
    ijj = cd.data[ijj_tag].to_numpy()[idx_ringing:]
    vjj = cd.data[vjj_tag].to_numpy()[idx_ringing:]
    phase = cd.data[phase_tag].to_numpy()[idx_ringing:]
    print(np.average(phase))

    if not no_bodefano and draw_plots:
        plt.plot(time_array, phase)
        plt.xlabel("Time (s)")
        plt.ylabel("Josephson junction phase")
        plt.grid()
        plt.tight_layout()
        plt.savefig("ic_results/ic_phase.png")
        plt.show()

    if "voltage_source" in variation and draw_plots:
        plt.plot(time_array, vin, label="Input voltage")
        plt.plot(time_array, vbase, label="Inductance to cancel")
        plt.plot(time_array, vcancel, label="Inductance canceling device")
        # plt.plot(time_array, vout, label="Output voltage")
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage (V)")
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig("ic_results/ic_voltage.png")
        plt.show()

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

    # generic input/output current graph
    if draw_plots:
        plt.plot(time_array, iin - offset(iin), zorder=2, label="Input (centered)")  #  - offset(iin)
        plt.plot(time_array, iout - offset(iout), zorder=1, label="Output (centered)")  #  - offset(iout)
        plt.legend()
        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")
        plt.tight_layout()
        plt.grid()
        plt.savefig(f"ic_results/ic_current.png")
        plt.show()
    
    if get_fft_data:
        iin_freq, iin_current = do_fft(time_array, iin)
        try: comp_idx = np.where(np.round(iin_freq) == freq)[0][0]
        except IndexError: comp_idx = np.where(iin_current == np.max(iin_current))[0][0]
        iout_freq, iout_current = do_fft(time_array, iout)
        if draw_plots:
            start_idx = np.min([1, comp_idx-5])
            plt.plot(iin_freq[start_idx:comp_idx+5], iin_current[start_idx:comp_idx+5], label="Input current")
            plt.plot(iout_freq[start_idx:comp_idx+5], iout_current[start_idx:comp_idx+5], label="Output current")
            plt.xlabel("Frequency (Hz)")
            plt.ylabel("FFT strength")
            plt.yscale("log")
            plt.legend()
            plt.savefig("ic_results/ic_fft.png")
            plt.show()
        iloop_freq, iloop_current = do_fft(time_array, iloop)
        
        ifreq_fft = iin_freq[comp_idx]
        iin_fft = iin_current[comp_idx]
        iout_fft = iout_current[comp_idx]
        iout_prev_fft = iloop_current[comp_idx]
        print(f"Max FFT freq at: {ifreq_fft}")

    # energy calcs  # currently need to turn on fft, maybe fix later
    leff = lin  # amplitude(vbase)/amplitude(iloop)/omega
    if "couple_out" in variation:
        vinput = cd.data[vinput_tag].to_numpy()[idx_ringing:]
        leff = amplitude(vinput)/amplitude(iin)/omega
    
    if "rout" in variation:  # variation ends with resistive output
        leff = cd.params["r1_mag"]  # input also has to be resistive to match dimensions
        loutput = cd.params["r1_mag"]  # i think this is correct, the 1/2 factor is due to ac amplitude
    
    energy_in = 1/2 * leff * amplitude(iin)**2
    energy_out = 1/2 * loutput * amplitude(iout)**2
    energy_in_fft = 1/2 * leff * iin_fft**2
    energy_out_fft = 1/2 * loutput * iout_fft**2
    print(f"INPUT INDUCTANCE COMPARISON - ACTUAL INPUT: {lin}, EFFECTIVE INPUT: {leff}\n\n")
    print(f"CURRENT COMP - INPUT: {amplitude(iin)}, OUTPUT: {amplitude(iout)}, \
            \nTRANSFER: {amplitude(iout)/amplitude(iin)}, TRANSFER BY FFT: {iout_fft / iin_fft}\n\n")
    print(f"ENERGY COMP - INPUT: {energy_in}, OUTPUT: {energy_out}, \
            \nTRANSFER: {energy_out/energy_in}, TRANSFER BY FFT: {energy_out_fft/energy_in_fft}\n\n")
    if cd.params["la_mag"] > 0:
        lout_prev_stage = cd.params["la_mag"]
    else:
        lout_prev_stage = cd.params["l0_mag"]
        
    energy_out_prev = 1/2 * lout_prev_stage * amplitude(iloop)**2
    energy_out_prev_fft = 1/2 * lout_prev_stage * iout_prev_fft**2
    print(f"(IF REQUIRED) ENERGY COMP @ PREV OUTPUT - INPUT: {energy_in}, OUTPUT: {energy_out_prev}, \
            \nTRANSFER: {energy_out_prev/energy_in}, TRANSFER BY FFT: {energy_out_prev_fft/energy_in_fft}\n\n")
    
    if loop_param is None: results_dict = None  # there's probably a better way to do this
    else: 
        results_dict[loop_param].append(loop_val)
        results_dict["phase_avg"].append(np.average(phase))
        results_dict["ijj_mag"].append(amplitude(ijj))
        results_dict["vjj_mag"].append(amplitude(vjj))
        results_dict["iin_mag"].append(amplitude(iin))
        results_dict["iloop_mag"].append(amplitude(iloop))
        results_dict["iout_mag"].append(amplitude(iout))
        results_dict["vinput_mag"].append(amplitude(vinput))
        results_dict["vin_mag"].append(amplitude(vin))
        results_dict["vbase_mag"].append(amplitude(vbase))
        results_dict["vcancel_mag"].append(amplitude(vcancel))
        results_dict["vout_mag"].append(amplitude(vout))
        results_dict["voutput_mag"].append(amplitude(voutput))
        results_dict["leff_input"].append(leff)
        results_dict["energy_in"].append(energy_in)
        results_dict["energy_out"].append(energy_out)
        results_dict["energy_out_prev"].append(energy_out_prev)
        results_dict["ifreq_fft"].append(ifreq_fft)
        results_dict["iin_fft"].append(iin_fft)
        results_dict["iout_fft"].append(iout_fft)
        results_dict["iout_prev_fft"].append(iout_prev_fft)
        results_dict["energy_in_fft"].append(energy_in_fft)
        results_dict["energy_out_fft"].append(energy_out_fft)
        results_dict["energy_out_prev_fft"].append(energy_out_prev_fft)

t2 = time.time()
print(f"Full loop took {(t2 - t1)/60} minutes.")

if save_results:
    if results_dict is None:
        results_dict = cd.data
        mode += "_onepoint"
    for key in results_dict.keys():  # remove empty
            if len(results_dict[key]) == 0: del results_dict[key]
    print("Saving Results.")
    data_pd = pd.DataFrame(results_dict)  # if loops
    data_pd.to_csv(f"ic_results/IndCancel_{variation}{mode}_data.csv", header=True, index=False)
    
if get_fft_data and loop_param is not None:
    # vna_like look at fft data
    vna_freq = np.array(results_dict["ifreq_fft"])
    vna_in = np.array(results_dict["energy_in_fft"])
    vna_out = np.array(results_dict["energy_out_fft"])
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Frequency (Hz)')
    ax1.set_ylabel('Input Energy (AU)', color="blue")
    ax1.plot(vna_freq, vna_in, color="blue")
    ax1.tick_params(axis='y', labelcolor="blue")
    ax2 = ax1.twinx()
    ax2.set_ylabel('Output-Input Energy Ratio (dB)', color="orange")
    ax2.plot(vna_freq, power2dB(vna_out/vna_in), label="Output/Input", color="orange")
    ax2.tick_params(axis='y', labelcolor="orange")
    '''
    if "resonant" in variation:
        vna_fit = lorentz_model.fit(vna_out/vna_in, x=vna_freq[0], par0=np.max(vna_out/vna_in), par1=vna_freq[len(vna_freq)//2], par2=(vna_freq[1]-vna_freq[0]))
        vna_fit_freq = np.linspace(loop_list[0], loop_list[-1], 1001)
        vna_fit_dB = power2dB(lorentzian(vna_fit_freq, np.max(vna_out/vna_in)*1.8, 4.565e7, (vna_freq[1]-vna_freq[0])/1.2))
        # np.max(vna_out/vna_in)*1.8, (vna_freq[6]+vna_freq[7]*1.05)/2.05, (vna_freq[1]-vna_freq[0])/1.2))
        # vna_fit.params['par0'].value, vna_fit.params['par1'].value, vna_fit.params['par2'].value))
        # np.max(vna_out/vna_in), (vna_freq[6]+vna_freq[7])/2, (vna_freq[1]-vna_freq[0])))
        print(vna_fit.fit_report())
        ax2.plot(vna_fit_freq, vna_fit_dB, "r--", label="Lorentzian Fit")
    '''
    ax2.legend()
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.grid()
    plt.savefig(f"ic_results/ic_fft_sweep_{variation}{mode}.png")
    plt.show()
    # plt.cla()
    
'''
# plt.plot(loop_list, power2dB(vna_out/vna_in))
# np.max(vna_out/vna_in)*1.8, (vna_freq[6]+vna_freq[7]*1.05)/2.05, (vna_freq[1]-vna_freq[0])/1.2)

# with hi res data

df1 = pd.read_csv("ic_results/IndCancel_voltage_source_resonant_nobodefano_data.csv")
df2 = pd.read_csv("ic_results/IndCancel_voltage_source_resonant_nobodefano_data_hires.csv")
freq1 = df1["i1_freq"].to_numpy()
vna1 = df1["energy_out_fft"].to_numpy() / df1["energy_in_fft"].to_numpy()
freq2 = df2["i1_freq"].to_numpy()
vna2 = df2["energy_out_fft"].to_numpy() / df2["energy_in_fft"].to_numpy()
plt.plot(freq1, power2dB(vna1), "b.", label="Resonant")
plt.plot(freq2, power2dB(vna2), "b.")
e1 = df1["energy_in_fft"].to_numpy()
e2 = df2["energy_in_fft"].to_numpy()  

df1 = pd.read_csv("ic_results/IndCancel_voltage_source_output_data.csv")
df2 = pd.read_csv("ic_results/IndCancel_voltage_source_output_nobodefano_data.csv")
freq1_ = df1["i1_freq"].to_numpy()
vna1 = df1["energy_out_fft"].to_numpy() / df1["energy_in_fft"].to_numpy()
freq2_ = df2["i1_freq"].to_numpy()
vna2 = df2["energy_out_fft"].to_numpy() / df2["energy_in_fft"].to_numpy()

# df3 = pd.read_csv("ic_results/IndCancel_couple_out_data.csv")
# freq3_ = df3["i1_freq"].to_numpy()
# vna3 = df3["energy_out_fft"].to_numpy() / df3["energy_in_fft"].to_numpy()
# e3_ = df3["energy_in_fft"].to_numpy()  

plt.plot(freq1_, power2dB(vna1), "g", label="Bode-Fano")
plt.plot(freq2_, power2dB(vna2), "r--", label="No Bode-Fano")

plt.plot(freq3_, power2dB(vna3), "c--", label="Bode-Fano, manual adjustment")

e1_ = df1["energy_in_fft"].to_numpy()
e2_ = df2["energy_in_fft"].to_numpy()  

plt.xlabel("Frequency (Hz)")
# plt.xlim(1e7, 1e8)
plt.ylabel("Output/Input (dB)")
# plt.yscale("log")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("ic_results/vna_comp_vsource2.png")
plt.show()

plt.plot(freq1, e1, "b.", label="Resonant")
plt.plot(freq2, e2, "b.")
plt.plot(freq1_, e1_, "g", label="Bode-Fano")
plt.plot(freq2_, e2_, "r--", label="No Bode-Fano")

# plt.plot(freq3_, e3_, "c--", label="Bode-Fano, manual adjustment")

plt.xlabel("Frequency (Hz)")
# plt.xlim(1e7, 1e8)
# plt.ylim(top=7e-10)
plt.ylabel("Power (Linear AU)")
plt.yscale("log")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("ic_results/e_comp_vsource2.png")
plt.show()

# customizable

df1 = pd.read_csv("ic_results/IndCancel_voltage_source_output_data.csv")
df2 = pd.read_csv("ic_results/IndCancel_voltage_source_output_nobodefano_data.csv")

df1 = pd.read_csv("ic_results/IndCancel_couple_out_data.csv")
df2 = pd.read_csv("ic_results/IndCancel_couple_out_nobodefano_data.csv")
freq1_ = df1["i1_freq"].to_numpy()
vna1 = df1["energy_out_fft"].to_numpy() / df1["energy_in_fft"].to_numpy()
freq2_ = df2["i1_freq"].to_numpy()
vna2 = df2["energy_out_fft"].to_numpy() / df2["energy_in_fft"].to_numpy()
e1_ = df1["energy_in_fft"].to_numpy()
e2_ = df2["energy_in_fft"].to_numpy()  

plt.plot(freq1_, power2dB(vna1), "g", label="Bode-Fano")
plt.plot(freq2_, power2dB(vna2), "r--", label="No Bode-Fano")
plt.xlabel("Frequency (Hz)")
# plt.xlim(1e7, 1e8)
plt.ylabel("Output/Input (dB)")
# plt.yscale("log")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("ic_results/vna_comp_coupleout.png")
plt.show()

plt.plot(freq1_, e1_, "g", label="Bode-Fano")
plt.plot(freq2_, e2_, "r--", label="No Bode-Fano")
plt.xlabel("Frequency (Hz)")
# plt.xlim(1e7, 1e8)
plt.ylabel("Power (Linear AU)")
plt.yscale("log")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("ic_results/e_comp_coupleout.png")
plt.show()

fig, ax1 = plt.subplots()

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Input Current (A)', color="blue")
ax1.plot(time_array, iin, color="blue")
ax1.tick_params(axis='y', labelcolor="blue")

ax2 = ax1.twinx()
ax2.set_ylabel('Output Current (A)', color="orange")
ax2.plot(time_array, iout, color="orange")
# ax2.plot(time_array, ics*np.sin(phase), "r--")
ax2.tick_params(axis='y', labelcolor="orange")

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.grid()
plt.savefig(f"ic_results/ic_current_comp.png")
plt.show()
plt.cla()

'''
