from circuit_reader import *
from circuit_calcs import *
from scipy.optimize import fsolve

def official_plots():
    plt.rcParams.update({'font.size': 18, "text.usetex": False, "font.family": "sans-serif", "figure.figsize": "9, 7.5"})
    
ltot_stable_phase = True
ltot_high = True
ltot_low = True
ltot_negative_tvi = True
ltot_positive_tvi = True
if any([ltot_stable_phase, ltot_high, ltot_low, ltot_negative_tvi, ltot_positive_tvi]):
    print("Using official plot settings.")
    official_plots()
    ic = 1.6e-7
    lj = calc_lj(ic)
    l1 = 1e-9  # is multiplicative constant
    # Loop phase equation (made to equal 2*pi)
    loop_eq = lambda phi : lg/lj*np.sin(np.pi - phi) + phi + np.pi
    phi_list = []
    lg_list_full = np.linspace(0e-9, 5e-9, 501)
    lg_list_full = np.round(lg_list_full * 1e9, 2) / 1e9

    for lg in lg_list_full:
        phi_initial_guess = np.pi/2
        phi_solution = fsolve(loop_eq, phi_initial_guess)[0]
        phi_solution = phi_solution % (2*np.pi)
        if phi_solution > np.pi: phi_solution = 2*np.pi - phi_solution
        # print(lg, phi_solution)
        phi_list.append(phi_solution)

    phi_list = np.array(phi_list)

if ltot_stable_phase:
    plt.plot(lg_list_full/lj, phi_list)
    plt.plot([0, 5e-9/lj], [np.pi/2, np.pi/2], "r--") 
    plt.xlabel("$L_{G}$/$L_{J}$")  # Geometric Inductance (H)")  # (JJ POV, H)")
    plt.ylabel("Stable Josephson junction phase (rad)")
    plt.ylim(0, 4.0)
    plt.locator_params(axis='y', nbins=4)
    # plt.grid()
    plt.tight_layout()
    plt.savefig("ltot_stable_phase.png")
    print("Plotted stable JJ phase a function of LG/LJ.")
    plt.show()

if ltot_low:
    lg_list = np.linspace(0e-9, 2e-9, 3)
    l2 = lambda k: lg / (1 - k**2)
    ltot = lambda k: l1 * (lg + lj_)*(1 - k**2) / ((lg + lj_) - lj_*k**2)
    ytot = lambda k: ((lg + lj_) - lj_*k**2) / (l1 * (lg + lj_)*(1 - k**2))

    k_list = []
    k_linspace = np.linspace(0, 1, 101)
    l2_list = []

    for idx in range(len(lg_list)):
        idx_full = np.where(lg_list_full == lg_list[idx])[0][0]
        lg = lg_list_full[idx_full]
        phase_ = phi_list[idx_full]
        lj_ = calc_lj(1.6e-7, phase_)
        k_initial_guess = 0.5
        k_solution = fsolve(ytot, k_initial_guess)[0]
        k_solution = np.abs(k_solution)
        if k_solution == 0.5: k_solution = np.nan
        # print(lg, lj_, k_solution)
        k_list.append(k_solution)
        l2_list.append(l2(k_solution))
        plt.plot(k_linspace, ltot(k_linspace), ".", label=f"$L_G$ = {round(lg*1e9)/1e9}")

    plt.plot([0, 1], [l1, l1], "k--")
    plt.plot([0, 1], [-l1, -l1], "k--")
    plt.ylim(-10*l1, 10*l1)
    # plt.grid()
    plt.legend()
    plt.xlabel("Coupling Constant")
    plt.ylabel("Total Inductance (H)")
    plt.tight_layout()
    plt.savefig("ltot_low_lg.png")
    print("Plotted total inductance for smaller L_G.")
    plt.show()
    
if ltot_high:
    lg_list = np.linspace(2.5e-9, 3.5e-9, 3)
    l2 = lambda k: lg / (1 - k**2)
    ltot = lambda k: l1 * (lg + lj_)*(1 - k**2) / ((lg + lj_) - lj_*k**2)
    ytot = lambda k: ((lg + lj_) - lj_*k**2) / (l1 * (lg + lj_)*(1 - k**2))

    k_list = []
    k_linspace = np.linspace(0, 1, 101)
    l2_list = []

    for idx in range(len(lg_list)):
        idx_full = np.where(lg_list_full == lg_list[idx])[0][0]
        print(idx_full)
        lg = lg_list_full[idx_full]
        phase_ = phi_list[idx_full]
        lj_ = calc_lj(1.6e-7, phase_)
        k_initial_guess = 0.5
        k_solution = fsolve(ytot, k_initial_guess)[0]
        k_solution = np.abs(k_solution)
        if k_solution == 0.5: k_solution = np.nan
        # print(lg, lj_, k_solution)
        k_list.append(k_solution)
        l2_list.append(l2(k_solution))
        plt.plot(k_linspace, ltot(k_linspace), ".", label=f"$L_G$ = {round(lg*1e9, 1)/1e9}")

    plt.plot([0, 1], [l1, l1], "k--")
    plt.plot([0, 1], [-l1, -l1], "k--")
    plt.ylim(-10*l1, 10*l1)
    # plt.grid()
    plt.legend()
    plt.xlabel("Coupling Constant")
    plt.ylabel("Total Inductance (H)")
    plt.tight_layout()
    plt.savefig("ltot_high_lg.png")
    print("Plotted total inductance for larger L_G.")
    plt.show()

if ltot_negative_tvi:
    data_pd = pd.read_csv("ind_cancel_specific_negative.csv")
    current = data_pd["i(lb)"].to_numpy()[1000:]
    voltage = data_pd["v(1)-v(0)"].to_numpy()[1000:]
    voltage_l0 = data_pd["v(2)-v(1)"].to_numpy()[1000:]
    voltage_l1 = data_pd["v(0)-v(2)"].to_numpy()[1000:]
    time_array2 = data_pd["time"][1000:]
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
    # plt.grid()
    plt.savefig(f"ltot_specific_negative.png")
    print("Plotted current-phase lag for inductance cancelled case.")
    plt.show()
    
if ltot_positive_tvi:
    data_pd = pd.read_csv("ind_cancel_specific_positive.csv")
    current = data_pd["i(lb)"].to_numpy()[1000:]
    voltage = data_pd["v(1)-v(0)"].to_numpy()[1000:]
    voltage_l0 = data_pd["v(2)-v(1)"].to_numpy()[1000:]
    voltage_l1 = data_pd["v(0)-v(2)"].to_numpy()[1000:]
    time_array2 = data_pd["time"][1000:]
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
    # plt.grid()
    plt.savefig(f"ltot_specific_positive.png")
    print("Plotted current-phase lag for NON-inductance cancelled case.")
    plt.show()
    
def window_size(is_smaller_screen=True):
    # set the window size defaults
    if is_smaller_screen:
        plt.rcParams["figure.figsize"] = [8,6]
        plt.rcParams.update({'font.size': 16})
    else:
        plt.rcParams["figure.figsize"] = [12,9]
        plt.rcParams.update({'font.size': 20})

def default_plt_settings(xlabel="Time (s)", ylabel="", legend=True, grid=True, tight_layout=True, xlims=None, ylims=None):
    # just since I used these all the time
    # these things could be done with some *args setup? not sure
    plt.xlabel(f"{xlabel}")
    plt.ylabel(f"{ylabel}")
    if xlims is not None:
        plt.xlim(xlims)
    if ylims is not None:
        plt.xlim(ylims)
    if legend: plt.legend()
    if grid: plt.grid()
    if tight_layout: plt.tight_layout()
        
def plot_cd_measurables(cd, input_params=None, show_plots=True, save_plots=False, remove_offsets=False, offset_idx=0, xlims=None, ylims=None):
    # cd: circuit data object   
    # default labels for axes
    labels = {"phases": "Phase (deg)", "currents": "Current (A)", "voltages": "Voltage (V)",
              "others": "A.U.", "single_phases": "Phase (deg)", "inductances": "Inductance (H)",
              "capacitances": "Capacitance (F)", "resistances": "Resistance (ohm)"}
    if cd is not None: 
        time_data = cd.data["time"]  # should always exist
        # case 1 (no specified input_params): plot all measurables stored in the CircuitData object
        if input_params is None:
            input_param_list = ["phases", "currents", "voltages"]
            for input_param in input_param_list:
                if len(cd.measurables[input_param]) > 0:  # is this check required?
                    for meas in cd.measurables[input_param]:
                        meas_val = cd.data[f"{meas}"]
                        # removing offsets in phases is not required
                        if remove_offsets and input_param != "phases": meas_val = remove_offset(meas_val)
                        plt.plot(time_data, meas_val, label=f"{meas}")
                default_plt_settings(ylabel=labels[input_param], xlims=xlims, ylims=ylims)
                if save_plots: plt.savefig(f"{cd.filename}_{input_param}")
                if show_plots: plt.show()
        # case 2 (specified input_params): plot all measurables stored in one of the measurable types
        elif input_params in ["phases", "currents", "voltages", "others", "single_phases"]:
            for meas in cd.measurables[input_params]:
                meas_val = cd.data[f"{meas}"]
                if remove_offsets and input_params != "phases": meas_val = remove_offset(meas_val)
                plt.plot(time_data, meas_val, label=f"{meas}")
            default_plt_settings(ylabel=labels[input_params], xlims=xlims, ylims=ylims)
            if save_plots: plt.savefig(f"{cd.filename}_{input_params}")
            if show_plots: plt.show()

def remove_offset(meas, idx=0):
    # removes the constant that creates an offset (such as a DC current) in the measurable (average it out)
    # if provided an idx, the offset calculation begins from that index
    # meas = meas.to_numpy()
    offset = np.average(meas[idx:])
    return meas - offset

def plot_two_axes_measurables():
    pass      
