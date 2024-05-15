from circuit_reader import *
from circuit_calcs import *

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
