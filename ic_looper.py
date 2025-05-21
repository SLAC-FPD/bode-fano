from circuit_reader import *
from circuit_calcs import *

cd = CircuitData()

# to do list?
'''
- Send results to separate folder
- Make loopable
'''

# what do we want to do
mode = "runner"  # "looper"  # 
template = "ind_cancel"  # "kent"  # "floquet"  # "kent_equiv"  # "parallel"  # 
param_file = "params/ind_cancel_params_250430.txt"  # None  # 
variation = "couple_out_load"  # "couple_out"  # None  # 
no_bodefano = False  # True  # 
phase_bias = False

# initialize
cd.read_template(template, param_file, variation)
results_to_watch = ["v(101)"]
save_results = False

# set base parameters
ic = 1.09e-6  # 1.6e-7  # 
lj0 = calc_lj(ic)
lj = lj0  # calc_lj(ic)
omega = 1.e6  # 
freq = omega/(2*np.pi)
idx_ringing = 1000  # if needed to remove

l0 = 3e-10  # 5e-10  # 
kin = np.sqrt(1./3.)  # used when coupling in and out
kout = kin  # these are free-ish parameters
lin = l0/2/(1 - kin**2)
lout = l0/2/(1 - kout**2)

linput = 1e-8  # 0  # input inductance
loutput = 1e-8  # 0  # output inductance in readout

k = 0.2 # between coupled circuits
cin = linput/lin
cout = loutput/lout
c = kout**2/(cout + 1) - 2
c_ = 1

a1 = (c_ + 1)*(2 - kin**2/(cin + 1) - kout**2/(cout + 1))
a2 = c_*(2 - kin**2/(cin + 1) - kout**2/(cout + 1))
b = c_ - 1 + k**2

det_x0 = (k**2*a2 - a1 + b*c)**2 - 4*c*a1*(c_ - b)
x0 = ((k**2*a2 - a1 + b*c) - np.sqrt(det_x0) ) / (2*c*a1)
xj = (c*x0 + k**2 - 1)/(c*x0 - 1)
print(f"Lin/L1 = {x0}, LJ/L2 = {xj}")

l1 = lin / x0  # 3e-10  # 5e-10  # 
l2 = lj / xj  # lj * 1.0296298147012206  # lj/(1-k**2/2)  # 2.1178582977171555e-09  # 

# do a solver later with a lg, k, and ltot that you want
lg = l2 * (1 - k**2*l1/(l0 + l1))
ltot = l0 + l1 * (1 - k**2*l2/(l2-lj))

print(f"LJ: {lj}, LG: {lg}, LTOT: {ltot}, L0: {l0}")

# change parameters from base
# cd.change_param("cpic_mag", 0)
# cd.change_param("rtype", 0)

time_scale = 1e0
cd.change_param("tran_step", 1/(omega*1e4)/time_scale)
cd.change_param("tran_stop", 7.5/omega/time_scale)  # 2e-10)
cd.change_param("tran_start", 0e-10)
cd.change_param("i1_freq", freq)
cd.change_param("icrit_mag", ic)
cd.change_param("ics1_mag", ic)
cd.change_param("l0_mag", l0)
cd.change_param("l1_mag", l1)
cd.change_param("l2_mag", l2)
cd.change_param("k1_mag", k)
cd.change_param("filename", f"IndCancel_{variation}")

iout_tag = "i(lb)"
iin_tag = "i(la)"
iloop_tag = "i(lb)"
vin_tag = "v(1)-v(0)"
vbase_tag = "v(2)-v(1)"
vcancel_tag = "v(0)-v(2)"
phase_tag = "v(101)"

if variation == "couple_out" or variation == "couple_out_load":
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
    iout_tag = "i(ld)"  # "i(lout)"

if phase_bias:
    cd.change_param("phi1_mag", np.pi)
    cd.change_param("l3_mag", 0)
    cd.change_param("k2_mag", 0)
    cd.change_param("l4_mag", 0)
    cd.change_param("idc_mag", 0)

if no_bodefano:
    cd.change_param("l1_mag", 0)
    cd.change_param("k1_mag", 0)
    cd.change_param("l2_mag", 0)
    cd.change_param("filename", f"IndCancel_{variation}_nobodefano")

cd.simulation_cycle(template, None, variation)

time_array = cd.data["time"].to_numpy()[idx_ringing:]
iin = cd.data[iin_tag].to_numpy()[idx_ringing:]
iout = cd.data[iout_tag].to_numpy()[idx_ringing:]
iloop = cd.data[iloop_tag].to_numpy()[idx_ringing:]
voltage = cd.data[vin_tag].to_numpy()[idx_ringing:]
voltage_l0 = cd.data[vbase_tag].to_numpy()[idx_ringing:]
voltage_l1 = cd.data[vcancel_tag].to_numpy()[idx_ringing:]
phase = cd.data[phase_tag].to_numpy()[idx_ringing:]

def amplitude(data):
    return (np.max(data) - np.min(data))/2

def offset(data):
    return np.max(data) - (np.max(data) - np.min(data))/2

plt.plot(time_array, phase)
plt.xlabel("Time (s)")
plt.ylabel("Josephson junction phase")
plt.grid()
plt.tight_layout()
# plt.savefig("ic_phase.png")
plt.show()

plt.plot(time_array, iin - offset(iin), zorder=2, label="Input")  #  - offset(iin)
plt.plot(time_array, iout - offset(iout), zorder=1, label="Output")  #  - offset(iout)
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Output current (A)")
plt.tight_layout()
plt.grid()
# plt.savefig("ic_current.png")
plt.show()
print(f"IN: {amplitude(iin)}, OUT: {amplitude(iout)}, TRANSFERRED: {amplitude(iout)/amplitude(iin)}")
print(f"INPUT VOLTAGE: {amplitude(voltage)}, CANCEL VOLTAGE: {amplitude(voltage_l1)}, CANCEL INDUCTANCE: {amplitude(voltage_l1)/amplitude(iloop)/omega}")

phase_avg = np.average(phase)
lj_ = np.abs(calc_lj(ic, phase_avg))
ltot_ = l0 + l1 * (1 - k**2*l2/(l2-lj_))
print(f"LJ: {lj} --> {lj_}, LG: {lg}, LTOT: {ltot} --> {ltot_}, L0: {l0}")
