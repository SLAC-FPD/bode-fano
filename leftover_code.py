# includes (ctrl+f for to go to section):
## fit function models
## fft stuff
## parameter space plots

## fit function models
# definition of functions/fitting models using lmfit
def exp_decay(t, a, b, c):
    return a*np.exp(-t/b)+c

def sine_wave(t, a, freq, phase):
    return a*np.sin(2*np.pi*freq*t + phase)

def resistance_phase(phi, r, alpha):
    return r/(1+alpha*np.cos(phi))

def lorentzian(x, par0, par1, par2):  # linear case
    # par0: amplitude, par1: res freq, par2: bw
    return par0 / (1. + pow(((x - par1) / (par2 / 2.)), 2.))

exp_decay_model = Model(exp_decay)
sine_model = Model(sine_wave)
rp_model = Model(resistance_phase)
lorentz_model = Model(lorentzian)

## fft stuff
def do_fft(time, data, round_stepsize=True, round_duration=True):
    stepsize = time[1] - time[0]  # higher maxfreq when data is finer
    print(stepsize)
    if round_stepsize:  # let's find the first zero to truncate, not 100% guaranteed
        stepsize = round_value(stepsize)
        print(f"FFT: Rounded stepsize to {stepsize}.")
    duration = time[-1] - time[0]  # better rbw when data is longer
    print(duration)
    if round_duration:  # let's find the first zero to truncate, not 100% guaranteed
        duration = round_value(duration)
        print(f"FFT: Rounded duration to {duration}.")
    # rbw = 1/duration
    # sample_rate = 1/stepsize
    npts = round(duration/stepsize)+1  # must be integer
    # print(stepsize, duration, npts)
    # xf = np.linspace(0.0, npts*duration, npts, endpoint=False)
    # fft will remove last point in fft because of how fftfreq works
    if npts == len(data):
        xf = fftfreq(npts - 1, stepsize)  # [:-1]
    else: xf = fftfreq(len(data) - 1, stepsize)
    yf = fft(data[:-1])  # [:npts//2]
    return xf, yf

def round_value(value):  #, nine_count=3):  # don't do for > 1 yet
    decimals_to_round = 0
    value_mult = value
    # nine_counts = 0  # rounds values like 9.99995
    first_nonzero = False
    finished_rounding = False
    while not finished_rounding:
        # current_digit_value = int(value_mult)%10
        current_digit_value = round(value_mult)%10
        #print("NextDigit", current_digit_value, value_mult, decimals_to_round)
        #if nine_counts > nine_count:
        #    print("TooManyNines")
        #    decimals_to_round -= nine_count
        #    finished_rounding = True
        decimals_to_round += 1
        value_mult = value_mult*10
        if current_digit_value != 0:
            first_nonzero = True
        #    print("FirstNonzero")
        #if current_digit_value == 9:
        #    print("AddNineCount")
        #    nine_counts += 1
        if first_nonzero and current_digit_value == 0:
        #    print("FinishedRounding")
            finished_rounding = True
        #input("...")
    return round(value, decimals_to_round)

def get_zeff(t, v_l, i_mag, f, timestep, skip_idx=0, no_r=False):
    if skip_idx > 0:
        t = t[skip_idx:]
        v_l = v_l[skip_idx:]
    omega = 2*np.pi*f
    z_mag = (np.max(v_l) - np.min(v_l))/2/i_mag
    if no_r: leff, radd = z_mag / omega, 0
    else:
        idx_max = np.where(v_l == np.max(v_l))[0][0] + skip_idx
        period_idx = 1/f/timestep
        theta = (idx_max / period_idx - 1/4) * 2*np.pi
        ratio = np.tan(theta)
        radd = z_mag / np.sqrt(1 + ratio**2)
        leff = - ratio * radd / omega  # for not no_r cases only
    return(leff, radd)

## parameter space plots
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from CircuitReader import *

lt = 1.6e-7
omega = 30.e6
r = 1
k = 0.09125
# cur_l = 0.2
cur_l_list = np.linspace(-8.e-10, 8.e-10, 801)
# cur_l_list = np.linspace(0, 80, 81)
# l_list = np.linspace(0, 1, 11)

hist_2d = []
l_list = []

pts = np.array([-3.0e-10, 3.2e-10])  # choose specific point

for i in range(1000):
    # l = 2*i
    l = 1e-12 * (i+1)
    l_list.append(l)
    # leff_list = np.linspace(0+i, 10+i, 11)
    # leff_list = cur_l_list + l
    leff_list = calc_leff(lt, omega, r, l, cur_l_list, k)
    hist_2d.append(leff_list)

plt.imshow(hist_2d, extent=(cur_l_list[0], cur_l_list[-1], l_list[-1], l_list[0]), interpolation='none', aspect="auto")  # , norm=colors.LogNorm())

# plt.title(f"Effective Inductance for {omega/1000} kHz, k = {k}, R = {r} Ohm")
plt.xlabel("Josephson Junction Inductance")
plt.ylabel("Paired Circuit Inductance")

cbar = plt.colorbar()
plt.clim(-5e-6, 5e-6)
cbar.set_label("Total Effective Inductance", rotation=90)

plt.scatter(pts[0], pts[1], marker=".", color="red", s=200) # scatter plot after color bar (incorrect values)
# plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig(f"Leff_parameter_space_k{k}_R{r}.png")
plt.show()

print(calc_leff(lt, omega, r, pts[1], pts[0], k))
