import numpy as np
import os
import fast_henry_data as fhd
import pandas as pd
from scipy.optimize import fsolve

lj0 = 3.019320903444526e-10
loop_eq = lambda phi : lg/lj0*np.sin(np.pi - phi) + phi + np.pi

# fasthenry params
w = 1.e-6  # width of traces
h = 1e-7  # height of traces
sigma = 5.8e7

loop1_x_list = [179e-6]  # np.linspace(178.5e-6, 179.5e-6, 11) # np.linspace(100e-6, 200e-6, 11) # 
loop1_y_list = np.linspace(100e-6, 200e-6, 101)  # np.linspace(134.5e-6, 135.5e-6, 11) # [136e-6]  # 
loop2_x_list = np.linspace(20e-6, 60e-6, 41)  # np.linspace(32.5e-6, 33.5e-6, 11) # [32e-6]  # 
loop12_list = np.linspace(1e-6, 2e-6, 2)  # np.linspace(0.95e-6, 1.05e-6, 11) # [1e-6]  # 
loop13_list = [2e-6]  # np.linspace(1.95e-6, 2.05e-6, 11) # 
results = {"loop1_x": [], "loop1_y": [], "loop2_x": [], "loop12": [], "loop13": [], "l_main": [], "l_jj": [], 
           "m_12": [], "l_in": [], "l_out": [], "m_13": [], "m_14": [], "m_23": [], "m_24": [], "m_34": [],
           "l1": [], "l2": [], "k12": [], "lg": [], "lj": [], "jj_phase": [], "l_bfe": []}

filename = "ind_cancel_circuit"
f_ext = "_ind_cancel"  # for zbuffile
index = 0

for loop1_x in loop1_x_list:
    for loop1_y in loop1_y_list:
        for loop2_x in loop2_x_list:
            for loop12 in loop12_list:
                for loop13 in loop13_list:
                    loop2_y = loop1_y
                    loop14 = loop13
                    loop3_x = loop1_x  # length of input line
                    loop4_x = loop1_x  # length of output line
                    fasthenry_block = f'''
Na0 x=0 y=0
Na1 x=0 y={loop1_y}
Na2 x={-loop1_x} y={loop1_y}
Na3 x={-loop1_x} y=0
Na4 x={-w/2} y=0

Nb0 x={loop12} y=0
Nb1 x={loop12} y={loop2_y}
Nb2 x={loop2_x + loop12} y={loop2_y}
Nb3 x={loop2_x + loop12} y=0
Nb4 x={loop12 + w/2} y=0

Nin0 x=0 y={-loop13}
Nin1 x={-loop3_x} y={-loop13}

Nout0 x=0 y={loop1_y + loop13}
Nout1 x={-loop4_x} y={loop1_y + loop13}

Ea0 Na0 Na1 w={w} h={h}
Ea1 Na1 Na2 w={w} h={h}
Ea2 Na2 Na3 w={w} h={h}
Ea3 Na3 Na4 w={w} h={h}
Eb0 Nb0 Nb1 w={w} h={h}
Eb1 Nb1 Nb2 w={w} h={h}
Eb2 Nb2 Nb3 w={w} h={h}
Eb3 Nb3 Nb4 w={w} h={h}
Ein0 Nin0 Nin1 w={w} h={h}
Eout0 Nout0 Nout1 w={w} h={h}
'''
                    f = open(f"fh_results/{filename}.inp", "w")
                    f.write(f"**fasthenry example**\n.Default z=0 sigma={sigma}\n\n")

                    f.write(fasthenry_block)
                    f.write(".external Na0 Na4\n.external Nb0 Nb4\n.external Nin0 Nin1\n.external Nout0 Nout1\n.freq fmin=1e4 fmax=1e5 ndec=1\n\n.end")
                    f.close()

                    os.system(f"fasthenry fh_results/{filename}.inp")
                    # os.system(f"fasthenry fh_results/{filename}.inp -f simple -S {f_ext}") # -b to end wrspice after run
                    # os.system(f"zbuf zbuffile{f_ext}") # -b to end wrspice after run

                    fh_data = fhd.FastHenryData()
                    # print(fh_data.inds)
                    l_mat = fh_data.l_mat
                    results["loop1_x"].append(loop1_x)
                    results["loop1_y"].append(loop1_y)
                    results["loop2_x"].append(loop2_x)
                    results["loop12"].append(loop12)
                    results["loop13"].append(loop13)
                    results["l_main"].append(l_mat[0][0])
                    results["l_jj"].append(l_mat[1][1])
                    results["m_12"].append(l_mat[0][1])
                    results["l_in"].append(l_mat[2][2])
                    results["l_out"].append(l_mat[3][3])
                    results["m_13"].append(l_mat[0][2])
                    results["m_14"].append(l_mat[0][3])
                    results["m_23"].append(l_mat[1][2])
                    results["m_24"].append(l_mat[1][3])
                    results["m_34"].append(l_mat[2][3])
                    l1 = l_mat[0][0] - np.abs(l_mat[2][2]) - np.abs(l_mat[3][3])
                    l2 = l_mat[1][1]
                    k12 = l_mat[0][1] / np.sqrt(l1*l2)
                    lin = l_mat[2][2]
                    lout = l_mat[3][3]
                    k13 = np.abs(l_mat[0][2]/l_mat[2][2])
                    k14 = np.abs(l_mat[0][3]/l_mat[3][3])
                    lg = l2*(1 - k12**2*l1/(l1 + lin*(1 - k13**2) + lout*(1 - k14**2)))
                    phi_initial_guess = np.pi/2
                    phi_solution = fsolve(loop_eq, phi_initial_guess)[0]
                    phi_solution = phi_solution % (2*np.pi)
                    if phi_solution > np.pi: phi_solution = 2*np.pi - phi_solution
                    lj = lj0 / np.cos(phi_solution)
                    l_bfe = l1*(1 - k12**2*l2/(l2+lj))
                    results["l1"].append(l1)
                    results["l2"].append(l2)
                    results["k12"].append(k12)
                    results["lg"].append(lg)
                    results["lj"].append(lj)
                    results["jj_phase"].append(phi_solution)
                    results["l_bfe"].append(l_bfe)

df = pd.DataFrame(results)
df.to_csv("fh_results/fast_henry_results.csv", index=True, header=True)  # _closeup

'''
import matplotlib.pyplot as plt
plt.plot(results["l_bfe"], ".")  # "lg" "lj" "jj_phase" "l_bfe"
bin_list = np.linspace(-8e-10, -2.8e-10, 27)
plt.hist(results["l_bfe"], bins=bin_list)
plt.xlabel("LCANCEL Inductance (H)")
plt.ylabel("Count")
plt.show()

'''
