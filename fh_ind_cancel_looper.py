import numpy as np
import os
import fast_henry_data as fhd
import pandas as pd
from scipy.optimize import fsolve

lj0 = 3.019320903444526e-10
loop_eq = lambda phi : lg/lj0*np.sin(np.pi - phi) + phi + np.pi

lin_target = 3.183126393426777e-10
l0_target = 3.183126393426777e-10
l1_target = 3.29730458e-10
l2_target = 3.44893e-10
k0_target = 0.3162307230947436
k_target = 0.3

# fasthenry params
w = 1.e-6  # width of traces
h = 1e-7  # height of traces
sigma = 5.8e7
freq = 5e7

'''
x_in = 210e-6
y_in = 27e-6
y_m_in = 14e-6
y_l1 = 32e-6
y_m_jj = 1.2e-6
y_l2 = 13e-6
x_out = 80e-6
'''

x_in_list = [210e-6]  # np.linspace(200e-6, 220e-6, 21)  # 
y_in_list = [27e-6]  # np.linspace(25e-6, 35e-6, 11) # 
y_m_in_list = [14e-6]  # np.linspace(11e-6, 17e-6, 7) # 
y_l1_list = [32e-6]  # np.linspace(25e-6, 40e-6, 16) # 
y_m_jj_list = [1.2e-6]  # [1.2e-6, 2.e-6]  # np.linspace(1e-6, 3e-6, 3) # 
y_l2_list = [13e-6]  # np.linspace(10e-6, 15e-6, 6) # 
x_out_list = [80e-6]
results = {"index": [],
           "x_in": [], "y_in": [], "y_m_in": [], "y_l1": [], "y_m_jj": [], "y_l2": [], "x_out": [],
           "l_in": [], "l_main": [], "l_jj1": [], "l_jj2": [], 
           "m_12": [], "m_13": [], "m_14": [], "m_23": [], "m_24": [], "m_34": [],
           # "l1": [], "l2": [], "k12": [], "lg": [], "lj": [], "jj_phase": [], "l_bfe": []
           "score": []
           }
notable_indices = []

filename = "bfe_parallel_circuit"
f_ext = "_bfe_parallel"  # for zbuffile
index = 0

for x_in in x_in_list:
    for y_in in y_in_list:
        for y_m_in in y_m_in_list:
            for y_l1 in y_l1_list:
                for y_m_jj in y_m_jj_list:
                    for y_l2 in y_l2_list:
                        results["index"].append(index)
                        print(index)
                        x_out = 80e-6  # not relevant to the other parts
                        fasthenry_block = f'''
Nin0 x=0 y=0
Nin1 x=0 y={y_in}
Nin2 x={x_in} y={y_in}
Nin3 x={x_in} y=0

Na0 x=0 y={y_in + y_m_in}
Na1 x=0 y={2*y_in + y_m_in}
Na2 x=0 y={2*y_in + y_m_in + y_l1}
Na3 x={x_in} y={2*y_in + y_m_in + y_l1}
Na4 x={x_in} y={2*y_in + y_m_in}
Na5 x={x_in} y={2*y_in + y_m_in}
Na6 x={x_in} y={y_in + y_m_in}

Nb0 x=0 y={2*y_in + y_m_in + y_l1 + y_m_jj}
Nb1 x=0 y={2*y_in + y_m_in + y_l1 + y_m_jj + y_l2}
Nb2 x={x_in/2 - w/2} y={2*y_in + y_m_in + y_l1 + y_m_jj + y_l2}
Nb3 x={x_in/2 + w/2} y={2*y_in + y_m_in + y_l1 + y_m_jj + y_l2}
Nb4 x={x_in} y={2*y_in + y_m_in + y_l1 + y_m_jj + y_l2}
Nb5 x={x_in} y={2*y_in + y_m_in + y_l1 + y_m_jj}

Nout0 x={-x_out} y={2*y_in + y_m_in}
Nout1 x=0 y={2*y_in + y_m_in}
Nout2 x={x_in} y={2*y_in + y_m_in}
Nout3 x={x_in + x_out} y={2*y_in + y_m_in}

Ein0 Nin0 Nin1 w={w} h={h}
Ein1 Nin1 Nin2 w={w} h={h}
Ein2 Nin2 Nin3 w={w} h={h}
Ea0 Na0 Na1 w={w} h={h}
Ea1 Na1 Na2 w={w} h={h}
Ea2 Na2 Na3 w={w} h={h}
Ea3 Na3 Na4 w={w} h={h}
Ea4 Na5 Na6 w={w} h={h}
Ea5 Na6 Na0 w={w} h={h}
Eb0 Nb0 Nb1 w={w} h={h}
Eb1 Nb1 Nb2 w={w} h={h}
Eb2 Nb3 Nb4 w={w} h={h}
Eb3 Nb4 Nb5 w={w} h={h}
Eb4 Nb5 Nb0 w={w} h={h}
Eout0 Nout0 Nout1 w={w} h={h}
Eout1 Nout2 Nout3 w={w} h={h}
'''
                        f = open(f"fh_results/{filename}.inp", "w")
                        f.write(f"**non-foster circuit**\n.Default z=0 sigma={sigma}\n\n")

                        f.write(fasthenry_block)
                        f.write(f".external Nin0 Nin3\n.external Na1 Na5\n.external Na1 Na4\n.external Nb2 Nb3\n.freq fmin={freq} fmax={freq*10} ndec=1\n\n.end")
                        f.close()

                        os.system(f"fasthenry fh_results/{filename}.inp")
                        # os.system(f"fasthenry fh_results/{filename}.inp -f simple -S {f_ext}") # -b to end wrspice after run
                        # os.system(f"zbuf zbuffile{f_ext}") # -b to end wrspice after run
                        
                        fh_data = fhd.FastHenryData(target_freq=freq)
                        # print(fh_data.inds)
                        l_mat = fh_data.l_mat
                        k_mat = fh_data.k_mat
                        results["x_in"].append(x_in)
                        results["y_in"].append(y_in)
                        results["y_m_in"].append(y_m_in)
                        results["y_l1"].append(y_l1)
                        results["y_m_jj"].append(y_m_jj)
                        results["y_l2"].append(y_l2)
                        results["x_out"].append(x_out)
                        results["l_in"].append(l_mat[0][0])
                        results["l_main"].append(l_mat[1][1])
                        results["l_jj1"].append(l_mat[2][2])
                        results["l_jj2"].append(l_mat[3][3])
                        results["m_12"].append(l_mat[0][1])
                        results["m_13"].append(l_mat[0][2])
                        results["m_14"].append(l_mat[0][3])
                        results["m_23"].append(l_mat[1][2])
                        results["m_24"].append(l_mat[1][3])
                        results["m_34"].append(l_mat[2][3])
                        score = np.sqrt(((lin_target-l_mat[0][0])/lin_target)**2+((l0_target-l_mat[1][1])/l0_target)**2+
                                        ((l1_target-l_mat[2][2])/l1_target)**2+((l2_target-l_mat[3][3])/l2_target)**2+
                                        ((k0_target-k_mat[0][1])/k0_target)**2+((k_target-k_mat[2][3])/k_target)**2)
                        results["score"].append(score)
                        if score < 0.02:
                            notable_indices.append(index)
                            print(f"Closeness: {score}")
                        index += 1
                        '''
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
                        '''

print(notable_indices)
df = pd.DataFrame(results)
df.to_csv(f"fh_results/fast_henry_results{f_ext}.csv", index=True, header=True)  # _closeup

'''
import matplotlib.pyplot as plt
plt.plot(results["l_bfe"], ".")  # "lg" "lj" "jj_phase" "l_bfe"
bin_list = np.linspace(-8e-10, -2.8e-10, 27)
plt.hist(results["l_bfe"], bins=bin_list)
plt.xlabel("LCANCEL Inductance (H)")
plt.ylabel("Count")
plt.show()

'''
