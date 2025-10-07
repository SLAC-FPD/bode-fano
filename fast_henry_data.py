import numpy as np

class FastHenryData:
    def __init__(self, filename="Zc.mat", target_freq=10000):  # reads WRSpice output data
        self.reset()
        self.filename = filename
        self.target_freq = target_freq
        self.read_mat()
        # self.inductance_summary()
        
    def reset(self):  # reads fasthenry data
        self.inds = {}
        self.r_mat = None
        self.l_mat = None
        self.k_mat = None
            
    def read_mat(self):  # reads basic .mat file structure
        infile = open(self.filename, "r")
        read_impedance = False
        for line in infile:
            # print(line)
            line_spl = line.split()  # [i for i in line.split(" ") if i !=""]  # remove all blanks from splits
            if line_spl[0] == "Row":
                print("Reading names of inductors")
                self.inds[int(line_spl[1].strip(":"))] = [0, f"{line_spl[2]} {line_spl[3]} {line_spl[4]}"]
            elif line_spl[0] == "Impedance":
                n_inds = len(self.inds.keys())
                freq = float(line_spl[-4])
                if freq == self.target_freq:
                    print(f"Preparing matrix for freq = {freq} Hz")
                    self.r_mat = np.zeros([n_inds, n_inds])
                    self.l_mat = np.zeros([n_inds, n_inds])
                    self.k_mat = np.zeros([n_inds, n_inds])
                    y_cnt = 0  # rows
                    read_impedance = True
                else: read_impedance = False
            elif read_impedance:  # let's read the frequencies
                for x_cnt in range(n_inds):
                    self.r_mat[x_cnt][y_cnt] = float(line_spl[2*x_cnt]) # resistances
                    self.l_mat[x_cnt][y_cnt] = float(line_spl[2*x_cnt + 1][:-1]) / (2*np.pi*freq) # inductances
                y_cnt += 1
            else: pass
        for cnt in range(n_inds):
            self.inds[cnt + 1][0] = self.l_mat[cnt][cnt]
        infile.close()
    
    def inductance_summary(self):
        if self.l_mat is None: print("Impedance Matrix is Empty.")
        else:
            n_inds = len(self.inds.keys())
            for y_cnt in range(n_inds):
                for x_cnt in range(n_inds):
                    self.k_mat[x_cnt][y_cnt] = self.l_mat[x_cnt][y_cnt] / np.sqrt(self.l_mat[x_cnt][x_cnt] * self.l_mat[y_cnt][y_cnt])
                    if x_cnt == y_cnt: print(f"{x_cnt + 1}: {self.inds[x_cnt + 1][1]} is {self.inds[x_cnt + 1][0]} H")
                    elif x_cnt > y_cnt: print(f"k_{x_cnt+1}{y_cnt+1} is {self.k_mat[x_cnt][y_cnt]}, M_{x_cnt+1}{y_cnt+1} is {self.l_mat[x_cnt][y_cnt]} H")


