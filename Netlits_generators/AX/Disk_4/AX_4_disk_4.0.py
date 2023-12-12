"""
Generating a netlist for 1 disk

"""

import numpy as np

d = 2.3622      # diameter of conductor (in mm)
t = 1.143       # thickness of insulator (in mm)
D_core = 110     # air core diameter (in mm)
N_d = 10        # number of disks
N_t = 6         # number of turns in each disk
e_r = 4         # relative permittivity of silicon rubber insulator
e_0 = 8.854e-12 # permittivity of free space
r_dc_m = 0.00574# DC resistance (Ohm/m)
r_dc_1000ft = 1.75  # DC resistance (Ohm/1000feet)
ft_to_m = 0.3048# meter to feet conversion

d_sep = 4       # number of separated disks

D_wire = d + 2*t                        # diameter of wire
d_dist = 4.0                              # distance between centers of conductors for axial displacement

R_turn = np.zeros((1,N_t))              # initializing turn radius vector (in mm)
R_mean = (N_t/2)*(2*t + d) + D_core/2   # mean radius of disk (in mm)

# Calculation of radii of turns
for i in range(N_t):
    R_turn[:,N_t-1-i] = (D_core/2) + (2*i + 1)*(t + d/2)
    
# %% Inductance calculation
 
# Constants for mutual inductance calculation
# Table of constant f for k_sq > 0.1
k_sq_L = np.round(np.arange(0.01, 1.01, 0.01), 2)
f_L = np.array([0.021474, .017315, .014937, .013284, .012026, .011017, .010179, .009464, .008843, .008297, .007810, .007371, .006974, .006611, .006278, .005970, .005685, .005420, .005173, .004941, .004723, .004518, .004325, .004142, .003969, .003805, .003649, .003500, .003359, .003224, .003095, .002971, .002853, .002740, .0026317, .0025276, .0024276, .0023315, .0022391, .0021502, .0020646, .0019821, .0029026, .0018259, .0017519, .0016805, .0016116, .00015451, .0014808, .0014186, .0013585, .0013004, .0012443, .0011900, .0011374, .0010865, .0010373, .0009897, .0009436, .0008990, .0008558, .0008141, .0007736, .0007345, .0006966, .0006600, .0006246, .0005903, .0005571, .0005251, .0004941, .0004642, .0004353, .0004074, .0003805, .0003545, .0003295, .0003054, .0002823, .00025998, .00023859, .00021806, .00019840, .00017959, .00016162, .00014450, .00012821, .00011276, .00009815, .00008438, .00007146, .00005940, .00004824, .00003798, .00002866, .00002035, .00001312, .00000708, .00000249, .0])
diff_L1 = f_L[1:] - f_L[:-1]
diff_L2 = diff_L1[1:] - diff_L1[:-1]
# Table of constant f for k_sq <= 0.1
log_k_sq_S = np.round(np.arange(-6.0, -0.9, 0.1), 1)
k_sq_S = 10**(log_k_sq_S)
f_S = np.array([.079093, .077647, .076200, .074753, .073306, .071860, .070413, .068966, .067520, .066073, .064626, .063180, .061733, .060287, .058840, .057394, .055947, .054500, .053055, .051609, .050163, .048717, .047272, .045827, .044382, .042938, .041494, .040051, .038608, .037167, .035727, .034288, .032851, .031416, .029984, .028554, .027128, .025707, .024291, .022881, .021478, .020084, .018700, .017329, .015972, .014632, .013311, .012013, .010742, .009502, .008297])
diff_S1 = f_S[1:] - f_S[:-1]
diff_S2 = diff_S1[1:] - diff_S1[:-1]

# function that determines the constant "f" for mutual inductance calculation
def f_calc(k_sq):
    if k_sq <= 1e-6:
        # estimating f
        f = 0.014468 * (np.log10(1/k_sq) - 0.53307)
    elif k_sq <= 0.01:
        # calculating log10 of k_sq
        log_k_sq = np.log10(k_sq)
        # flooring the number to the first decimal
        log_k_sq_floor = np.floor(log_k_sq * 10) / 10
        # finding the corresponding value of f
        f1 = f_S[log_k_sq_floor*10 == log_k_sq_S*10]
        # finding the first order difference
        diff1 = diff_S1[(log_k_sq_floor == log_k_sq_S)[:-1]]
        # finding the second-order difference
        diff2 = diff_S2[(log_k_sq_floor == log_k_sq_S)[:-2]]
        # calculating the value of "u"
        u = (log_k_sq - log_k_sq_floor) * 10
        # calculating u2 = (1-u)/2
        u2 = (1 - u) / 2
        # estimating f
        f = f1 + u * (diff1 - u2*diff2)
    elif k_sq < 1:
        # flooring the number to the second decimal
        k_sq_floor = np.floor(k_sq * 100) / 100
        # finding corresponding value of f
        f1 = f_L[k_sq_floor == k_sq_L]
        # finding the first-order difference
        diff1 = diff_L1[(k_sq_floor == k_sq_L)[:-1]]
        # finding the second-order difference
        diff2 = diff_L2[(k_sq_floor == k_sq_L)[:-2]]
        # calculating the value "u"
        u = (k_sq - k_sq_floor) * 100
        # calculating u2 = (1-u)/2
        u2 = (1 - u) / 2
        # estimating f
        f = f1 + u * (diff1 - u2*diff2)
    elif k_sq == 1:
        f = 0
    elif k_sq > 1:
        f = []
        print("Error in k_sq")
    return f        
    
    
# Self-induction calculation using Kirchoff equation [59] Grover
L = 4*np.pi*(R_turn*1e-3)*(np.log(8*R_turn/(d/2))-1.75)*1e-7
#L_micro = np.round(L * 1e6, 4)
L_micro = L * 1e6

# Mutual inductance calculation for filament (3 or more disks)
M = np.zeros((N_t*N_d,N_t*N_d))     # initializing mutual inductance matrix
y = []                              # vertical distance of center of coils from top coil
R = R_turn                          # initialize the radial distance vector
Rflip = R_turn                      # flipping radial dimensions
for i in range(N_d):
    y += [i * (d + 2*t)]*N_t   # distances between center of coils
y[d_sep*N_t:] = list(np.asarray(y[d_sep*N_t:]) + (d_dist - D_wire)) # adding disk separation 


for i in range(N_d-1):
    R = np.concatenate((R, np.flip(Rflip,axis=1)),axis=1)     # radial distances from central axis
    Rflip = np.flip(Rflip,axis=1)
    
    
    
for i in range(N_t*N_d):
    for j in range(i+1,N_t*N_d):
        # finding k^2 from Grover
        k_sq = ((y[i] - y[j])**2 + (R[:,i] - R[:,j])**2)/((y[i] - y[j])**2 + (R[:,i] + R[:,j])**2)
        # calculation of f
        f = f_calc(k_sq)
        # calculation of mutual inductance
        M[i,j] = f*(R[:,i]*R[:,j]*1e-2)**0.5*1e-6
#M_micro = np.round(M * 1e6, 4)      # converting to microhenries
M_micro = M * 1e6

M_trans = np.transpose(M_micro)     # transposing the inductances
M_tot = M_micro + M_trans           # adding matrices
index_flip = np.arange(N_t)         # index of self-inductance 
cur_disk = 1                        # current disk
cc = 0                              # counter of turns
for h in range(N_d*N_t): 
    M_tot[h,h] = L_micro[:,index_flip[h-N_t*(cur_disk-1)]]     # adding diagonal elements
    cc += 1                         # increment turn number
    if cc == N_t:
        cur_disk += 1               # increment disk counter
        cc = 0                      # nullify turn counter
        index_flip = np.flip(index_flip, axis=0)        # flip index vector
    

M_fin = M_tot * 1e-6               # from micro to true value 2 disks
K_rel = np.linalg.inv(M_fin)      # 2 disks reluctance

#M1_fin = M_tot[:N_t,:N_t] * 1e-6    # inductance matrix for 1 disk
#K1_rel = np.linalg.inv(M1_fin)      # reluctance matrix for 1 disk
Lnew = 1 / np.diag(K_rel) * 1e6      # new inductances in microHenries

alpha = np.zeros((N_t*N_d,N_t*N_d))
for i in range(N_t*N_d):
    alpha[i,:] = K_rel[i,:]/K_rel[i,i]

# %% Capacitance calculation
    
# Calculation of radie of middle of conductors (in mm)
R_mid = R_turn[:,:N_t-1] - (d/2) - t

# Calculation of length of middle of conductors lt (in mm)
l_mid = 2 * np.pi * R_mid

# Calculation of turn-to-turn capacitances between turns of same disk
theta = np.arccos(1 - np.log(D_wire/d)/e_r)     # crossing angle Theta*
Ctt = e_0 * (l_mid*1e-3) * (e_r * theta / np.log(D_wire/d) + 1/np.tan(theta/2) - 1/np.tan(np.pi/12))
Ctt_piko = Ctt * 1e12
Ctt_p = Ctt_piko

# Calculation of disk-to-disk capacitances between turns of adjacent disks
Cdd = e_0 * (2*np.pi*R_turn*1e-3) * (e_r * theta / np.log(D_wire/d) + 1/np.tan(theta/2) - 1/np.tan(np.pi/12))
Cdd_piko = Cdd * 1e12

Cdd_flip = Cdd_piko         # vector of disk-to-disk capacitance
Cdd_print = Cdd_piko        # first 6 Cdds
for i in range(N_d-2):
    Cdd_print = np.concatenate((Cdd_print, np.flip(Cdd_flip, axis=1)),axis=1)
    Cdd_flip = np.flip(Cdd_flip, axis=1)

# Calculation of disk-to-disk capacitances for displaced disks
a = d_dist/D_wire
if a > 1:
    Cdd_disp = e_0 * (2*np.pi*R_turn*1e-3) * (1/np.sqrt(a**2-1)) * (np.arctan(np.sqrt((a+1)/(a-1))*np.tan(np.pi/4)) - np.arctan(np.sqrt((a+1)/(a-1))*np.tan(-np.pi/4)))
    Cdd_print[:,(d_sep-1)*N_t:d_sep*N_t] = Cdd_disp * 1e12


# %% Calculation of AC resistance coefficient

# Calculation of AC resistance of coils (coef*sqrt(f))
Res_ac_coef = 0.027678/(r_dc_1000ft)**0.5
Res_ac_wire = ((2*np.pi*R_turn*1e-3) * (r_dc_1000ft/(ft_to_m*1e3))) * Res_ac_coef
r_ac = Res_ac_wire

# %% Generating a SPICE netlist for 10 disks
filename = "AX_4_disk_4.0.py.txt"  # name of the txt file
N_d = 10                         # number of disks
cc = 0                          # counter of turns of disk
cur_disk = 1                    # current disk number
with open(filename, "w") as nt:
    for i in range(N_t*N_d):
        nt.write("\nL{} N{}P1 N{}P2 {:0.4f}u\n".format(i+1, i+1, i+1, Lnew[i]))
        pn = 2                  # node point N1P2
        for j in range(N_t*N_d):
            if j!=i:
                pn += 1         # move through node
                nt.write("E{}_{} N{}P{} N{}P{} N{}P{} N{}P{} {:0.4f}\n".format(i+1, j+1, i+1, pn, i+1, pn-1, j+1, 1, j+1, N_t*N_d+1, alpha[i,j]))
        
        nt.write("R%d N%dP%d N%dP%d R=1 LAPLACE=+1/{K%d}/(-S*S/4/3.14^2)^0.25\n" % (i+1, i+1, pn, i+2, 1, i+1))
        nt.write(".param K%d=%0.4e\n" % (i+1, r_ac[:,i-N_t*(cur_disk-1)]))
        
        cc += 1                 # increment counter
        if cc < N_t:
            nt.write("Ctt%d_%d N%dP%d N%dP%d %0.2fp\n" % (i+1, i+2, i+1, 1, i+2, 1, Ctt_p[:,i-N_t*(cur_disk-1)])) #
            
        elif cc == N_t:
            cc = 0              # nullify the counter (new disk)
            Ctt_p = np.flip(Ctt_p, axis=1)
            cur_disk += 1       # increment disk number
            r_ac = np.flip(r_ac, axis=1)
            
    nt.write("\n")
    
    cur_disk = 1                # current disk number
    cc = 0                      # counter of turn
    for i in range(N_t*(N_d-1)):
        nt.write("Cdd%d_%d N%dP%d N%dP%d %0.2fp\n" % (i+1, 2*cur_disk*N_t-i, i+1, 1, 2*cur_disk*N_t-i, 1, Cdd_print[:,i]))
        cc += 1
        if cc == N_t:
            cur_disk += 1
            cc = 0
            
        
    nt.write("\nV1 In 0 AC 10\n")
    nt.write("Rin In N1P1 50\n")
    nt.write("Rout N{}P1 0 50\n".format(N_d*N_t+1))
    nt.write(".print V(N{}P1)\n".format(N_d*N_t+1))
    nt.write(".print V(In)\n")
    nt.write("\n.ac dec 1000 10 200000k\n.backanno\n.end\n")