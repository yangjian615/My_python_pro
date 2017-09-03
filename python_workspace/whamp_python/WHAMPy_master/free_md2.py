#Taking free energy values directly from the files

import math
import os
import linecache
import time
MAXITS = 1000
k_boltz = 1.3806503 * 6.0221415 / 4184.0 # Boltzmann constant in kcal/mol/K
TOLERANCE = 1.0e-5
deltaU = 0

#--------------------------------------------------------------------------------------------------------
#Importing models & variables to be used
#--------------------------------------------------------------------------------------------------------
def importmodules():
    global H_mk, filecount, deltaU, M, temper, K, U_mids, H_m
    f = open("histoStore.txt", "r")
    for line in f:
        line = line.split()
        if line[0] == 'M':
            M = int(line[1])
            U_mids = [0 for x in range(M)]
            H_m = [0 for x in range(M)]
        elif line[0] == 'K':
            K = int(line[1])
            temper = [0 for x in range(K)]
            filecount = [0 for x in range(K)]
            H_mk = [[0 for x in range(K)] for x in range(M)]
        elif line[0] == 'temper':
            for idx in range(K):
                temper[idx] = int(float(line[idx+1]))
        elif line[0] == 'filecount':
            for idx in range(K):
                filecount[idx] = int(line[idx+1])
        elif line[0] == 'H_m':
            for idx in range(M):
                H_m[idx] = int(line[idx+1])
        elif line[0] == 'U_mids':
            for idx in range(M):
                U_mids[idx] = float(line[idx+1])
        elif line[0] == 'H_mk':
            for idx in range(M):
                for idx2 in range(K):
                    #print("index = ", idx*K+idx2+1)
                    H_mk[idx][idx2] = int(line[idx*K+idx2+1])
        elif line[0] == 'delta_U':
            deltaU = float(line[1])
    f.close()
    print(deltaU)
    
#--------------------------------------------------------------------------------------------------------
#Defining data structures required in the program
#--------------------------------------------------------------------------------------------------------
def definingArrays():
    global log_omega_m, log_numer_m, log_denom_m, f_l, f_l_old, delta_f_l, temp_f_l
    log_omega_m = [0 for x in range(M)]
    log_numer_m = [0 for x in range(M)]
    log_denom_m = [0 for x in range(M)]
    f_l = [0 for x in range(K)]
    f_l_old = [0 for x in range(K)]
    delta_f_l = [0 for x in range(K)]
    temp_f_l = [0 for x in range(K)]

#--------------------------------------------------------------------------------------------------------
#Defining beta from the given temperatures
#--------------------------------------------------------------------------------------------------------
def calculateBeta():
    global beta_t
    beta_t = []
    for i in range(0, K):
        beta_t.insert(i, 1/(temper[i]*k_boltz))

def computeBinStatIneff():
    print("~~In BinStatIneff Module~~")
    global C_mk, g_mk
    C_mk = [[0 for x in range(K)] for x in range(M)]
    Epsi_mk = [[0 for x in range(K)] for x in range(M)]
    has_crossed_zero = [[False for x in range(K)] for x in range(M)]
    g_mk = [[1 for x in range(K)] for x in range(M)]
    tlag = 1
    increment = 1
    N = min(filecount)
    print("N=", N)
    while(tlag < N-1):
        for i in range(0, M):
            for j in range(0, K):
                C_mk[i][j] = 0
        for simuno in range(0, K):
            filename = "results//fhisto" + str(simuno) + ".txt"
            f = open(filename, "r")
            f.readline()
            for timeori in range(0, N-tlag):
                test = linecache.getline(filename, timeori + 2).split()
                m1 = int(test[2])
                m2 = int(linecache.getline(filename, timeori + tlag + 2).split()[2])
                if m1 == m2:
                    C_mk[m1][simuno] = C_mk[m1][simuno] + 1
        for binno in range(0, M):
            for simuno in range(0, K):
                C_mk[binno][simuno] = C_mk[binno][simuno] / float(N-tlag)
        for binno in range(0, M):
            for simuno in range(0, K):
                if(Epsi_mk[binno][simuno] != 0):
                    C_mk[binno][simuno] = (C_mk[binno][simuno] - Epsi_mk[binno][simuno] ** 2) / (Epsi_mk[binno][simuno] - Epsi_mk[binno][simuno]**2)
                if C_mk[binno][simuno] < 0:
                    has_crossed_zero[binno][simuno] = True
        for binno in range(0, M):
            for simuno in range(0, K):
                if has_crossed_zero[binno][simuno] != True:
                    g_mk[binno][simuno] = g_mk[binno][simuno] + 2.0 * C_mk[binno][simuno] * (1.0- tlag/N) * increment
        flag = 0
        for i in range(0, M):
            for j in range(0, K):
                if(has_crossed_zero[i][j] ==False):
                    flag = 1
                    break
        if flag == 0:
            return 0
        tlag = tlag + increment
        increment = increment + 1
    print("~~Exit BinStatIneff Module~~")
    #print(g_mk)
    return 1
    
def computeEffectiveCounts():
    print("~~In EffectiveCounts Module~~")
    log_zero = -7
    global log_Heff_m, log_Neff_tm
    Heff_m = [0 for x in range(M)]
    log_Heff_m = [log_zero for x in range(M)]
    Neff_tm = [[0 for x in range(K)] for x in range(M)]
    log_Neff_tm = [[log_zero for x in range(K)] for x in range(M)]
    for binno in range(0, M):
        for simuno in range(0, K):
            Heff_m[binno] = Heff_m[binno] + H_mk[binno][simuno]/g_mk[binno][simuno]
        if(Heff_m[binno] > 0):
            log_Heff_m[binno] = math.log(Heff_m[binno])
    for binno in range(0, M):
        for tempno in range(0, K):
            Neff_tm[binno][tempno] = filecount[tempno]/g_mk[binno][tempno]
            if(Neff_tm[binno][tempno] > 0):
                log_Neff_tm[binno][tempno] = math.log(Neff_tm[binno][tempno])
    #print("heff=\n",log_Heff_m, "\n\nNeff=\n",log_Neff_tm)
    print("~~Exit EffectiveCounts Module~~")
    
def logsum(array):
    maxarg = max(array)
    temparray =[]
    for i in array:
        temparray.append(math.e**(i-maxarg))
    logsum = math.log(sum(temparray)) + maxarg
    return logsum

#--------------------------------------------------------------------------------------------------------
#Compute the log density of states from the current dimensionless free energies. 
#--------------------------------------------------------------------------------------------------------
def compute_DensityOfStates():
    arrayK = [0 for x in range(K)]
    global log_omega_m 
    for binno in range(0, M):
        log_numer_m[binno] = log_Heff_m[binno]
        for simuno in range(0, K):
            arrayK[simuno] = log_Neff_tm[binno][simuno] + f_l[simuno] + (-beta_t[simuno]*U_mids[binno]) 
        log_denom_m[binno] = logsum(arrayK)
        log_omega_m[binno] = log_numer_m[binno] - log_denom_m[binno]
        #omega_m[binno] = round(math.e ** log_omega_m[binno], 4) 

#--------------------------------------------------------------------------------------------------------
#Conmpute dimensionless free energies by self-consistent iteration
#--------------------------------------------------------------------------------------------------------
def computingFreeEnergies():
    print("~~In FreeEnergies Module~~")
    arrayK = [0 for x in range(M)]
    for iters in range(0, MAXITS):
        #Save the old values of the free energies
        for idx in range(0, K):
            f_l_old[idx] = f_l[idx]
        compute_DensityOfStates()
        #Initialise f_l to zero
        for idx in range(0, K):
            f_l[idx] = 0
        #Compute new f_l values for each temperature
        for idx in range(0, K):
            for idx1 in range(M):
                arrayK[idx1] = log_omega_m[idx1] + (-beta_t[idx]* U_mids[idx1])
            f_l[idx] = -logsum(arrayK)
            
        tempo = min(f_l)
        #Explanation for shifting energies?
        """for idx, val in enumerate(f_l):
            f_l[idx] = f_l[idx] - tempo"""
            #print(f_l[idx])

        #Calculation of max_delta
        for idx in range(0, K):
            delta_f_l[idx] = f_l[idx] - f_l_old[idx]
        for idx in range(0, K):
            if f_l[idx] != 0:
                temp_f_l[idx] = delta_f_l[idx]/f_l[idx]
        max_delta = max(temp_f_l)
        if max_delta < TOLERANCE:
            break
    if iters<=MAXITS:
        print("Converged to tolerance of ", max_delta, "in", iters+1, " iterations")
    else:
        print("WARNING: Did not converge to within the specified tolerance")
    for idx in range(0, K):
       f_l[idx] = round(f_l[idx], 4)
    """omega_m = [0 for x in range(M)]
    for i in range(M):
        omega_m[i] = math.e ** (log_omega_m[i])"""
    print("Free Energies are:", f_l)
    #print("\n\nLog of Density of States:", log_omega_m)
    #print("\nDensity of states:", omega_m)
    freeE = open("results//freeE.txt", "w+")
    Omega = open("results//Omega.txt", "w+")
    for idx in range(K):
        freeE.write(str(f_l[idx]))
        freeE.write("\n")
    for idx in range(M):
        Omega.write(str(log_omega_m[idx]))
        Omega.write("\n")
    freeE.close()
    Omega.close()
    print("~~Exit FreeEnergies Module~~")

def computingFreeEnergies_ONE():
    print("~Enter Free energy extraction~")
    global f_l
    f = open("diap//f_k1.dat", "r")
    for line in f:
        f_l.append(float(line.split()[0]))
    f.close()
    print("~~Exit FreeEnergies Module~~")
    
def main():
        importmodules()
        definingArrays()
        calculateBeta()
        computeBinStatIneff()
        computeEffectiveCounts()
        #computingFreeEnergies_ONE()
        compute_DensityOfStates()
        computingFreeEnergies()
        #computeWeights()
        print(time.time())
if __name__ == '__main__':
    main()
