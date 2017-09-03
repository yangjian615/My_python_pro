import math
k_boltz = 1.3806503 * 6.0221415 / 4184.0 # Boltzmann constant in kcal/mol/K


def calculateBeta():
    global beta_t
    beta_t = []
    for i in range(0, K):
        beta_t.insert(i, 1/(temper[i]*k_boltz))

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
    
    global log_omega_m
    log_omega_m = [0 for x in range(M)]
    f = open("results//Omega.txt", "r")
    idx = 0
    for line in f:
        if line != "":
            log_omega_m[idx] = float(line.strip())
            idx = idx + 1
    #print(log_omega_m)
    f.close()
    print(deltaU)


#--------------------------------------------------------------------------------------------------------
#Computing weights for each energy snapshot in every simulation 
#--------------------------------------------------------------------------------------------------------
def computeWeights():
#Calculating weights and setting them in the array w_kn
    maxval = -float("inf")
    wsum = 0
    for simuno in range(K):
        tot = filecount[simuno]
        temp = temper[simuno]
        simufile = open("results//fhisto" + str(simuno) + ".txt", "r")
        print simufile.readline()
        line = simufile.readline().strip()
        #resultfile = open("REMDResults//fres" + str(temp) + ".txt", "w+")
        line = line + "\tWeight\n"
        whtfile = open("results//fwttemp" + str(simuno) + ".txt", "w+")
        #resultfile.write(line)
        for snapno in range(tot):
            charline= simufile.readline().strip()
            line = charline.split()
            if line == []:
                break
            binno = int(line[2])
            if H_m[binno] != 0:
                #print(math.e ** (-beta_t[simuno] * U_mids[binno]))
                weight = log_omega_m[binno] - (beta_t[simuno] * U_mids[binno]) + math.log(H_m[binno])
                wsum = wsum + weight
                if weight > maxval:
                    maxval = weight
                weight = round(weight, 4)
            else:
                #print("In else:")
                weight = 0
            #charline = charline + "\t" + str(weight) + "\n"
            writestring = str(weight) + "\n"
            whtfile.write(writestring)
            #resultfile.write(charline)
        simufile.close()
        #resultfile.close()
        whtfile.close()
    
    for simuno in range(K):
        tot = filecount[simuno]
        temp = temper[simuno]
        resultfile = open("results//fres" + str(simuno) + ".txt", "w+")
        whtfile = open("results//fwttemp" + str(simuno) + ".txt", "r")
        for snapno in range(tot):
            charline = whtfile.readline().strip()
            line = charline.split()
            if line == []:
                break
            weight = float(line[0])
            weight = (math.e ** (weight - maxval)) / wsum
            charline = str(weight) + "\n"
            resultfile.write(charline)
        resultfile.close()
        whtfile.close()

def main():
        importmodules()
        calculateBeta()
        computeWeights()
        #print(time.time())
if __name__ == '__main__':
    main()
