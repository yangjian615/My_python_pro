#COPYRIGHT Ramani Kothadia 2015
#PLease maintain this header in all forms of this program
import glob
import math
import os
import time

U_mids = []
deltaU = 1
M=500
temp_dict = {}
temper = []
OBSCOL = 0

def defineArrays():
    global psi_mkn, m_kn, H_mk, H_m
    H_mk = [[0 for x in range(K)] for x in range(M)]
    H_m = [0 for x in range(M)]

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Reading from standard output file: *inputStore.txt* 
#-------------------------------------------------------------------------------------------------------------------------------------------------
def importmodules():
    f = open("inputStore.txt", "r")
    global K, filecount, temp_dict, temper, U_min, U_max
    for line in f:
        line = line.split()
        if line[0] == 'U_min':
            U_min = float(line[1])
        elif line[0] == 'U_max':
            U_max = float(line[1])
        elif line[0] == 'K':
            K = int(line[1])
            temper = [0 for x in range(K)]
            filecount = [0 for x in range(K)]
        elif line[0] == 'temper':
            for idx in range(K):
                temper[idx] = float(line[idx+1])
        elif line[0] == 'filecount':
            for idx in range(K):
                filecount[idx] = int(line[idx+1])
    f.close()
    for idx, temp in enumerate(temper):
        temp_dict[temp] = idx
    #print("filecount=", filecount)
    #print(U_min, U_max, K, temper, filecount, temp_dict)

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Computing M bins of the observable U 
#-------------------------------------------------------------------------------------------------------------------------------------------------
def computeBins():
    print("~~In Binning Module~~")
    global deltaU, U_min, U_mids, psi_mkn, m_kn, K_nm, H_mk, H_m
    
    deltaU = abs((U_max- U_min)/M)

    #Filling U_mids
    for bin_idx in range(1, M+1):
        U_mids.append(round(U_min + deltaU * (bin_idx-0.5), 4))
    print("~~Exit Binning Module~~")

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Calculating Histograms and Main Processing
#-------------------------------------------------------------------------------------------------------------------------------------------------
def readFromFiles():
    print("~~In Histogram Main Processing Module~~")
#Listing processed standard input files
    filelist = glob.glob("results//fobs*.txt")
    global H_mk
#Creating new output bin files
    for i in range(M):
        binfile = open("results//bin" + str(i) + ".txt", "w+")
        binfile.close()

#Processing each file in the listing
    for filename in filelist:
        H_k = [0 for x in range(M)]
        infile = open(filename, 'r+')
        temp_idx = int(infile.readline().strip())
        temp = temper[temp_idx]
        """temp = int(infile.readline().strip())
        temp_idx = temp_dict[temp]"""
#Ignoring column headers
        infile.readline()

#Creating new output histogram files corresponding to each temperature
        fhisto = open("results//fhisto" + str(temp_idx)+ ".txt", "w+")
        fhisto.write("SNO\t\tVal\tBIN\tStr\tStep")
        snapno = -1
        for line in infile:
            line = line.strip().split()
#Storing observable in variable 'val'
            val = float(line[OBSCOL])
            snapno = snapno + 1
#Calculating binno for that value[; taking care of boundary conditions]
            binno = int(math.floor((val - U_min)/deltaU))
            if binno == M:
                binno = M-1
#Writing that value in the corresponding bin file
            binfile = open("results//bin" + str(binno) + ".txt", "a+")
            writestring = str(snapno) + "\t" + str(temp_idx) + "\n"
            binfile.write(writestring)
            binfile.close()
#Update Histogram array for that temperature and bin number
            H_k[binno] = H_k[binno] + 1
#Write value in the standard histogram output file fhisto
            writestring = "\n" + str(snapno) + "\t" + line[0] + "\t" + str(binno)
            fhisto.write(writestring)
        infile.close()
        #os.remove(filename)
        fhisto.close()

#Update global H_mk array by copying temporary H_k array in the H_mk[][]'s temperature row
        for binno, Hcount in enumerate(H_k):
            H_mk[binno][temp_idx] = Hcount

#Calculating H_m
    for binno in range(M):
        H_m[binno] = sum(H_mk[binno])
        #print("binno:", binno, "=>", H_m[binno])
    print("H_m:", H_m)
    print("~~Exit Histogram Main Processing Module~~")

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Creating a standard output file for this program; to be used by other programs; format specified in the README
#-------------------------------------------------------------------------------------------------------------------------------------------------
def store():
    print("~~Creating output files~~")
    f = open("histoStore.txt", "w+")
    ws = "M\t" + str(M)
    f.write(ws)
    ws = "\nK\t" + str(K)
    f.write(ws)
    ws = "\nU_mids\t"
    f.write(ws)
    for idx in range(M):
        f.write(str(U_mids[idx]))
        f.write("\t")
    ws = "\ndelta_U\t" + str(deltaU)
    f.write(ws)
    ws = "\nH_m\t"
    f.write(ws)
    for idx in range(M):
        f.write(str(H_m[idx]))
        f.write("\t")
    ws = "\nH_mk\t"
    f.write(ws)
    for idx1 in range(M):
        for idx2 in range(K):
            f.write(str(H_mk[idx1][idx2]))
            f.write("\t")
    ws = "\ntemper\t"
    f.write(ws)
    for idx in range(K):
        f.write(str(temper[idx]))
        f.write("\t")
    ws = "\nfilecount\t"
    f.write(ws)
    for idx in range(K):
        f.write(str(filecount[idx]))
        f.write("\t")
    f.close()
    print("~~Successful creation of output files~~")

def main():
    print("***Hi, We're in Histogram Processing Program.***")
    importmodules()
    defineArrays()
    computeBins()
    readFromFiles()
    #writeFiles()
    #constructHistogram()
    store()
    print("***Hope, you didn't have to wait for long. Bye!***")

if __name__ == '__main__':
    main()
