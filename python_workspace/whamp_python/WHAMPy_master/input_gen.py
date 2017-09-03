import glob

nmaxE = -float("inf")
nminE = float("inf")

tlow = 300
tup = 300
delta = 1.0
K = 1

temper = []
temp_dict = {}

INPUTFOLDER = "/Users/yangjian/Documents/python_workspace/whamp_python/WHAMPy-master/inputfolder/"
TEMPFILENAME = "/Users/yangjian/Documents/python_workspace/whamp_python/WHAMPy-master/inputfolder/temps"
ENERCOL = 0
INPUTFILENAME = "Ex1.en"
INPUTEXT = "en"
REPSHIFT = 0

def temp_ip():
    global K, observable, maxvalArray, minvalArray, avgvalArray, filecount
    print("===TEMPERATURE MODULE===")
    print("Temp I/p Method: 1.Temp File 2.Temp Range [Enter 1/2] ::")
    choice = int(raw_input())
    if choice == 1:
        if TEMPFILENAME == "":
            print("\'tempfilename\' is empty. Please assign a value to it.")
        else:
            tempfile = open(TEMPFILENAME, "r")
            listt = tempfile.readline().split()
            for tempc in listt:
                temper.append(float(tempc))
            K = len(temper)
    elif choice == 2:
        print("Enter lower and upper limits [separated by a space]::")
        limits = raw_input().split()
        tlow = float(limits[0])
        tup = float(limits[1])
        print("Enter delta::")
        delta = float(raw_input())
        print("Total number of Replicas[K]::")
        K = int(raw_input())
        for idx in range(K):
            temper.append(float(tlow + idx*delta))
    if temper != []:
        for idx, temp in enumerate(temper):
            temp_dict[temp] = idx
    minvalArray = [float("inf") for x in range(K)]
    maxvalArray = [-float("inf") for x in range(K)]
    avgvalArray = [0 for x in range(K)]
    filecount = [0 for x in range(K)]
    print temper
    print("XXX.TEMPERATURE.XXX")

def obs_ip():
    print("===In Obs_IP===")
    global minvalArray, maxvalArray, avgvalArray
    searchquery = INPUTFOLDER + "*." + INPUTEXT
    print(searchquery)
    filelist = glob.glob(searchquery)
    print("filelist", filelist)
    for inputfile in filelist:
        minE = float("inf")
        maxE = -float("inf")
        fip = open(inputfile, "r")
        rep = int(fip.readline()) - REPSHIFT
        print rep
        """temp = int(fip.readline())
        rep = temp_dict[temp]"""
        fip.readline() #Header
        for line in fip:
            val = float(line.strip().split()[ENERCOL])
            if val < minE:
                minE = val
            if val > maxE:
                maxE = val
            filecount[rep] += 1
        minvalArray[rep] = minE
        maxvalArray[rep] = maxE
        avgvalArray[rep] = (minE + maxE)/2
        fip.close()
    print("XXX.Obs_ip.XXX")

def normalize():
    global nminE, nmaxE
    searchquery = INPUTFOLDER + "*." + INPUTEXT
    filelist = glob.glob(searchquery)
    for inputfile in filelist:
        fip = open(inputfile, "r")
        rep = int(fip.readline()) - REPSHIFT
        """temp = int(fip.readline())
        rep = temp_dict[temp]"""
        avg = avgvalArray[rep]
        outputfile = "results//fobs" + str(rep) + ".txt"
        fop = open(outputfile, "w+")
        fop.write(str(rep))
        fop.write("\nNorm-Energy\n")
        fip.readline() #Header
        for line in fip:
            val = float(line.strip().split()[ENERCOL])
            newval = val - avg
            fop.write(str(newval))
            fop.write("\n")
            if newval < nminE:
                nminE = newval
            if newval > nmaxE:
                nmaxE = newval
        fip.close()
        fop.close()

def input_store():
    #Write inputstore.txt for use in the next module: histo module
    fis = open("inputStore.txt", "w+")
    ws = "U_min\t" + str(nminE) + "\n"
    fis.write(ws)
    ws = "U_max\t" + str(nmaxE) + "\n"
    fis.write(ws)
    ws = "K\t" + str(len(temper)) + "\n"
    fis.write(ws)
    ws = "temper\t"
    for temp in temper:
        ws = ws + str(temp) + "\t"
    ws = ws + "\n"
    fis.write(ws)
    ws = "filecount\t"
    for count in filecount:
        ws = ws + str(count) + "\t"
    ws = ws + "\n"
    fis.write(ws)
    fis.close()

def main():
    temp_ip()
    obs_ip()
    normalize()
    input_store()

if __name__ == '__main__':
    main()
