import glob
import math

temp_dict = {}
filecount = []
psi_mids = []
temper = []
M = 20
STARTIDX = 28
OBSCOL = 1
TEMPFILENAME = "/home/krb/New CDAC/Vinod_US/temps"
FOLDERPATH = "/home/krb/New CDAC/Vinod_US/"
EXT = "rmsd"
TEMPNOSHIFT = 1
PMFTEMP = 307

def tempIp():
    global K, observable
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
            for idx, temp in enumerate(temper):
                temp_dict[temp] = idx
            K = len(temper)
    elif choice == 2:
        print("Enter lower and upper limits [separated by a space]::")
        limits = raw_input().split()
        tlow = float(limits[0])
        tup = float(limits[1])
        print("Enter delta::")
        delta = int(raw_input())
        print("Total number of Replicas[K]::")
        K = int(raw_input())
        for idx in range(K):
            temper.append(float(tlow + idx*delta))
    print("XXX.TEMPERATURE.XXX")

#tors[].txt = xyz[].txt
def obsIp():
    print("===OBSERVABLE EXTRACTION===")
    global filecount, minvalpsi, maxvalpsi, minvalphi, maxvalphi
    for idx in range(K):
        filename = "results//xyz" + str(idx) + ".txt"
        f = open(filename, "w+")
        writestring = "Value\n"
        f.write(writestring)
        filecount.append(0)
        f.close()
    minvalpsi = float("inf")
    maxvalpsi = -float("inf")
    minvalphi = float("inf")
    maxvalphi = -float("inf")

    searchquery = FOLDERPATH + "//*." + EXT
    filelist = glob.glob(searchquery)
    for filename in filelist:
        print(filename[25], filename[26], filename[27], filename[28:])
        ftors = open(filename, "r")
        idx = STARTIDX
        ch = filename[idx]
        tempno = ''
        while(ch!='.'):
            tempno = tempno + ch
            idx = idx + 1
            ch = filename[idx]
        tempno = int(tempno) - TEMPNOSHIFT
        print("tempno=", tempno)
        f = open("results//xyz" + str(tempno) + ".txt", "a+")

        
        for line in ftors:
            psi2 = float(line.split()[OBSCOL])
            if psi2 < minvalpsi :
                minvalpsi = psi2
            if psi2 > maxvalpsi:
                maxvalpsi = psi2 
            writestring = str(psi2) + "\n"
            f.write(writestring)
            #print(filecount)
            filecount[tempno] = filecount[tempno] + 1
        f.close()
        ftors.close()

def statineff(nameA, nameB, N, ch, muA, muB):
    statin = 1
    sum_prod = 0
    if ch == 1 or ch == 3:
        fA = open(nameA, "r")
        for line in fA:
            sum_prod = sum_prod + (float(line.strip().split()[0])-muA)**2
        fA.close()
    elif ch == 2:
        fA = open(nameA, "r")
        fB = open(nameB, "r")
        for line in fA:
            sum_prod = sum_prod + (float(line.strip().split()[0])-muA)*(float(fB.readline().strip().split()[0]) - muB)
        fA.close()
        fB.close()
    sigma2AB = sum_prod/(N-1.0)
    if sigma2AB == 0:
        return 1
    if ch == 1 or ch == 3:
        fB = open("copy.txt", "w+")
        fA = open(nameA, "r")
        for line in fA:
            fB.write(line)
        fB.close()
        fA.close()
    fA = open(nameA, "r")
    if ch == 1 or ch == 3:
        fB = open("copy.txt", "r")
    elif ch == 2:
        fB = open(nameB, "r")
    list1 = fA.readlines()
    list2 = fB.readlines()
    N = min(len(list1), len(list2))
    fA.close()
    fB.close()
    t = 1
    increment = 1
    while(t < N-1):
        fA = open(nameA, "r")
        if ch == 1 or ch == 3:
            fB = open("copy.txt", "r")
        elif ch == 2:
            fB = open(nameB, "r")
        for idx in range(t):
            fB.readline()
        C = 0
        """list1 = fB.readlines()
        print(len(list1))
        print(N-t)
        exit()"""
        for idx in range(N-t):
            valA = float(fA.readline().strip().split()[0])
            valB = float(fB.readline().strip().split()[0])
            C = C + valA * valB
        fA.close()
        fB.close()
        fA = open(nameA, "r")
        if ch == 1 or ch == 3:
            fB = open("copy.txt", "r")
        elif ch == 2:
            fB = open(nameB, "r")
        if ch == 1 or ch == 3:
            C = C/(N-t)
        elif ch == 2:
            for idx in range(t):
                fA.readline()
            for idx in range(N-t-1):
                C = C + float(fA.readline().strip().split()[0])* float(fB.readline().strip().split()[0])
            C = C/(2.0 * (N-t))
        if ch ==1 or ch == 3:
            muB = muA
        C = (C - muA*muB)/sigma2AB
        if C<0:
            return 1
        statin = statin + 2.0 * C * (1.0 - float(t)/float(N)) * increment
        t = t + increment
        increment = increment + 1
    return statin

def computeExpectation():
    X = 0
    Y = 0
    d2x = 0
    d2y = 0
    dxdy = 0
    for simuno in range(K):
        N = filecount[simuno]
        name_wt_val = "results//fweightedvalues" + str(simuno) + ".txt"
        name_wt = "results//fwttemp" + str(simuno) + ".txt"
        name_rmsd = "results//xyz" + str(simuno) + ".txt"
        fwttempval = open(name_wt_val, "w+")
        fwttemp = open(name_wt , "r")
        frmsd = open(name_rmsd, "r")
        frmsd.readline()
        sum_vals = 0
        sum_wt = 0
        for line in fwttemp:
            obsval = float(frmsd.readline().strip().split()[0])
            weight = float(line.strip().split()[0])
            newval = obsval*weight
            writestring = str(newval) + "\n"
            fwttempval.write(writestring)
            sum_vals = sum_vals + newval
            sum_wt = sum_wt + weight
        mu_Aw = sum_vals/N
        mu_w = sum_wt/N
        fwttempval.close()
        fwttemp.close()
        frmsd.close()
        
        fwttempval = open(name_wt_val, "r")
        fwttemp = open(name_wt, "r")
        #following line probably not necessary
        #frmsd = prmsd
        sum_vals = 0
        sum_prod = 0
        sum_wt = 0
        for line in fwttempval:
            c_wtval = float(line.strip().split()[0])
            c_wt = float(fwttemp.readline().strip().split()[0])
            sum_vals = sum_vals + (c_wtval - mu_Aw)**2
            sum_prod = sum_prod + (c_wtval-mu_Aw)*(c_wt - mu_w)
            sum_wt = sum_wt + (c_wt - mu_w)**2
        sigma2_Aw_Aw = sum_vals/(N-1)
        sigma2_Aw_w = sum_prod/(N-1)
        sigma2_w_w = sum_wt/(N-1)
        fwttempval.close()
        fwttemp.close()

        g_Aw_Aw = statineff(name_wt_val, 0, N, 1, mu_Aw, 0)
        g_Aw_w = statineff(name_wt_val, name_wt, N, 2, mu_Aw, mu_w)
        g_w_w = statineff(name_wt, 0, N, 3, mu_w, 0)
        
        print(g_Aw_Aw, g_w_w)
        if sigma2_Aw_w!=0:
            max_g_Aw_w = math.sqrt(g_Aw_Aw*g_w_w) * abs(math.sqrt(sigma2_Aw_Aw * sigma2_w_w)/sigma2_Aw_w)
        else:
            max_g_Aw_w = 0
        if (g_Aw_w > max_g_Aw_w) and (max_g_Aw_w > 1):
            g_Aw_w = max_g_Aw_w

        X = X + mu_Aw
        Y = Y + mu_w

        d2x = d2x + sigma2_Aw_Aw /(N/g_Aw_Aw)
        d2y = d2y + sigma2_w_w / (N/g_w_w)
        dxdy = dxdy + sigma2_Aw_w / (N/g_Aw_w)

    A = X/Y
    dA = math.sqrt(A**2 * (d2x /(X**2) + d2y/(Y**2) - 2*dxdy/(X*Y)))
    print("Expected value of observable:", A, "\nStatistical Uncertainty:", dA)

def computeBins():
    range_psi = maxvalpsi - minvalpsi
    #delta_psi = abs(range_obs/M)
    delta_psi = abs(range_psi/M)
    for simuno in range(K):
        temp = temper[simuno]
        frmsd = open("results//xyz" + str(simuno) + ".txt", "r")
        frmsdbin = open("results//xyzbin" + str(simuno) + ".txt", "w+")
        frmsd.readline()
        for line in frmsd:
            value = float(line.strip().split()[0])
            binno = int(math.floor((value - minvalpsi)/delta_psi))
            if binno == M:
                binno = M-1
            writestring = str(binno) + "\n"
            frmsdbin.write(writestring)
        frmsd.close()
        frmsdbin.close()
    mids = open("results//obs_mids.txt", "w+")
    mids.write("psi-mid\tstate-no\n")
    for bin_idx in range(1, M+1):
        valpsi = round(minvalpsi + delta_psi * (bin_idx-0.5), 4)
        psi_mids.append(valpsi)
    #print(psi_mids)
    count = -1
    for psimid in psi_mids:
        count = count + 1
        writestring = str(psimid) + "\t" + str(count) + "\n"
        mids.write(writestring)
    mids.close()
    

def computePMF(temp_of_int):
    beta = 1/(temp_of_int* 0.0098)
    
    X_i = [0 for x in range(M)]
    Y_i = [0 for x in range(M)]
    d2x_i = [0 for x in range(M)]
    d2y_i = [0 for x in range(M)]
    dxdy_i = [0 for x in range(M)]

    max_g_Aw_w_i = [0 for x in range(M)]
    P_i = [0 for x in range(M)]
    dP_i = [0 for x in range(M)]
    F_i = [0 for x in range(M)]
    dF_i = [0 for x in range(M)]
    
    for simuno in range(K):
        temp = temper[simuno]
        wtfilename = "results//fwttemp" + str(simuno) + ".txt"
        rmsdbinname = "results//xyzbin" + str(simuno) + ".txt"
        fwttemp = open(wtfilename, "r")
        fwttempcopy = open("wtcopy.txt", "w+")
        for line in fwttemp:
            fwttempcopy.write(line)
        fwttemp.close()
        fwttempcopy.close()
        fwttemp = open(wtfilename, "r")
        sum_wt = 0
        sum_sq = 0
        N = 0
        for line in fwttemp:
            value = float(line.strip())
            sum_wt = sum_wt + value
            sum_sq = sum_sq + value**2
            N = N + 1
        mu_w = sum_wt/N
        sigma2_w_w = sum_sq/N - mu_w**2
        fwttemp.close()

        #Test statements
        #print("mu_w", mu_w)
        #print("sigma_w_w", sigma2_w_w)
        
        sigma2_Aw_Aw_i = [0 for x in range(M)]
        sigma2_Aw_w_i = [0 for x in range(M)]
        mu_Aw_i = [0 for x in range(M)]
        
        fbin = open(rmsdbinname, "r")
        fcopy = open("copybin.txt", "w+")
        for line in fbin:
                fcopy.write(line)
        fbin.close()
        fcopy.close()

        fbinno = open(rmsdbinname, "r")
        fwttemp = open(wtfilename, "r")
        for snapno in range(N):
            wt = float(fwttemp.readline().strip())
            line1 = fbinno.readline().split()
            #print(line1, snapno)
            if line1 == []:
                break
            iIndex = int(line1[0])
            #print(iIndex)
            mu_Aw_i[iIndex] = mu_Aw_i[iIndex] + wt
            sigma2_Aw_Aw_i[iIndex] = sigma2_Aw_Aw_i[iIndex] + wt**2
            sigma2_Aw_w_i[iIndex] = sigma2_Aw_w_i[iIndex] + wt**2
        fbinno.close()
        fwttemp.close()
        for idx in range(M):
            mu_Aw_i[idx] = mu_Aw_i[idx] / N
            sigma2_Aw_Aw_i[idx] = sigma2_Aw_Aw_i[idx]/N - mu_Aw_i[idx]**2
            sigma2_Aw_w_i[idx] = sigma2_Aw_Aw_i[idx]/N - mu_Aw_i[idx]*mu_w

        g_w_w = statineff(wtfilename, 0, N, 3, mu_w, 0)
        #Test
        #print("g_w_w", g_w_w)
        #print("mu_Aw_i", mu_Aw_i)
        
        g_Aw_Aw_i = [1 for x in range(M)]
        g_Aw_w_i = [1 for x in range(M)]
        t = 1
        increment = 1
        has_crossed_zero_Aw_Aw_i = [False for x in range(M)]
        has_crossed_zero_Aw_w_i = [False for x in range(M)]


        while(t <N-1):
            C_Aw_Aw_i = [0 for x in range(M)]
            C_Aw_w_i = [0 for x in range(M)]

            #Test
            #Works fine. All C functions getting initialized to zero. No problem here
            fbin = open(rmsdbinname, "r")
            fbcopy = open("copybin.txt", "r")
            for idx in range(t):
                fbcopy.readline()
            fwttemp = open(wtfilename, "r")
            fwcopy = open("wtcopy.txt", "r")
            for idx in range(t):
                fwcopy.readline()
            for t0 in range(N-t):
                row1 = fbin.readline().split()
                row2 = fbcopy.readline().split()
                if row1 == [] or row2 == []:
                    break
                iIndex = int(row1[0])
                jIndex = int(row2[0])
                
                wt = float(fwttemp.readline().strip())*float(fwcopy.readline().strip())
                if iIndex == jIndex:
                    C_Aw_Aw_i[iIndex] = C_Aw_Aw_i[iIndex] + wt

                C_Aw_w_i[iIndex] = C_Aw_w_i[iIndex] + 0.5 * wt
                C_Aw_w_i[jIndex] = C_Aw_w_i[jIndex] + 0.5 * wt
                
            fwttemp.close()
            fwcopy.close()
            fbin.close()
            fcopy.close()

            for idx in range(M):
                C_Aw_Aw_i[idx] = C_Aw_Aw_i[idx]/(N-t)
                C_Aw_w_i[idx] = C_Aw_w_i[idx] /(N-t)
            #Test
            #print("C_Aw_Aw_i", C_Aw_Aw_i)
            #print("C_Aw_w_i", C_Aw_w_i)

            for idx in range(M):
                if sigma2_Aw_Aw_i[idx] != 0:
                    C_Aw_Aw_i[idx] = (C_Aw_Aw_i[idx] - mu_Aw_i[idx]**2)/ sigma2_Aw_Aw_i[idx]
                if sigma2_Aw_w_i[idx] != 0:
                    C_Aw_w_i[idx] = (C_Aw_w_i[idx] - mu_Aw_i[idx]*mu_w)/ sigma2_Aw_w_i[idx]

            #Test
            #print("C_Aw_Aw_i", C_Aw_Aw_i)
            #print("C_Aw_w_i", C_Aw_w_i)
            
            for idx in range(M):
                if C_Aw_Aw_i[idx] <= 0:
                    has_crossed_zero_Aw_Aw_i[idx] = True
                if C_Aw_w_i[idx] <= 0:
                    has_crossed_zero_Aw_w_i[idx] = True

            #print(has_crossed_zero_Aw_Aw_i, has_crossed_zero_Aw_w_i)
            for idx in range(M):
                if has_crossed_zero_Aw_Aw_i[idx] == False:
                    g_Aw_Aw_i[idx] = g_Aw_Aw_i[idx] + 2 * C_Aw_Aw_i[idx] * (1.0 - float(t)/float(N)) * increment
                if has_crossed_zero_Aw_w_i[idx] == False:
                    g_Aw_w_i[idx] = g_Aw_w_i[idx] + 2 * C_Aw_w_i[idx] * (1.0 - float(t)/float(N)) * increment
            
            

            
            flag1 = 0
            flag2 = 0
            for idx in range(M):
                if has_crossed_zero_Aw_Aw_i[idx] == False:
                    flag1 = 1
                    break
            for idx in range(M):
                if has_crossed_zero_Aw_w_i[idx] == False:
                    flag2 = 1
                    break
            #print(flag1, flag2)
            if flag1 == 0 and flag2 == 0:
                break
            t = t + increment
            increment = increment + 1

        for idx in range(M):
            if sigma2_Aw_w_i[idx] !=0:
                max_g_Aw_w_i[idx] = math.sqrt(g_Aw_Aw_i[idx]*g_w_w) * abs(math.sqrt(sigma2_Aw_Aw_i[idx]* sigma2_w_w)/sigma2_Aw_w_i[idx])
            else:
                max_g_Aw_w_i[idx] = 0
        for idx in range(M):
            if g_Aw_w_i[idx] > max_g_Aw_w_i[idx] and max_g_Aw_w_i[idx] > 1:
                g_Aw_w_i[idx] = max_g_Aw_w_i[idx]

        for idx in range(M):
            X_i[idx] = X_i[idx] + mu_Aw_i[idx]
            Y_i[idx] = Y_i[idx] + mu_w

        print("X_i[", idx, "]", X_i[idx])
        print("Y_i[", idx, "]", Y_i[idx])

        for idx in range(M):
            d2x_i[idx] = d2x_i[idx] + sigma2_Aw_Aw_i[idx]/(N/g_Aw_Aw_i[idx])
            d2y_i[idx] = d2y_i[idx] + sigma2_w_w /(N/g_w_w)
            dxdy_i[idx] = dxdy_i[idx] + sigma2_Aw_w_i[idx] / (N/g_Aw_w_i[idx])

    for idx in range(M):
        if Y_i[idx] == 0:
            P_i[idx] = -1
        else:
            P_i[idx] = X_i[idx] / Y_i[idx]
        if X_i[idx] == 0 or Y_i[idx] == 0:
            dP_i[idx] = -1
        else:
            dP_i[idx] = math.sqrt(P_i[idx]**2 * (d2x_i[idx]/(X_i[idx]**2) + d2y_i[idx]/(Y_i[idx]**2) - 2*dxdy_i[idx]/(X_i[idx]*Y_i[idx])))

    for idx in range(M):
        if P_i[idx] <= 0:
            F_i[idx] = float("inf")
            dF_i[idx] = float("inf")
        else:
            F_i[idx] = -(1/beta) * math.log(P_i[idx])
            dF_i[idx] = (1/beta) * (dP_i[idx]/P_i[idx])

    minfree = min(F_i)
    for idx in range(M):
        F_i[idx] = F_i[idx] - minfree

    print(sum(P_i))
    fpmf = open("results//pmfProb.txt", "w+")
    fpmf.write("Probabilities\tUncertainty:\n")
    for idx in range(M):
        writestring = str(round(P_i[idx],6)) + "\t" + str(round(dP_i[idx], 6)) + "\n"
        fpmf.write(writestring)
    fpmf.close()
    fpmf = open("results//pmfFree", "w+")
    fpmf.write("Free Energies\tUncertainty:\n")
    for idx in range(M):
        writestring = str(round(F_i[idx], 6)) + "\t" + str(round(dF_i[idx], 6)) + "\n"
        fpmf.write(writestring)
    fpmf.close()

def resultfile():
    f = open("results//obs_mids.txt", "r")
    prob = open("results//pmfProb.txt", "r")
    free = open("results//pmfFree", "r")
    prob.readline()
    free.readline()
    final1 = open("results//final1.txt", "w+")
    final2 = open("results//final2.txt", "w+")
    writestring = "psi_mid\tprob\tuncertainty"
    final1.write(writestring)
    writestring = "psi_mid\tfreeval\tuncertainty"
    final2.write(writestring)
    f.readline()
    for line in f:
        psi_mid = line.split()[0]
        line1 = prob.readline().split()
        probval = line1[0]
        probun = line1[1]
        line2 = free.readline().split()
        print(line2)
        freeval = line2[0]
        freeun = line2[1]
        writestring1 = psi_mid + "\t" + probval + "\t" + probun + "\n"
        writestring2 = psi_mid + "\t" + freeval + "\t" + freeun + "\n"
        final1.write(writestring1)
        final2.write(writestring2)
    final1.close()
    final2.close()
    prob.close()
    free.close()
    f.close()

    
                
def main():
	tempIp()
	obsIp()
	#computeExpectation()
	computeBins()
	computePMF(PMFTEMP)
	resultfile()

if __name__ == '__main__':
	main()
