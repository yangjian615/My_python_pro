#COPYRIGHT Ramani Kothadia 2015
#PLease maintain this header in all forms of this program
import glob
import math

#Defininf structures
M = 20
TLOW = 286
TUP = 373
Tdel = 3.0
REPS = 30
OBSCOL = 1
PMFTEMP = 307
STARTIDX1 = 54
STARTIDX2 = 56

CVPREFIX1= "Helix1rms18mic.out."
CVPREFIX2 = "Helix2_3rms18mic.out."
CVFOLDER = "/home/krb/New CDAC/remd/shruti_inp/"
CVEXT = "rms"
TEMPFILENAME = "/home/krb/New CDAC/remd/shruti_inp/temps"
TEMPNOSHIFT = 1

temp_dict = {}
filecount1 = []
filecount2 = []
temper = []
psi_mids = []
phi_mids = []
nstates = 0
countArray = [[0 for x in range(M)] for x in range(M)]

#Taking temperature and other inputs from the user
def tempIp():
	global K, observable, filecount1, filecount2
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
	filecount1 = [0 for x in range(K)]
	filecount2 = [0 for x in range(K)]

def obsIpAll():
	global filecount1, filecount2, minvalpsi, maxvalpsi, minvalphi, maxvalphi
	minvalpsi = float('inf')
	maxvalpsi = -float('inf')
	minvalphi = float('inf')
	maxvalphi = -float('inf')
	searchquery = CVFOLDER + "/" + CVPREFIX1 + "*." + CVEXT
	print searchquery
	filelist = glob.glob(searchquery)
	for filename in filelist:
		#print(filename)
		print(filename[STARTIDX1], filename[STARTIDX1:])
		ftors = open(filename, "r")
		idx = STARTIDX1
		ch = filename[idx]
		tempno = ''
		while(ch!='.'):
			tempno = tempno + ch
			idx = idx + 1
			ch = filename[idx]
		tempno = int(tempno) - TEMPNOSHIFT
		print("tempno=", tempno)

		ftors = open(filename, 'r')
		temp = temper[tempno]
		#print("Temp=", temp)
		for line in ftors:
			psi1 = float(line.split()[OBSCOL])
			if psi1  < minvalpsi:
				minvalpsi = psi1
			if psi1 > maxvalpsi:
				maxvalpsi = psi1
			filecount1[tempno] = filecount1[tempno] + 1

		#print(filecount)
		ftors.close()

	searchquery = CVFOLDER + "/" + CVPREFIX2 + "*." + CVEXT 
	filelist = glob.glob(searchquery)
	for filename in filelist:
		#print(filename)
		print(filename[STARTIDX2], filename[STARTIDX2:])
		ftors = open(filename, "r")
		idx = STARTIDX2
		ch = filename[idx]
		tempno = ''
		while(ch!='.'):
			tempno = tempno + ch
			idx = idx + 1
			ch = filename[idx]
		tempno = int(tempno) - TEMPNOSHIFT
		print("tempno=", tempno)
		#print("Tempno=", tempno)
		ftors = open(filename, 'r')
		temp = temper[tempno]
		#print("Temp=", temp)
		for line in ftors:
			phi1 = float(line.split()[OBSCOL])
			if phi1  < minvalphi:
				minvalphi = phi1
			if phi1 > maxvalphi:
				maxvalphi = phi1
			filecount2[tempno] = filecount2[tempno] + 1

		#print(filecount)
		ftors.close()
	#Here in this case filecounts of psi and phi are equal for corresponding values so, we don't mind using a single filecount array for both filecount1 & filecount2

def statineff(nameA, nameB, N, ch, muA, muB):
	statin = 1
	sum_prod = 0
	if ch == 1 or ch == 3:
		fA = open(nameA, 'r')
		for line in fA:
			sum_prod = sum_prod + (float(line.strip().split()[0])- muA) ** 2
		fA.close()
	elif ch == 2:
		fA = open(nameA, 'r')
		fB = open(nameB, 'r')
		for line in fA:
			sum_prod = sum_prod + (float(line.strip().split()[0])-muA)*(float(fB.readline().strip().split()[0])-muB)
		fA.close()
		fB.close()
	sigma2AB = sum_prod/(N-1.0)
	if sigma2AB == 0 :
		return 1
	if ch ==1 or ch == 3:
		fB = open("copy.txt", 'w+')
		fA = open(nameA, 'r')
		for line in fA:
			fB.write(line)
		fB.close()
		fA.close()
	fA = open(nameA, 'r')
	if ch == 1 or ch == 3:
		fB = open('copy.txt', 'r')
	elif ch == 2:
		fB = open(nameB, 'r')
	list1 = fA.readlines()
	list2 = fB.readlines()
	N = min(len(list1), len(list2))
	fA.close()
	fB.close()
	t = 1
	increment = 1
	while(t < N-1):
		fA = open(nameA, 'r')
		if ch == 1 or ch ==3:
			fB = open('copy.txt', 'r')
		elif ch == 2:
			fB = open(nameB, 'r')
		for idx in range(t):
			fB.readline()
		C = 0 
		for idx in range(N-t):
			valA = float(fA.readline().strip().split()[0])
			valB = float(fB.readline().strip().split()[0])
			C = C + valA * valB
		fA.close()
		fB.close()
		fA = open(nameA, 'r')
		if ch == 1 or ch == 3:	
			fB = open('copy.txt', 'r')
		elif ch == 2:
			fB = open(nameB, 'r')
		if ch == 1 or ch == 3:
			C = C/(N-t)
		elif ch == 2:
			for idx in range(t):
				fA.readline()
			for idx in range(N-t-1):
				C = C + float(fA.readline().strip().split()[0]) * float(fB.readline().strip().split()[0])
			C = C/(2.0 * (N-t))
		if ch == 1 or ch == 3:
			muB = muA
		C = (C-muA*muB)/sigma2AB
		if C < 0:	
			return 1
		statin = statin + 2.0 * C * (1.0 - float(t)/float(N)) * increment
		t = t + increment
		increment  = increment + 1
	return statin			
		 	
def computeBins():
	global nstates, nstatearray
	nstatearray = [[0, 0] for x in range(M*M)]
	range_psi = maxvalpsi - minvalpsi
	range_phi = maxvalphi - minvalphi
	delta_psi = abs(range_psi/M)
	delta_phi = abs(range_phi/M)
	for simuno in range(K):
		temp = temper[simuno]
		psifilename = CVFOLDER + "/" + CVPREFIX1 + str(simuno + TEMPNOSHIFT) + "." + CVEXT
		phifilename = CVFOLDER + "/" + CVPREFIX2 + str(simuno + TEMPNOSHIFT) + "." + CVEXT
		fpsi = open(psifilename, "r")
		fphi  = open(phifilename, "r")
		frmsdbin = open("results//cvbin" + str(int(simuno)) + ".txt", "w+")
		for line in fpsi:	
			psival = float(line.strip().split()[OBSCOL])
			phival = float(fphi.readline().strip().split()[OBSCOL])
			binno1 = int(math.floor((psival - minvalpsi)/delta_psi))
			binno2 = int(math.floor((phival- minvalphi)/delta_phi))
			if binno1 ==M:
				binno1 = M-1
			if binno2 ==M:
				binno2 = M-1
			writestring = str(psival) + "\t" + str(binno1) + '\t' + str(phival) + '\t' + str(binno2) + '\n'
			frmsdbin.write(writestring)
		fpsi.close()
		fphi.close()
		frmsdbin.close()

	mids = open('results//cv_mids.txt', 'w+')
	mids.write('cv-1-mid\tcv-2mid\tstate-no\n')
	for bin_idx in range(1, M+1):
		valpsi = round(minvalpsi + delta_psi * (bin_idx-0.5), 4)
		psi_mids.append(valpsi)
		valphi = round(minvalphi + delta_phi * (bin_idx-0.5), 4)
		phi_mids.append(valphi)
	print(phi_mids, psi_mids, 'mids')
	count = -1
	for psimid in psi_mids:
		for phimid in phi_mids:
			count = count + 1
			writestring = str(psimid) + '\t' + str(phimid) + '\t' + str(count) + '\n'
			mids.write(writestring)
	mids.close()

	nstates = M ** 2
	filelist = glob.glob('results//cvbin*.txt')
	for filename in filelist:
		idx = 13
		tempno = ''
		print(filename)
		while(filename[idx]!='.'):
			tempno = tempno + filename[idx]
			idx = idx + 1 
		tempno = int(tempno)
		print(temp, tempno)
		f = open('results//binCV' + str(tempno) + '.txt', 'w+')
		filepointer = open(filename, 'r+')
		for line in filepointer:
			arrayval = line.split()
			psibin = int(arrayval[1])
			phibin = int(arrayval[3])
			countArray[psibin][phibin] = countArray[psibin][phibin] + 1
			state = M * (psibin) + phibin
			#nstatearray[state] = [psi_mids[psibin], phi_mids[phibin]]
			
			writestring = arrayval[0] + '\t' + arrayval[1] + '\t' + arrayval[2] + '\t' + arrayval[3] + '\t' + str(state) + '\n'
			f.write(writestring)
		f.close()
		filepointer.close()
	idx0 = 0	
	for idx1 in range(M):
		for idx2 in range(M):
			nstatearray[idx0] = [psi_mids[idx1], phi_mids[idx2]]
			idx0 = idx0 + 1 
	for idx in countArray:
		print(idx)

def computePMF(temp_of_int):
	beta = 1/(temp_of_int * 0.0098)
	X_i = [0 for x in range(nstates)]	
	Y_i = [0 for x in range(nstates)]
	d2x_i = [0 for x in range(nstates)]
	d2y_i = [0 for x in range(nstates)]
	dxdy_i = [0 for x in range(nstates)]

	max_g_Aw_w_i = [0 for x in range(nstates)]
	P_i = [0 for x in range(nstates)]
	dP_i = [0 for x in range(nstates)]
	F_i = [0 for x in range(nstates)]
	dF_i = [0 for x in range(nstates)]

	for simuno in range(K):
		temp = temper[simuno]
		print(temp)
		wtfilename = 'results//fres' + str(simuno) + '.txt'
		rmsdbinname = 'results//binCV' + str(simuno) + '.txt'
		
		fwt = open(wtfilename, 'r')
		fwtcopy = open('wtcopy.txt', 'w+')
		for line in fwt:
			fwtcopy.write(line)
		fwt.close()
		fwtcopy.close()
		
		fwt = open(wtfilename, 'r')
		sum_wt = 0
		sum_sq = 0
		N = 0
		for line in fwt:
			value = float(line.strip())
			sum_wt = sum_wt + value
			sum_sq = sum_sq + value**2
			N = N +1 
		mu_w = sum_wt/N
		sigma2_w_w = sum_sq/N - mu_w ** 2
		fwt.close()

		sigma2_Aw_Aw_i = [0 for x in range(nstates)]
		sigma2_Aw_w_i = [0 for x in range(nstates)]
		mu_Aw_i = [0 for x in range(nstates)]
	
		fbin = open(rmsdbinname, "r")
		fcopy = open('copybin.txt', 'w+')
		for line in fbin:
			fcopy.write(line)
		fbin.close()
		fcopy.close()

		fbinno = open(rmsdbinname, 'r')
		fwt = open(wtfilename, 'r')
		for snapno in range(N):
			wt = float(fwt.readline().strip())
			line1 = fbinno.readline().split()
			if line1 == []:
				break
			iIndex = int(line1[4])
			mu_Aw_i[iIndex] = mu_Aw_i[iIndex] + wt
			sigma2_Aw_Aw_i[iIndex] = sigma2_Aw_Aw_i[iIndex] + wt**2
			sigma2_Aw_w_i[iIndex] = sigma2_Aw_w_i[iIndex] + wt**2
		fbinno.close()
		fwt.close()
		

		for idx in range(nstates):
			mu_Aw_i[idx] = mu_Aw_i[idx] / N
			sigma2_Aw_Aw_i[idx] = sigma2_Aw_Aw_i[idx]/N - mu_Aw_i[idx]**2
			sigma2_Aw_w_i[idx] = sigma2_Aw_Aw_i[idx]/N - mu_Aw_i[idx]*mu_w
	
		g_w_w = statineff(wtfilename, 0, N, 3, mu_w, 0)
		g_Aw_Aw_i = [1 for x in range(nstates)]
		g_Aw_w_i = [1 for x in range(nstates)]
		t = 1
		increment = 1
		has_crossed_zero_Aw_Aw_i = [False for x in range(nstates)]
		has_crossed_zero_Aw_w_i = [False for x in range(nstates)]
		
		while(t < N-1):
			C_Aw_Aw_i = [0 for x in range(nstates)]
			C_Aw_w_i = [0 for x in range(nstates)]
			fbin = open(rmsdbinname, 'r')
			fbcopy = open('copybin.txt', 'r')
			for idx in range(t):
				fbcopy.readline()
			fwt = open(wtfilename, 'r')
			fwcopy = open('wtcopy.txt', 'r')
			for idx in range(t):	
				fwcopy.readline()
			for t0 in range(N-t):
				iIndex = int(fbin.readline().split()[4])
				one = fbcopy.readline().split()
				if one == []:
					break
				jIndex = int(one[4])
				wt =float(fwt.readline().strip())*float(fwcopy.readline().strip())
				if iIndex == jIndex:
					C_Aw_Aw_i[iIndex] = C_Aw_Aw_i[iIndex] + wt
				C_Aw_w_i[iIndex] = C_Aw_w_i[iIndex] + 0.5 * wt
				C_Aw_w_i[jIndex] = C_Aw_w_i[jIndex] + 0.5 * wt
			fwt.close()
			fwcopy.close()
			fbin.close()
			fcopy.close()

			for idx in range(nstates):	
				C_Aw_Aw_i[idx] = C_Aw_Aw_i[idx] / (N-t)
				C_Aw_w_i[idx] = C_Aw_w_i[idx] / (N-t)
			
			for idx in range(nstates):
				if sigma2_Aw_Aw_i[idx] != 0:
					C_Aw_Aw_i[idx] = (C_Aw_Aw_i[idx] - mu_Aw_i[idx]**2)/sigma2_Aw_Aw_i[idx]
				if sigma2_Aw_w_i[idx] != 0:
					C_Aw_w_i[idx] = (C_Aw_w_i[idx] - mu_Aw_i[idx]*mu_w)/sigma2_Aw_w_i[idx]		

			for idx in range(nstates):	
				if C_Aw_Aw_i[idx] <=0:
					has_crossed_zero_Aw_Aw_i[idx] = True
				if C_Aw_w_i[iIndex] <=0:
					has_crossed_zero_Aw_w_i[idx] = True

			for idx in range(nstates):	
				if has_crossed_zero_Aw_Aw_i[idx] == False:
					g_Aw_Aw_i[idx] = g_Aw_Aw_i[idx] + 2*C_Aw_Aw_i[idx] * (1.0 - float(t)/float(N)) * increment
				if has_crossed_zero_Aw_w_i[idx] == False:
					g_Aw_w_i[idx] = g_Aw_w_i[idx] + 2 * C_Aw_w_i[idx] * (1.0 - float(t)/float(N)) * increment
			
			flag1 = 0
			flag2 = 0 
			for idx in range(nstates):
				if has_crossed_zero_Aw_Aw_i[idx] == False:
					flag1 = 1
					break
			for idx in range(nstates):
				if has_crossed_zero_Aw_w_i[idx] == False:
					flag2 =1
					break
			if flag1 == 0 and flag2 == 0:
				break
			t = t + increment
			increment = increment + 1
	
		for idx in range(nstates):
			if sigma2_Aw_w_i[idx] != 0:
				max_g_Aw_w_i[idx] = math.sqrt(g_Aw_Aw_i[idx] * g_w_w) * abs(math.sqrt(sigma2_Aw_Aw_i[idx] * sigma2_w_w)/sigma2_Aw_w_i[idx])
			else:
				max_g_Aw_w_i[idx] = 0
		for idx in range(nstates):
			if g_Aw_w_i[idx] > max_g_Aw_w_i[idx] and max_g_Aw_w_i[idx] >1:
				g_Aw_w_i[idx] = max_g_Aw_w_i[idx]
		#print("X_i", X_i)
		for idx in range(nstates):
			X_i[idx] = X_i[idx] + mu_Aw_i[idx]
			Y_i[idx] = Y_i[idx] + mu_w
		#print('mu_Aw')
		print(mu_Aw_i)
		#print(N)
		for idx in range(nstates):
			d2x_i[idx] = d2x_i[idx] + sigma2_Aw_Aw_i[idx]/(N/g_Aw_Aw_i[idx])
			d2y_i[idx] = d2y_i[idx] + sigma2_w_w / (N/g_w_w)
			dxdy_i[idx] = dxdy_i[idx] + sigma2_Aw_w_i[idx]/(N/g_Aw_w_i[idx])
	for idx in range(nstates):
		if Y_i[idx] == 0:
			P_i[idx]  = -1
		else:
			P_i[idx] = X_i[idx] / Y_i[idx]
		if X_i[idx] == 0 or Y_i[idx] == 0:
			dP_i[idx] = -1
		elif X_i[idx]**2 == 0 or Y_i[idx]**2 == 0:
		    dP_i[idx] -1
		else:
			print(X_i[idx], Y_i[idx])
			dP_i[idx] = math.sqrt(P_i[idx]**2 * (d2x_i[idx]/(X_i[idx]**2) + d2y_i[idx]/(Y_i[idx]**2) - 2*dxdy_i[idx]/(X_i[idx]*Y_i[idx])))

	for idx in range(nstates):
		if P_i[idx] <= 0:
			F_i[idx] = float('inf')
			dF_i[idx] = float('inf')
		else:
			F_i[idx] = -(1/beta) * math.log(P_i[idx])
			dF_i[idx] = (1/beta) * (dP_i[idx]/P_i[idx])

	minfree = min(F_i)
	for idx in range(nstates):
		F_i[idx] = F_i[idx] - minfree

	for idx in range(nstates):
		P_i[idx] = round(P_i[idx], 4)
		dP_i[idx] = round(dP_i[idx], 4)
		F_i[idx] = round(F_i[idx], 4)
		dF_i[idx] = round(dF_i[idx], 4)
	print(sum(P_i))
	print('nstatearray:', nstatearray)
	fpmf = open('results//2dpmfprob.txt', 'w+')
	fpmf.write('Probabilities & Uncertainties: \n')
	for idx in range(nstates):
		writestring = str(nstatearray[idx][0]) + '\t' + str(nstatearray[idx][1]) + '\t' + str(P_i[idx]) + '\t' + str(dP_i[idx]) + '\n'
		fpmf.write(writestring)
	fpmf.close()
	fpmf = open('results//2dpmffree.txt', 'w+')
	fpmf.write('Free Energies & Uncertainties:\n')
	for idx in range(nstates):	
		writestring = str(nstatearray[idx][0]) + '\t' + str(nstatearray[idx][1]) + '\t' + str(F_i[idx]) + '\t' + str(dF_i[idx]) + '\n'
		fpmf.write(writestring)
	fpmf.close()

def main():
	tempIp()
	obsIpAll()
	computeBins()
	computePMF(PMFTEMP)

if __name__ == '__main__':
	main()
