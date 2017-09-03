import numpy as np
import spacepy.LANLstar as LS
import spacepy.time as spt
import spacepy.omni as om
import datetime as datetime

data=np.loadtxt('/Users/yangjian/Desktop/temporary/QinDenton_20130323_1min.txt')

inputdict = {}
inputdict['Kp']     = data[:,0]           # Kp index
inputdict['Dst']    = data[:,1]            # Dst index (nT)
inputdict['dens']   = data[:,2]           # solar wind density (/cc)
inputdict['velo']   = data[:,3]          # solar wind velocity (km/s)
inputdict['Pdyn']   = data[:,4]           # solar wind dynamic pressure (nPa)
inputdict['ByIMF']  = data[:,5]           # GSM y component of IMF magnetic field (nT)
inputdict['BzIMF']  = data[:,6]             # GSM z component of IMF magnetic field (nT)
inputdict['G1']     = data[:,7]           # as defined in Tsganenko 2003
inputdict['G2']     = data[:,8]
inputdict['G3']     = data[:,9]
inputdict['W1']     = data[:,10]          # as defined in Tsyganenko and Sitnov 2005
inputdict['W2']     = data[:,11]
inputdict['W3']     = data[:,12]
inputdict['W4']     = data[:,13]
inputdict['W5']     = data[:,14]
inputdict['W6']     = data[:,15]
inputdict['Year']   = np.zeros(len(data[:,0]))+2013
inputdict['DOY']    = np.zeros(len(data[:,0]))+82
inputdict['Hr']     = data[:,16]

inputdict['PA']     = np.zeros(len(data[:,0]))+90            # pitch angle [deg]



# data_out=LS.LANLstar(data, ['OPDYN','OPQUIET','T01QUIET','T01STORM','T89','T96','T05','RAMSCB'])

data_out=LS.LANLmax(inputdict,'T05')

print data_out
np.savetxt(r'/Users/yangjian/Desktop/l_star_max.txt',data_out['T05'])

import pylab as pl
pl.subplot(511)
pl.plot(inputdict['BzIMF'])
pl.xlim(0,1440)
pl.ylabel('BzIMF')
pl.xticks([0,100,200],['',''])


pl.subplot(512)
pl.plot(inputdict['ByIMF'])
pl.ylabel('ByIMF')
pl.xlim(0,1440)

pl.subplot(513)
pl.plot(inputdict['Dst'])
pl.ylabel('Dst')
pl.xlim(0,1440)

pl.subplot(514)
pl.plot(inputdict['velo'])
pl.ylabel('velo')
pl.xlim(0,1440)

pl.subplot(515)
pl.plot(data_out['T05'])
pl.ylabel('L*')
pl.xlim(0,1440)

pl.show()