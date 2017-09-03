import numpy as np
import spacepy.LANLstar as LS
import spacepy.time as spt
import spacepy.omni as om
import datetime as datetime

st = datetime.datetime(2013,3,23)
en = datetime.datetime(2013,3,24)
delta = datetime.timedelta(minutes=5)
ticks=spt.tickrange(st,en,delta,'UTC')

data=om.get_omni(ticks)

omni_t=data['ticks']
import numpy as np
from spacepy import pycdf

# cdf = pycdf.CDF('/Users/yangjian/Desktop/rbsp-a_magnetometer_4sec-gsm_emfisis-l3_20130323_v1.3.3.cdf')
# pos_data=cdf['coordinates'][...]
# epoch=cdf['Epoch'][...]
#
# pos_data_new=np.interp(omni_t,epoch,pos_data)

inputdict = {}
inputdict['Kp']     = data['Kp']           # Kp index
inputdict['Dst']    = data['Dst']            # Dst index (nT)
inputdict['dens']   = data['dens']            # solar wind density (/cc)
inputdict['velo']   = data['velo']           # solar wind velocity (km/s)
inputdict['Pdyn']   = data['Pdyn']            # solar wind dynamic pressure (nPa)
inputdict['ByIMF']  = data['ByIMF']            # GSM y component of IMF magnetic field (nT)
inputdict['BzIMF']  = data['BzIMF']             # GSM z component of IMF magnetic field (nT)
inputdict['G1']     = data['G'][:,0]           # as defined in Tsganenko 2003
inputdict['G2']     = data['G'][:,1]
inputdict['G3']     = data['G'][:,2]
inputdict['W1']     = data['W'][:,0]          # as defined in Tsyganenko and Sitnov 2005
inputdict['W2']     = data['W'][:,1]
inputdict['W3']     = data['W'][:,2]
inputdict['W4']     = data['W'][:,3]
inputdict['W5']     = data['W'][:,4]
inputdict['W6']     = data['W'][:,5]
inputdict['Year']   = data['Year']
inputdict['DOY']    = data['DOY']
inputdict['Hr']     = data['Hr']

inputdict['PA']     = np.zeros(len(omni_t))+90            # pitch angle [deg]

# inputdict['Lm']     = data['Hr']             # McIllwain L
# inputdict['Bmirr']  = data['Hr']             # magnetic field strength at the mirror point
# inputdict['rGSM']   = data['Hr']            # radial coordinate in GSM [Re]
# inputdict['lonGSM'] = data['Hr']            # longitude coodrinate in GSM [deg]
# inputdict['latGSM'] = data['Hr']             # latitude coordiante in GSM [deg]
# inputdict['PA']     = data['Hr']            # pitch angle [deg]
# inputdict['SMx']    = data['Hr']
# inputdict['SMy']    = data['Hr']
# inputdict['SMz']    = data['Hr']

# data_out=LS.LANLstar(data, ['OPDYN','OPQUIET','T01QUIET','T01STORM','T89','T96','T05','RAMSCB'])

data_out=LS.LANLmax(inputdict,'T05')

print data_out
np.savetxt(r'/Users/yangjian/Desktop/l_star_max.txt',data_out['T05'])

import pylab as pl
pl.plot(data['BzIMF'])
pl.show()