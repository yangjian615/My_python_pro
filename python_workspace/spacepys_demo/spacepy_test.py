
# region Description
# coding=utf-8
# from spacepy import coordinates as coord
# cvals = coord.Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
# cvals.x # returns all x coordinates
#
# from spacepy.time import Ticktock
# cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticks
# newcoord = cvals.convert('GSM', 'sph')
# print newcoord

# import spacepy.time as spt
# import spacepy.coordinates as spc
# import spacepy.irbempy as ib
# t = spt.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
# y = spc.Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
# result=ib.get_Bfield(t,y)
# print result
# endregion


# region Description Fetch and Plot OMNI Data
import spacepy.omni as om
import spacepy.time as spt
import matplotlib.pyplot as plt

ticks = spt.tickrange('2007-03-07T00:00:00', \
'2007-05-16T00:00:00', 1./24.)
data = om.get_omni(ticks)

plt.Figure()
ax0 = plt.subplot(211)
plt.plot(data['UTC'], data['Dst'], 'r-')

ax1 = plt.subplot(212)
plt.plot(data['UTC'],data['velo'], 'k-')

ax0.set_ylabel('D$_{st}$ [nT]')
ax1.set_ylabel('V$_{sw}$ [km s$^{1}$]')
ax1.set_xlabel('UTC')
plt.show()

# endregion



#
# # region Description:
# import spacepy.pybats as bats
# obj = bats.Bats2d('filename')
# obj.regrid(0.25, [40, 15], [30,30])
# fig = figure()
# ax = fig.add_subplot(111)
# obj.contourf(ax, 'x', 'y', 'p')
# obj.add_body(ax)
# obj.add_planet_field(ax)
# # endregion

# from spacepy import radbelt as rb
# import datetime as dt
# r = rb.RBmodel()
# starttime = dt.datetime(2003,10,20)
# endtime = dt.datetime(2003,12,5)
# delta = dt.timedelta(minutes=60)
# r.setup_ticks(starttime, endtime, delta)
# r.evolve()
#
# plt=r.plot(clims=[4,11])





import spacepy.time as spt
import spacepy.omni as om
import datetime as datetime
st = datetime.datetime(2005,1,1)
en = datetime.datetime(2009,1,1)
delta = datetime.timedelta(hours=1)
ticks=spt.tickrange(st,en,delta,'UTC')
data=om.get_omni(ticks)




import spacepy.seapy as se
import spacepy.omni as om
import spacepy.toolbox as tb
import datetime as datetime
import numpy as np

epochs =se.readepochs('/Library/Python/2.7/site-packages/spacepy/data/SEA_epochs_OMNI.txt')

st = datetime.datetime(2005,1,1)
en = datetime.datetime(2009,1,1)
einds,oinds = tb.tOverlap([st,en],data['UTC'])
omni1hr = np.array(data['UTC'])[oinds]
delta = datetime.timedelta(hours=1)
window= datetime.timedelta(days=3)
sevx = se.Sea(data['velo'][oinds],omni1hr, epochs, window, delta)
# sevx = se.Sea(data['Pdyn'][oinds],omni1hr, epochs, window, delta)

sevx.sea()
sevx.plot(epochline=True,yquan='V$_{sw}$',xunits='days',yunits='Km s$^{-1}$')






