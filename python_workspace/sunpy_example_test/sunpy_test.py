# region Description
import sunpy.data.sample
import sunpy.map
aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
aia.peek()
# endregion

# region Description
import numpy as np
import sunpy.data.sample
from sunpy.lightcurve import LightCurve
times = np.arange(1000) * 2.0
signal = np.sin(np.arange(1000)*0.02 ) + np.random.random(1000)
light_curve = LightCurve.create({"signal": signal},index = times)
light_curve.peek()
# endregion

# region Description
import matplotlib.pyplot as plt
import sunpy.spectra
plt.ion
import sunpy.data.sample
from sunpy.spectra.sources.callisto import CallistoSpectrogram
image = CallistoSpectrogram.read(sunpy.data.sample.CALLISTO_IMAGE)
image.peek()
# endregion



# region Description
import sunpy.map
import matplotlib.pyplot as plt
import sunpy.data.sample
plt.ion
aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
fig = plt.figure()
ax = plt.subplot(111)
aia.plot()
aia.draw_limb()
aia.draw_grid()
plt.colorbar()
aia.draw_limb()
plt.show()
# endregion






