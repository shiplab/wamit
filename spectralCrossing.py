from seaspectrum import seaspec, legenda
import numpy as np
import matplotlib.pyplot as plt
import analysis_wamit


#plt.style.use('grayscale')

print(plt.style.available)

Tp = np.arange(3.5, 10.5, .5)
Hs = 2*np.ones(Tp.shape)
Sw, w = seaspec(Tp, Hs)

if isinstance(Tp, np.ndarray):
    y = w.T*np.ones((len(Tp), 1))
    Hs_check = 4 * np.sqrt(np.trapz(Sw.T, y))
else:
    y = w
    Hs_check = 4 * np.sqrt(np.trapz(Sw.T, y.T))

print(Hs_check)

fig = plt.figure()
plt.plot(2*np.pi/w, Sw)
plt.xlim([0, 25])
plt.legend(legenda(Hs_check))
plt.xlabel('Period (s)')
plt.ylabel('PSD [m.s²]')

plt.show()

rao = analysis_wamit.raos(1)