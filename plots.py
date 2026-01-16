import numpy as np
import matplotlib.pyplot as plt

def plotdat(dat):
  data = np.loadtxt(dat, skiprows=6)
  age = data[:,2]
  mass = data[:,4]
  plt.plot(age, mass)
  #plt.semilogy(age, mass)

f1 = "LOGS_rho4em17_y80_mod/history.data"
f2 = "LOGS_rho4em17_y80_base/history.data"

plotdat(f1)
plotdat(f2)
plt.show()
