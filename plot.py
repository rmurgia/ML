import numpy as np
import matplotlib.pyplot as plt

x,y = np.loadtxt("prova", unpack='True', usecols=[0,1])

plt.loglog(x,y)
plt.show()
