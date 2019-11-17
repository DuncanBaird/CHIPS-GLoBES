import numpy as np
import matplotlib.pyplot as plt

x,y = np.loadtxt("histsim.dat",unpack=True)

plt.figure()
plt.hist(x,len(x),weights=y)

plt.savefig('quickhistCHIPS.pdf',)