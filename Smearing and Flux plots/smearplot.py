import numpy as np
import matplotlib.pyplot as plt


a = np.loadtxt("/home/duncan/Desktop/pydata.txt",delimiter=',')

plt.imshow(a)
plt.title('Heat map of CC nu sk2 matrix')
plt.colorbar()

plt.savefig("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/smearingplot.pdf",format='pdf')
plt.show()