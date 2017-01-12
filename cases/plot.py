#! /usr/bin/env python3
import sys
import numpy as np
from matplotlib import pyplot as plt

if(len(sys.argv) < 2):
	print("Error. Please provide input file name.")
	sys.exit(-1)
	
fname = sys.argv[1]
scheme = sys.argv[2]

data = np.genfromtxt(fname)
n = data.shape[0]

#plt.plot(data[:,0],data[:,1],'o-', label="Density")
plt.plot(data[:,0],data[:,2], 's-', label="Mach")
plt.plot(data[:,0],data[:,3], 'v-', label="Pres")
#plt.axis((0,1,0.2,1.6))
plt.title(scheme)
plt.xlabel("x")
plt.ylabel(" ")
plt.grid("on")
plt.legend(loc = "upper left")
plt.savefig(scheme+".png", format="png")
plt.show()

