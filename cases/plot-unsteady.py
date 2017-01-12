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

plt.plot(data[:,0],data[:,1],'-', label="Density")
plt.plot(data[:,0],data[:,2], '-', label="Mach")
plt.plot(data[:,0],data[:,3], '-', label="Pressure")
plt.axis((0,1,-0.1,1.1))
plt.title(scheme)
plt.xlabel("x")
plt.ylabel(" ")
plt.grid("on")
plt.legend(loc = "upper right")
plt.savefig(scheme+".png", format="png")
plt.show()

