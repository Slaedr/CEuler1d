#! /usr/bin/env python
import sys
import numpy as np

fname = "cross-section-profile-49.dat"
h = 0.15
t1 = 0.8
t2 = 3.0
N = 49
L = 1.0;

pi = 3.14159265359;
dx = L/(float)(N)
x = np.linspace(0.0+dx/2.0,L-dx/2.0,num=N)
S = np.ones(N) - h*np.power(np.sin(pi*np.power(x, t1)), t2)

np.savetxt(fname,S)
