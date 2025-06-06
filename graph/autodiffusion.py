import sys
import numpy as np
run_dir = sys.argv[1]
filepath = open('./'+run_dir+'/autodiffusion.dat','r')
data = np.loadtxt(filepath,skiprows=1)

tv = data[:,0]
vx = data[:,1]
vy = data[:,2]
vz = data[:,3]

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1,figsize=(8,5))
ax.plot(tv,1./6.*(vx+vy+vz))
ax.set_xlabel('t')
ax.set_ylabel('1/6 var(r^2)')
fig.savefig('./'+run_dir+'/autodiffusion.pdf')
