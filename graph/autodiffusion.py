f = open('./data/autodiffusion.dat','r')
tv = []
vx = []
vy = []
vz = []
ii = 0
for line in f:
    line = line.split(' ')
    if ii>0:
        tv.append(float(line[0]))
        vx.append(float(line[1]))
        vy.append(float(line[2]))
        vz.append(float(line[3]))
    ii += 1


import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1,figsize=(8,5))
ax.plot(tv,vx,label='var x')
ax.plot(tv,vy,label='var y')
ax.plot(tv,vz,label='var z')
ax.legend()
ax.set_xlabel('t')
fig.savefig('./data/autodiffusion.pdf')
