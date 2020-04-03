f = open('./data/corr_func.dat','r')
xx = []
cf = []
ii = 0
for line in f:
    line = line.split(' ')
    if ii>0:
        xx.append(float(line[0]))
        cf.append(float(line[1]))
    ii += 1

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1,figsize=(8,5))
ax.plot(xx,cf,label='pair_correlation')
ax.legend()
ax.set_xlabel('r')
fig.savefig('./data/corr_func.pdf')
