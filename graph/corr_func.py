f = open('./data/corr_func.dat','r')
xx = []
cf = []
pp = []
ii = 0
for line in f:
    line = line.split(' ')
    if ii>0:
        xx.append(float(line[0]))
        cf.append(float(line[1]))
        pp.append(float(line[2]))
    ii += 1

import matplotlib.pyplot as plt
fig, axs = plt.subplots(1,2,figsize=(15,5))
axs[0].plot(xx,cf,label='pair correlation function')
axs[0].plot(xx,pp,label='4pi rˆ2 delta r')
axs[0].legend()
axs[0].set_xlabel('r')
axs[1].plot(xx,[cf[i]/pp[i] for i in range(len(pp))],label='pair correlation function density')
axs[1].set_xlabel('r')
axs[1].legend()
fig.savefig('./salvati/corr_func.pdf')
