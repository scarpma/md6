
# coding: utf-8

# In[12]:


f = open('./data/stat.dat','r')
tv = []
px = []
py = []
pz = []
ek = []
ep = []
et = []
ii = 0
for line in f:
    line = line.split(' ')
    if  ii>0:
        tv.append(float(line[0]))
        px.append(float(line[1]))
        py.append(float(line[2]))
        pz.append(float(line[3]))
        ek.append(float(line[4]))
        ep.append(float(line[5]))
        et.append(float(line[6]))
    ii += 1


# In[13]:


import matplotlib.pyplot as plt

fig, axs = plt.subplots(1,2,figsize=(10,5))
axs[0].plot(tv,px,label='px')
axs[0].plot(tv,py,label='py')
axs[0].plot(tv,pz,label='pz')
axs[1].plot(tv,ek,label='ek')
axs[1].plot(tv,ep,label='ep')
axs[1].plot(tv,et,label='et')

axs[0].legend()
axs[1].legend()
axs[0].set_xlabel('t')
axs[1].set_xlabel('t')
fig.savefig('./data/stat.pdf')


# In[7]:


px

