import matplotlib.pyplot as plt
import numpy as np

print("Testing ...")

data = np.loadtxt("stat.dat")
data_ref = np.loadtxt("ref_files/stat.dat")
diff = data - data_ref

n_cols = data.shape[1] - 1
for i in range(1, n_cols):
  print(f"col {i} max diff: {diff[:,i].max():.2g}")

print("===========================")
if diff.max()==0.:
  print("ALL TEST PASSED")
  print("YOU'RE GOOD TO GO")
else:
  print("SOME TEST DID NOT PASS")
print("===========================\n\n")


#for i in range(1, n_cols):
#  plt.plot(data[:,0], data[:,i], label='new')
#  plt.plot(data_ref[:,0], data_ref[:,i], label='ref')
#  plt.title(f'col {i}')
#  plt.show()
