import matplotlib.pyplot as plt
import numpy as np


offsets = {
  "DB": 28.28,
  "DG": 27.16,
  "DL": 27.35,
  "DT": 27.05,
  "DX": 28.77,
  "DY": 26.57,
  "PO": 26.5
}
mem = ["up","lw"]
species = ["DT","DL","DB","DX"]
fig = plt.figure(figsize = (5,5))

for spec in species:
  for leafs in mem:
    data = np.genfromtxt((spec+'.'+leafs+'.height.dat'),missing_values='np.nan')
    height = data[:,2:] + offsets[spec]

    x = np.arange(2,71,4)
    y_vals_lower = np.zeros(18)



    for r in range(height.shape[0]):
      tot = 0
      for y in height[r]:
        tot = tot + y 
      tot = tot / 30
      y_vals_lower[r] = tot
    plt.plot(x,y_vals_lower)

plt.show()

