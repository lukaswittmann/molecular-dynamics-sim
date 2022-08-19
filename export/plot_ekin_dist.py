from cProfile import label
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import os  
  
velo = np.loadtxt("Ekin.txt")

bins = 50

plt.grid()
plt.hist(velo[-1,:], bins = bins)
plt.xlabel("Ekin [J]")
plt.ylabel("N")
plt.tight_layout()  
plt.show()
