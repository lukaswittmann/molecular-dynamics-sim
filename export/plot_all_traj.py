import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import os  
  
# read text file into pandas DataFrame
#r1 = pd.read_csv("r_1.txt", sep="        ")
#r2 = pd.read_csv("r_2.txt", sep="        ")

#pd.set_option('display.max_columns', None)


  
mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')


 
for x in os.listdir():
    if x.endswith(".txt"):
        r1 = np.loadtxt(x, delimiter='\t')
        ax.scatter(r1[:,0],r1[:,1],r1[:,2])

ax.set_xlim3d(0,10)
ax.set_ylim3d(0,10)
ax.set_zlim3d(0,10)

plt.show()
