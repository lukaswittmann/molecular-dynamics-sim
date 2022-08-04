from cProfile import label
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np
import os  
  
time = np.loadtxt("time.txt")
P = np.loadtxt("P.txt")
V = np.loadtxt("V.txt")
T = np.loadtxt("T.txt")


t0 = 0
t1 = 1000

plt.rcParams['figure.figsize'] = [12, 12]
figure, axis = plt.subplots(3, 1)
  
axis[0].plot(time[t0:t1], T[t0:t1])
axis[0].set_title("Temperature")
  
axis[1].plot(time[t0:t1], P[t0:t1])
axis[1].set_title("Pressure")

axis[2].plot(time[t0:t1], V[t0:t1])
axis[2].set_title("Volume")
  
plt.tight_layout()
plt.show()