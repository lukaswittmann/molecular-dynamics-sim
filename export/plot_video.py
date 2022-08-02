import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
import pandas as pd


r1 = np.loadtxt("r_1.txt", delimiter='\t')
#r2 = np.loadtxt("r_2.txt", delimiter='\t')
#r3 = np.loadtxt("r_3.txt", delimiter='\t')
#r4 = np.loadtxt("r_4.txt", delimiter='\t')
#r5 = np.loadtxt("r_5.txt", delimiter='\t')
t = np.arange(0,len(r1))

df = pd.DataFrame({"time": t ,"x" : r1[:,0], "y" : r1[:,1], "z" : r1[:,2]})



def update_graph(num):
    data=df[df['time']==num]
    graph.set_data (data.x, data.y)
    graph.set_3d_properties(data.z)
    title.set_text('3D Test, time={}'.format(num))
    return title, graph


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d(0,10)
ax.set_ylim3d(0,10)
ax.set_zlim3d(0,10)
title = ax.set_title('3D Test')

data=df[df['time']==0]
graph, = ax.plot(data.x, data.y, data.z, linestyle="", marker="o")

ani = matplotlib.animation.FuncAnimation(fig, update_graph, len(r1), interval=0, blit=True)
plt.show()

