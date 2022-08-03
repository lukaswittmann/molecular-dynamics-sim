from pydoc import doc
from turtle import shape
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os  

nparticles = 50
r = np.loadtxt("r_1.txt", delimiter='\t')
frames = int(len(r))

r_array = np.zeros([nparticles,frames,3],dtype=float)
i = 0

for x in os.listdir():
    if x.endswith(".txt") and x.startswith("r_"):
        r = np.loadtxt(x, delimiter='\t')
        r_array[i,:,:] = r
        i += 1

def update(t):
    ax.cla()

    for i in range(0,nparticles):
        ax.scatter(r_array[i,t,0], r_array[i,t,1], r_array[i,t,2], s = 30, marker = 'o',color="tab:blue")

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_zlim(0, 10)

fig = plt.figure(dpi=100)
ax = fig.add_subplot(projection='3d')

ani = FuncAnimation(fig = fig, func = update, frames=np.arange(0, frames, 10), repeat=True, interval = 0)
#ani.save('md.mp4', dpi=150, writer='ffmpeg', fps=60)
plt.show()