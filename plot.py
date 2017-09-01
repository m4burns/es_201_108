"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import threading
import sys

fig, ax = plt.subplots()

x = np.zeros([100,15])
xl = threading.Lock()
image = ax.imshow(x.T, clim=(0,10))

def push(line):
  global x
  with xl:
    x = np.concatenate((x, [line]), axis=0)
    x = np.delete(x, (0), axis=0)

def animate(i):
  with xl:
    image.set_data(x.T)  # update the data
  return image,

def init():
  with xl:
    image.set_data(x.T)
  return image,

def render():
  ani = animation.FuncAnimation(fig, animate, np.arange(1, 200), init_func=init,
                                interval=100, blit=True)
  plt.show()
  sys.exit(0)
