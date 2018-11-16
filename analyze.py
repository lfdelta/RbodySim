#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import argparse

colsPerObj = 2

# parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('file', metavar='f', type=str, help="CSV file to analyze")
parser.add_argument('-d', metavar='displaymode', type=int, nargs=1, default=[2], choices=(1,2))
arg = parser.parse_args()


# parsing data
data = np.loadtxt(arg.file, dtype=float, delimiter=",", unpack=True)
nobjs = len(data)//colsPerObj

if (arg.d[0] == 1): # plot x and y separately against t
  fig, ax = plt.subplots()
  for n in range(nobjs):
    ax.plot(data[n*colsPerObj], label=f"$x_{{{n+1}}}$") # x position
    ax.plot(data[n*colsPerObj + 1], label=f"$y_{{{n+1}}}$") # y position

  ax.set_title("System State Evolution")
  ax.set_xlabel("Timestep $t_i$")
  ax.set_ylabel("Position")
  ax.legend()


if (arg.d[0] == 2): # 2D plot
  markers = ['o', 's', 'v', 'D', '*', '^', '.', 'P', 'x']
  fig, ax = plt.subplots()
  fig.set_tight_layout(True)
  nsteps = len(data[0])
  for n in range(nobjs):
    xs = data[n*colsPerObj]
    ys = data[n*colsPerObj + 1]
    timeplot = ax.scatter(xs, ys, marker=markers[n%len(markers)], edgecolors=(0,0,0,0.5), c=range(nsteps), cmap=cmap.Reds, label=f"Object {n}")
  cbar = fig.colorbar(timeplot, ticks=[0, nsteps-1])
  cbar.ax.set_yticklabels(['$i=1$', f'$i={nsteps}$'])
  cbar.ax.set_ylabel("Timestep $t_i$", rotation=-90)
  ax.set_title("System State Evolution")
  ax.set_xlabel("$x$ Position")
  ax.set_ylabel("$y$ Position")
  ax.legend()

plt.show()