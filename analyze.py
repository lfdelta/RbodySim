#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import argparse

# parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('file', metavar='f', type=str, help="File prefix to analyze (.csv and .sys)")
parser.add_argument('-d', metavar='displaymode', type=int, nargs=1, default=[2], choices=(1,2))
arg = parser.parse_args()


# parsing static data
sysdata = np.loadtxt(arg.file+".sys", dtype=str, delimiter=",")
sys = {}
for i in range(len(sysdata[0])):
  sys[sysdata[0][i]] = sysdata[1][i]
nsteps = int(sys['nsteps'])
tstep = float(sys['tstep'])

sysdata = np.loadtxt(arg.file+".sys", dtype=float, delimiter=",", skiprows=2)
staticXs = []; staticYs = []
for i in range(0, len(sysdata), 2):
  xs = (sysdata[i]).tolist()
  ys = (sysdata[i+1]).tolist()
  xs = np.array(xs + [xs[0]])
  ys = np.array(ys + [ys[0]])

  staticXs += [xs] # wraparound to close the polygon
  staticYs += [ys]


# parsing simulated data
colsPerObj = 4
data = np.loadtxt(arg.file+".csv", dtype=float, delimiter=",", unpack=True)
nobjs = len(data)//colsPerObj


if (arg.d[0] == 1): # plot x and y separately against t
  fig, ax = plt.subplots()
  for n in range(nobjs):
    ax.plot(data[n*colsPerObj], label=f"$x_{{{n+1}}}$") # x position
    ax.plot(data[n*colsPerObj + 1], label=f"$y_{{{n+1}}}$") # y position

  ax.set_title(f"System State Evolution, $dt = {tstep}$")
  ax.set_xlabel("Timestep $t_i$")
  ax.set_ylabel("Position")
  ax.legend()


if (arg.d[0] == 2): # 2D plot
  markers = ['o', 's', 'v', 'D', '*', '^', '.', 'P', 'x']
  fig, ax = plt.subplots()
  fig.set_tight_layout(True)

  print(type(staticXs[0]))
  for s in range(len(staticXs)):
    ax.plot(staticXs[s], staticYs[s], 'k-')

  for n in range(nobjs):
    xs = data[n*colsPerObj]
    ys = data[n*colsPerObj + 1]
    timeplot = ax.scatter(xs, ys, marker=markers[n%len(markers)], edgecolors=(0,0,0,0.5), c=range(nsteps+1), cmap=cmap.Reds, label=f"Object {n+1}")

  cbar = fig.colorbar(timeplot, ticks=[0, nsteps])
  cbar.ax.set_yticklabels(['$i=0$', f'$i={nsteps}$'])
  cbar.ax.set_ylabel("Timestep $t_i$", rotation=-90)
  ax.set_title(f"System State Evolution, $dt = {tstep}$")
  ax.set_xlabel("$x$ Position")
  ax.set_ylabel("$y$ Position")
  ax.legend()

plt.show()