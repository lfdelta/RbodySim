#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import argparse

# parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('file', metavar='f', type=str, help="File prefix to analyze (.csv and .sys)")
parser.add_argument('-d', metavar='displaymode', type=int, nargs=1, default=[2])
arg = parser.parse_args()


# parse static data
sysdata = np.loadtxt(arg.file+".sys", dtype=str, delimiter=",") # global variables
sys = {}
for i in range(len(sysdata[0])):
  sys[sysdata[0][i]] = float(sysdata[1][i])
nsteps = int(sys['nsteps'])
tstep = sys['tstep']

dyndata = np.loadtxt(arg.file+".dyn", dtype=float, delimiter=",") # dynamic inertias
colsPerDyn = 2
masses = []; inertias = []
for i in range(0, len(dyndata), colsPerDyn):
  masses.append(dyndata[i])
  inertias.append(dyndata[i+1])

statdata = np.loadtxt(arg.file+".sta", dtype=float, delimiter=",") # static vertices
staticXs = []; staticYs = []
for i in range(0, len(statdata), 2):
  xs = (statdata[i]).tolist()
  ys = (statdata[i+1]).tolist()
  xs = np.array(xs + [xs[0]])
  ys = np.array(ys + [ys[0]])

  staticXs += [xs] # wraparound to close the polygon
  staticYs += [ys]


# parse simulated data
colsPerObj = 6 # x, y, vx, vy, angle, rotspeed
data = np.loadtxt(arg.file+".sim", dtype=float, delimiter=",", unpack=True)
nobjs = len(data)//colsPerObj




markers = ['o', 's', 'v', 'D', '*', '^', '.', 'P', 'x']

if (0 in arg.d):
  arg.d = [x for x in range(10)]

if (1 in arg.d): # plot x and y separately against t
  fig, ax = plt.subplots()
  for n in range(nobjs):
    a = 1 # alpha
    k1=0.5; k2=1.5 # min/max color value
    k = k1 + (k2-k1) * n/max(nobjs-1, 1)
    on = min(k, 1)
    off = max(0, k-1)
    r = [on,off,off,a]; g = [off,on,off,a]; b = [off,off,on,a]
    doLabel = (n == nobjs//2)

    ax.plot(data[n*colsPerObj], c=r, label=f"$x$" if doLabel else None) # x position
    ax.plot(data[n*colsPerObj + 1], c=g, label=f"$y$" if doLabel else None) # y position

  ax.set_title(f"System State Evolution, $dt = {tstep}$")
  ax.set_xlabel("Timestep $t_i$")
  ax.set_ylabel("Position")
  ax.legend()
  plt.show()


if (2 in arg.d): # 2D position plot
  fig, ax = plt.subplots()
  fig.set_tight_layout(True)

  for s in range(len(staticXs)):
    ax.plot(staticXs[s], staticYs[s], 'k-')

  for n in range(nobjs):
    xs = data[n*colsPerObj]
    ys = data[n*colsPerObj + 1]
    timeplot = ax.scatter(xs, ys, marker=markers[n%len(markers)], edgecolors=(0,0,0,0.5), c=range(nsteps+1), cmap=cmap.Reds, alpha=0.5)

  cbar = fig.colorbar(timeplot, ticks=[0, nsteps])
  cbar.ax.set_yticklabels(['$i=0$', f'$i={nsteps}$'])
  cbar.ax.set_ylabel("Timestep $t_i$", rotation=-90)
  ax.set_title(f"System State Evolution, $dt = {tstep}$")
  ax.set_xlabel("$x$ Position")
  ax.set_ylabel("$y$ Position")
  plt.show()


if (3 in arg.d): # conservation of energy and momentum
  fig, axs = plt.subplots(1, 3)
  fig.set_tight_layout(True)
  sysE = np.zeros(nsteps+1)
  sysPx = np.zeros(nsteps+1)
  sysPy = np.zeros(nsteps+1)

  for n in range(nobjs):
    m = masses[n]
    I = inertias[n]
    xs = data[n*colsPerObj]
    ys = data[n*colsPerObj + 1]
    xvels = data[n*colsPerObj + 2]
    yvels = data[n*colsPerObj + 3]
    rots = data[n*colsPerObj + 4]
    rvel = data[n*colsPerObj + 5]

    px = m * xvels
    py = m * yvels

    KE = 0.5 * (m*(xvels*xvels + yvels*yvels) + I*rvel*rvel) # mv^2 + Iω^2
    U = -m * ( sys['g_y']*(ys - ys.min()) + sys['g_x']*(xs - xs.min()) ) # mgh
    E = KE + U

    sysE += E
    sysPx += px
    sysPy += py

    a = 0.5 # alpha
    k1=0.4; k2=1.6 # min/max color value
    k = k1 + (k2-k1) * n/max(nobjs-1, 1)
    on = min(k, 1)
    off = max(0, k-1)
    r = [on,off,off,a]; g = [off,on,off,a]; b = [off,off,on,a]

    doLabel = (n == nobjs//2)
    axs[0].plot(xvels, c=r, label="$v_x$" if doLabel else None)
    axs[0].plot(yvels, c=g, label="$v_y$" if doLabel else None)

    axs[1].plot(KE, c=r, label="Kinetic $T$" if doLabel else None)
    axs[1].plot(U, c=g, label="Potential $U$" if doLabel else None)
    axs[1].plot(E, c=b, label="$E=T+U$" if doLabel else None)

  axs[2].plot(sysE, label="Total $E$")
  axs[2].plot(sysPx, label="Total $P_x$")
  axs[2].plot(sysPy, label="Total $P_y$")

  axs[0].set_title("Object Momenta")
  axs[0].legend()
  axs[1].set_title("Object Energies")
  axs[1].set_xlabel(f"Timestep $t_i$, $dt = {tstep}$")
  axs[1].set_ylabel("System Energy")
  axs[1].legend()
  axs[2].set_title("System-Wide Conservation")
  axs[2].legend()
  plt.show()
