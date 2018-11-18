#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import argparse

# parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('file', metavar='f', type=str, help="File prefix to analyze (.csv and .sys)")
parser.add_argument('-d', metavar='displaymode', type=int, nargs=1, default=[2], choices=(1,2,3))
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

if (arg.d[0] == 1): # plot x and y separately against t
  fig, ax = plt.subplots()
  for n in range(nobjs):
    ax.plot(data[n*colsPerObj], label=f"$x_{{{n+1}}}$") # x position
    ax.plot(data[n*colsPerObj + 1], label=f"$y_{{{n+1}}}$") # y position

  ax.set_title(f"System State Evolution, $dt = {tstep}$")
  ax.set_xlabel("Timestep $t_i$")
  ax.set_ylabel("Position")
  ax.legend()


if (arg.d[0] == 2): # 2D position plot
  fig, ax = plt.subplots()
  fig.set_tight_layout(True)

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


if (arg.d[0] == 3): # conservation of energy and momentum
  fig, axs = plt.subplots(1, 3)
  fig.set_tight_layout(True)
  sysE = np.ndarray(nsteps+1)

  for n in range(nobjs):
    m = masses[n]
    I = inertias[n]
    xs = data[n*colsPerObj]
    ys = data[n*colsPerObj + 1]
    xvels = data[n*colsPerObj + 2]
    yvels = data[n*colsPerObj + 3]
    rots = data[n*colsPerObj + 4]
    rvel = data[n*colsPerObj + 5]
    gh = ys - ys.min()

    KE = 0.5 * ((m * xvels*xvels + yvels*yvels) + I * rvel*rvel)
    U = m * gh
    E = KE+U
    sysE += E

    axs[0].plot(xvels, label=f"$v_{{x{n+1}}}$") # x velocity
    axs[0].plot(yvels, label=f"$v_{{y{n+1}}}$") # y velocity

    axs[1].plot(KE, label=f"Kinetic $T_{{{n+1}}}$")
    axs[1].plot(U, label=f"Potential $U_{{{n+1}}}$")
    axs[1].plot(E, alpha=0.5, label=f"Total $E_{{{n+1}}}$")

  axs[2].plot(sysE, label="Total $E$")

  axs[0].set_title("Object Momenta")
  axs[0].legend()
  axs[1].set_title("Object Energies")
  axs[1].set_xlabel(f"Timestep $t_i$, $dt = {tstep}$")
  axs[1].set_ylabel("System Energy")
  axs[1].legend()
  axs[2].set_title("System-Wide Conservation")
  axs[2].legend()

plt.show()