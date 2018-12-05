#!/usr/local/bin/python3

# plots the results from the speed analysis, which must be passed in as an argument

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as clr
import sys

continuous, scene, tstep, runtime = np.loadtxt(sys.argv[1], skiprows=1, delimiter=',', unpack=True)
nsteps = np.array(10 / tstep, dtype=int)
scene = np.array(scene, dtype=int)

tmap = {0.01:0, 0.05:1, 0.1:2, 0.8:3, 1.5:4}

ctimes = np.ndarray((8, 5))
dtimes = np.ndarray((8, 5))
negclr = (0,0.8,0.3,1)
tunnelclr = (0,0.3,0.8,1)

neg = cm.Reds_r
neg.set_bad(color=negclr)

for i in range(len(scene)//2):
  ctimes[scene[i]-1, tmap[tstep[i]]] = runtime[i]/nsteps[i]
for i in range(len(scene)//2, len(scene)):
  dtimes[scene[i]-1, tmap[tstep[i]]] = runtime[i]/nsteps[i]

neg_mask = np.ma.masked_where(ctimes < 0, ctimes)

ctunnels = [(6,x) for x in tmap if x != 0.1]
dtunnels = [(3,0.1), (3,0.8), (7, 0.8)]
dtunnels += [(x,1.5) for x in range(1,9)]
dtunnels += [(6,x) for x in tmap]

tunnel = clr.LinearSegmentedColormap.from_list('tunnel', [(0,0,0,0), (0,0,0,0)], N=2)
tunnel.set_bad(color=tunnelclr)
cmasktmp = ctimes.copy()
for (i,j) in ctunnels:
  cmasktmp[i-1][tmap[j]] = -2
ctunnel_mask = np.ma.masked_where(cmasktmp == -2, cmasktmp)

dmasktmp = dtimes.copy()
for (i,j) in dtunnels:
  dmasktmp[i-1][tmap[j]] = -2
dtunnel_mask = np.ma.masked_where(dmasktmp == -2, dmasktmp)

axsz = 15
headersz = 20
tunnel_patch = mpatches.Patch(color=tunnelclr, label='Tunneling')
neg_patch = mpatches.Patch(color=negclr, label='Infinite\nLoop')


fig, ax = plt.subplots()
fig.set_tight_layout(True)
cplot = ax.imshow(neg_mask, cmap=neg, vmin=0, extent=(-0.5,4.5, 8.5,0.5), aspect='equal')
ax.imshow(ctunnel_mask, cmap=tunnel, extent=(-0.5,4.5, 8.5,0.5), aspect='equal')
ax.set_ylabel("Test Scenario", size=axsz)
ax.set_xlabel("Timestep index", size=axsz)
cpbar = fig.colorbar(cplot, ax=ax)
cpbar.ax.set_ylabel("Real time (s) per timestep", size=axsz)
ax.set_title("Continuous Collision\nResponse", size=headersz)
fig.legend(handles=[tunnel_patch, neg_patch], loc="center left", prop={'size' : axsz})
plt.show()

fig, ax = plt.subplots()
fig.set_tight_layout(True)
dplot = ax.imshow(dtimes, cmap=neg, extent=(-0.5,4.5, 8.5,0.5), aspect='equal')
ax.imshow(dtunnel_mask, cmap=tunnel, extent=(-0.5,4.5, 8.5,0.5), aspect='equal')
ax.set_ylabel("Test Scenario", size=axsz)
ax.set_xlabel("Timestep index", size=axsz)
dpbar = fig.colorbar(dplot, ax=ax)
dpbar.ax.set_ylabel("Real time (s) per timestep", size=axsz)
ax.set_title("Discrete Collision\nResponse", size=headersz)
fig.legend(handles=[tunnel_patch, neg_patch], loc="center left", prop={'size' : axsz})

plt.show()