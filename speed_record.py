#!/usr/local/bin/python3

# runs a suite of simulations and prints the total simulation time, along with other CSV-formatted data, to stdout

import os

ts = [1.5, 0.8, 0.1, 0.05, 0.01]

infloops = [(3,0.1), (6,0.1), (2,0.05), (4,0.05), (4,0.01)]
infloops += [(7,x) for x in ts]
infloops += [(8,x) for x in ts if x != 0.01]

print("continuous,scene,tstep,runtime")

for cont in [True, False]:
  for tstep in ts:
    nsteps = int(10 / tstep)
    for scene in range(1,9):
      if (cont) and (scene, tstep) in infloops:
        print("1,{},{},-1".format(scene, tstep))
      else:
        testsimargs = "-s {} -t {}{} -n {}".format(scene, tstep, " -c" if cont else "", nsteps)
        runtime = os.popen("./testsim " + testsimargs).read()
        print("{},{},{},{}".format(int(cont), scene, tstep, float(runtime)))