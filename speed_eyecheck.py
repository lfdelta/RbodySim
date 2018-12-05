#!/usr/local/bin/python3

# produces an animation-like visualization of the system through time, for checking tunneling and other anomolous behavior by eye

import os

for cont in [True, False]:
  for tstep in [1.5, 0.8, 0.1, 0.05, 0.01]:
    nsteps = int(10 / tstep)
    for scene in range(1,9):
      testsimargs = "-s {} -t {}{} -n {} -f data/tmp".format(scene, tstep, " -c" if cont else "", nsteps)
      os.system("./testsim " + testsimargs)
      analyzeargs = "-p {}{} -d 2 data/tmp".format("c" if cont else "d", scene)
      os.system("./analyze.py " + analyzeargs)