#!/usr/local/bin/python3

import sys
import os
import argparse

# extract info from the argument string
# parser = argparse.ArgumentParser()
# parser.add_argument('-s', type=int, nargs=1)
# parser.add_argument('-c', action='store_true')
# argflag = parser.parse_args(sys.argv[1:])
argS = 0
argC = False
file = False

for i in range(1, len(sys.argv)-1):
  if sys.argv[i] == "-s":
    argS = int(sys.argv[i+1])
  elif sys.argv[i] == "-c":
    argC = True

# fill in the blanks to pass in filename info
if (argS > 0):
  prefixarg = "data/{}{}".format("c" if argC else "d", argS)
  file = True
else:
  prefixarg = ""

mainargs = " ".join(sys.argv[1:])
testsimargs = " ".join([mainargs, "-f", prefixarg]) if file else mainargs
analyzeargs = " ".join(["-p", prefixarg, prefixarg])

os.system("./testsim " + testsimargs)
if file:
  os.system("./analyze.py " + analyzeargs)