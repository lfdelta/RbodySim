#!/usr/local/bin/python3

import sys
import os
import argparse

# extract info from the argument string
argS = 0
argC = False
file = False

for i in range(1, len(sys.argv)):
  if sys.argv[i] == "-s":
    argS = int(sys.argv[i+1])
  elif sys.argv[i] == "-c":
    argC = True

# fill in the blanks to pass in filename info
if (argS > 0):
  prefixarg = "data/{}{}".format(argS, "c" if argC else "d")
  file = True
else:
  prefixarg = ""

mainargs = " ".join(sys.argv[1:])
testsimargs = " ".join([mainargs, "-f", prefixarg]) if file else mainargs
analyzeargs = " ".join(["-p", prefixarg, prefixarg])

os.system("./testsim " + testsimargs)
if file:
  os.system("./analyze.py " + analyzeargs)