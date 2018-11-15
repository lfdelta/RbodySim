#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', metavar='f', type=str, help='CSV file to analyze')
arg = parser.parse_args()

xs, ys = np.loadtxt(arg.file, dtype=float, delimiter=" ", unpack=True)

plt.plot(ys)
plt.show()