# RbodySim

## Setup
Requires [Eigen 3.3.5](http://eigen.tuxfamily.org/index.php) headers to be in the local "include" path.
Requires Python3 in order to run analysis scripts.

## Instructions
1. (First time only) Run `make` while inside the root directory. This will produce an executable file (`testsim`) in the same directory.
1. To compute and analyze a full simulation, call `./simulate.py [-c] [-s test_scenario] [-t tstep] [-n nsteps]`
  * `-c` tells the program to run with continuous simulation (default is discrete)
  * `test_scenario` is an integer between 1 and 8, corresponding to a predefined setup to simulate.
  * `tstep` is the time (in generic units) which each timestep will integrate over. Usually 0.01 to 0.1 is reasonable.
  * `nsteps` is the number of timesteps to simulate. Usually 100 to 1000 is reasonable.
1. You can also run a raw simulation using `./testsim [-c] [-s test_scenario] [-t tstep] [-n nsteps] [-f outfile_prefix]`
  * `outfile_prefix` is the identifier for a set of files which will store information for later analysis
1. The analysis can also be performed raw, via `./analyze.py [-d display] [-p img_prefix] infile_prefix`
  * `infile_prefix` must correspond to the `outfile_prefix` produced by `testsim`
  * `display` is an integer (1, 2, 3) corresponding to a particular visualization. Pass 0 or none to view all three.
  * `img_prefix` is the identifier for a set of visualization output files. If none is provided, visualizations are rendered to the screen.

## Samples
Sample visualization outputs are provided in the `data` directory. This is where data and image outputs are stored when `simulate.py` is called.