#!/bin/sh
#
# installs the necessary dependencies for building with GNU (`gfortran`) and OpenMPI
#
sudo apt-get install -y --no-install-recommends gfortran libopenmpi-dev libfftw3-dev libblas-dev liblapack-dev
