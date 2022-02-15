#!/bin/bash
gfortran -c -I srclib -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  $1
