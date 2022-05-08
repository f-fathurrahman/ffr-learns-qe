#!/bin/bash
basn=`basename $1 .f90`
rm -vf $basn.x

gfortran -I srclib -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  $1 -o  $basn.x srclib/libmain.a \
/home/efefer/mysoftwares/libxc-4.3.4/lib/libxcf03.a \
/home/efefer/mysoftwares/libxc-4.3.4/lib/libxc.a \
-lblas -llapack

echo "Executable: $basn.x"
