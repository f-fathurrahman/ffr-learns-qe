#!/bin/bash

INC="-I./srclib"
QE_HOME="/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib"
LIB="-Wl,-rpath,${QE_HOME} -L${QE_HOME} -lqemain -lblas -llapack"
LIBXC_HOME="/home/efefer/mysoftwares/libxc-4.3.4dyn"
LIB_LIBXC="-Wl,-rpath,${LIBXC_HOME} -L${LIBXC_HOME}/lib -lxcf90 -lxcf03 -lxc"

basn=`basename $1 .f90`
rm -vf $basn.x

gfortran -Wall $INC $1 $LIB $LIB_LIBXC -o $basn.x
echo "Executable: $basn.x"
