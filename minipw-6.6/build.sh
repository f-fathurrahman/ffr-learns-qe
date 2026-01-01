#!/bin/bash

#cd srclib
#make
#cd ../

INC="-I./srclib"
QE_HOME="/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib"

LIB="-Wl,-rpath ${QE_HOME} -L${QE_HOME} -lqemain -lblas -llapack"

LIBXC_HOME="/home/efefer/mysoftwares/libxc-4.3.4dyn"
LIB_LIBXC="-Wl,-rpath ${LIBXC_HOME}/lib -L${LIBXC_HOME}/lib -lxcf03 -lxc"

basn=`basename $1 .f90`
rm -vf $basn.x

# FIXME: rpath for Libxc is not

# -fbounds-check
# -fsanitize=address (for checking leaked memory)
gfortran -Wall $INC $1 $LIB_LIBXC $LIB -o $basn.x
echo "Executable: $basn.x"
