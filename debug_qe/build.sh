#!/bin/bash

basn=`basename $1 .f90`
rm -vf $basn.x

INC="-I./Modules -I./UtilXlib"
LIB="Modules/libqemod.a UtilXlib/libutil.a"

# -fbounds-check
# -fsanitize=address (for checking leaked memory)
mpif90 -Wall $INC $1 $LIB -o $basn.x
echo "Executable: $basn.x"

