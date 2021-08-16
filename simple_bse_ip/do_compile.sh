QEPATH=/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/
INCLUDES="-I$QEPATH/Modules -I$QEPATH/include "

#mpif90 -O3 -g -x f95-cpp-input -D__FFTW3 -D__MPI    -I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//include -I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//FoX/finclude -I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//S3DE/iotk/include/ -I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//iotk/src -I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//Modules -I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//FFTXlib -I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//LAXlib -I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//UtilXlib -I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//FoX/finclude -I../../PW/src -I../pw4gww -I../gww

echo $INCLUDES

mpif90 -c -O3 -g -x f95-cpp-input -D__FFTW3 -D__MPI \
-I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/include \
-I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/FoX/finclude \
-I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/S3DE/iotk/include \
-I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/iotk/src \
-I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/Modules \
-I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/FFTXlib \
-I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/LAXlib \
-I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/UtilXlib \
-I/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e/FoX/finclude \
-I$QEPATH/PW/src \
simple.f90

#mpif90 -g -o simple.x \
#simple.o libsimple.a ../../PW/src/libpw.a ../gww/libgww.a ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../KS_Solvers/libks_solvers.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a /home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//clib/clib.a  /home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//iotk/src/libiotk.a    -L/home/efefer/WORKS/QE/qe_simple_ip_bse/q-e//FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys  -lfftw3  -L/home/efefer/mysoftwares/LIBMKL -lmkl_gf_lp64  -lmkl_sequential -lmkl_core     ../minpack/minpacklib.a

mpif90 -g -o simple.x simple.o \
$QEPATH/PW/src/libpw.a \
$QEPATH/GWW/gww/libgww.a \
$QEPATH/Modules/libqemod.a \
$QEPATH/FFTXlib/libqefft.a \
$QEPATH/KS_Solvers/libks_solvers.a \
$QEPATH/LAXlib/libqela.a \
$QEPATH/UtilXlib/libutil.a \
$QEPATH/clib/clib.a \
$QEPATH/iotk/src/libiotk.a \
-L$QEPATH/FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys \
-lfftw3 \
-L/home/efefer/mysoftwares/LIBMKL -lmkl_gf_lp64  -lmkl_sequential -lmkl_core \
$QEPATH/GWW/minpack/minpacklib.a
