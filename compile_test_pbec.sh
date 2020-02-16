
#LIBXC_HOME=/home/efefer/mysoftwares/libxc-3.0.0/

#LIBXC_HOME=/home/efefer/mysoftwares/libxc-4.0.4/

LIBXC_HOME=/home/efefer/mysoftwares/libxc-4.3.4/

#LIBXC_HOME2=/home/efefer/.julia/packages/Libxc/ArRdo/deps/usr

mpif90 -I$LIBXC_HOME/include/ test_pbec.f90 qe/Modules/libqemod.a qe/UtilXlib/libutil.a \
-L$LIBXC_HOME/lib/ -lxcf90 -lxc -o test_pbec.x
