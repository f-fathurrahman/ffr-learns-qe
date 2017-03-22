QE_DIR = /home/efefer/WORKS/QE/qe-6.0

# import several variables from here
include $(QE_DIR)/make.inc


# The variable TOPDIR in make.inc of QE is not defined properly
# for our purpose, so we define again some important variables here

INC_QEMOD = $(QE_DIR)/Modules
INC_PW    = $(QE_DIR)/PW/src
INC_QEFFT = $(QE_DIR)/FFTXlib
INC_PP    = $(QE_DIR)/PP/src

INC_ALL   = -I$(INC_QEMOD) -I$(INC_PW) -I$(INC_QEFFT) -I$(INC_PP)


LIB_LINALG = $(BLAS_LIBS) $(LAPACK_LIBS)
LIB_FFT    = $(FFT_LIBS)

LIB_PW    = $(QE_DIR)/PW/src/libpw.a
LIB_PP    = $(QE_DIR)/PP/src/libpp.a
LIB_QEMOD = $(QE_DIR)/Modules/libqemod.a
LIB_QEFFT = $(QE_DIR)/FFTXlib/libqefft.a
LIB_QELA  = $(QE_DIR)/LAXlib/libqela.a
LIB_CLIB  = $(QE_DIR)/clib/clib.a
LIB_IOTK  = $(QE_DIR)/iotk/src/libiotk.a
LIB_QE    = $(LIB_PW) $(LIB_PP) $(LIB_QEMOD) $(LIB_QEFFT) $(LIB_QELA) \
            $(LIB_CLIB) $(LIB_IOTK)

# This is required if F90 is not defined (or commented out) in make.inc
# file of QE
F90 = $(MPIF90)

F90_OPTS = -Wall
# F90FLAGS is used when compiling with macros defined in QE installation
# F90_OPTS is used for other purpose

