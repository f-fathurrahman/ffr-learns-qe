MAKE = make
AR = ar

F77 = gfortran
F77_OPTS = -Wall -O3 -cpp

CC = gcc
CC_OPTS = -std=c99

F90 = gfortran
F90_OPTS = -Wall -O3 -cpp -ffree-line-length-none  -I /usr/include
# -g -fbacktrace -ffpe-trap=zero,overflow,underflow

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIBS       = $(LIB_LINALG) $(LIB_FFTW3)
