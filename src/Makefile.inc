# C and C++ compilers
CC = gcc
CXX = g++

# C and C++ compilation flags
# -DENABLE_S4_TRACE  -DENABLE_RCWA_TRACE -DDUMP_MATRICES
FLAGS = -m32 -O0 -ggdb -DUSE_LAPACK -DUSE_BLAS -DHAVE_PTHREADS -DHAVE_CHOLMOD -DHAVE_FFTW3 -Wall -Wextra
CFLAGS = $(FLAGS)
CXXFLAGS = $(FLAGS)
LDFLAGS = -m32 -O0 -ggdb

# Lua directories:
LUA_INC = -I/c/Program\ Files\ \(x86\)/Lua/5.1/include
LUA_LIB = /d/dev/libs/output/lua-5.1.4/i386/liblua.a

# PThread directories:
PTHREAD_INC = -DPTW32_STATIC_LIB
PTHREAD_LIB = -lpthreadGC2-static -lwsock32

# LAPACK directories:
LA_INC =
LAPACK_LIBS = /d/dev/libs/output/lapack-3.2.2-mingw32-vanilla/lapack_mingw32.a /d/dev/libs/output/lapack-3.2.2-mingw32-vanilla/blas_mingw32.a -lgfortran -static
#LAPACK_LIBS = /d/dev/libs/GotoBLAS2/libgoto2_barcelona-r1.13.lib -lgfortran
#LAPACK_LIBS = -static /d/dev/libs/output/numerical_mingw32/libgoto2_barcelonap-r1.13.lib /d/dev/libs/output/numerical_mingw32/libgoto2_barcelona-r1.13.lib -lgfortran
#LAPACK_LIBS = /c/AMD/acml4.4.0/ifort32/lib/libacml.lib

FFTW3_INC = -I/d/dev/libs/output/fftw3
FFTW3_LIBS = -L/d/dev/libs/output/fftw3 -llibfftw3-3

# CHOLMOD directories
CHOLMOD_INC = -I/d/dev/libs/output/numerical_mingw32/cholmod
CHOLMOD_LIBS = -L/d/dev/libs/output/numerical_mingw32 -lcholmod -lamd -lcamd -lcolamd -lccolamd

#DEBUG_LIBS = /d/dev/libs/output/duma/libduma.a

# All needed libraries
LIBS = $(LUA_LIB) $(FFTW3_LIBS) $(CHOLMOD_LIBS) $(PTHREAD_LIB) $(LAPACK_LIBS) $(DEBUG_LIB)

