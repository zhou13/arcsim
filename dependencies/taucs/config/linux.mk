#########################################################
# Linux                                                 #
#########################################################
OBJEXT=.o
LIBEXT=.a
EXEEXT= 
F2CEXT=.f
PATHSEP=/
DEFFLG=-D

FC        = gfortran
FFLAGS    = -O3 -g -fno-second-underscore -Wall 
FOUTFLG   =-o 

COUTFLG   = -o
CFLAGS    = -O3 -g -D_POSIX_C_SOURCE=199506L -Wall -pedantic -ansi -fPIC -fexceptions -D_GNU_SOURCE 
CFLAGS    = -g -O3 -Wall -Werror -pedantic -ansi 
# for some reason, -std=c99 -pedantic crashes 
# with the error message "imaginary constants are a GCC extension"
# (seems to be a gcc bug, gcc 3.3.1)
CFLAGS    = -O3 -Wall -Werror -std=c89 -pedantic
CFLAGS    = -O3 -Wall -Werror -std=c99 
CFLAGS    = -O3 -Wall -std=c99

LD        = $(CC) 
LDFLAGS   = 
LOUTFLG   = $(COUTFLG)

AR        = ar cr
AOUTFLG   =

RANLIB    = ranlib
RM        = rm -rf

LIBBLAS   = -lblas
LIBLAPACK = -llapack

LIBMETIS  = -lmetis

LIBF77 = -lgfortran
LIBC   = -lm 
