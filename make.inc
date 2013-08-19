#COMPILATION FLAGS:
FC=ifort
MOPT=-module 
STD=-O2 #more optimized options:-O3 -ftz -unroll -assume nobuffered_io 
DEB=-O0 -p -traceback -g -debug -fpe0 -check all

#FLAGS for the library
FFLAG+= $(STD) -static
DFLAG+= $(DEB) -static

##If you have MKL installed and sourced:
#IFORT 10: -lmkl_em64t -lmkl_core -liomp5  -lpthread
#IFORT 11,12:
MATHLIB = -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
#MATHLIB = -Wl,--start-group  $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a $(MKLROOT)/lib/em64t/libmkl_sequential.a $(MKLROOT)/lib/em64t/libmkl_core.a -Wl,--end-group -lpthread -lm

##Else comment previous line and comment out the following ones:
##In some systems LAPACK/BLAS are default, so there's no need to
##comment out their path directories (LAPACKDIR, BLASDIR)
#LAPACKDIR=  #e.g. /usr/lib or /opt/lapack or $(HOME)/lib/local/lapack sometimes this is not even necessary
#BLASDIR=    #e.g. /usr/lib or /opt/blas or $(HOME)/lib/local/blas sometimes this is not even necessary
#MATHLIB = -L$(LAPACKDIR) -L$(BLASDIR) -llapack -lblas

SFLIBS 	   = -L./library_src -lscifor ${MATHLIB}
SFMODS     = -I./library_src

SFLIBS_DEB = -L./library_src -lscifor_deb ${MATHLIB}
SFMODS_DEB = -I./library_src/debug

