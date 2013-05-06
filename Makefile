#=========================================================================
include make.inc
#=========================================================================
#FC=$(SFMPI)/mpif90
FC=mpif90
DIREXE=$(HOME)/.bin
DIR=drivers
EXE=ahm_matsubara_trap
#EXE=ahm_matsubara_disorder
#EXE=ahm_real_trap
#EXE=ahm_real_disorder

.SUFFIXES: .f90 
OBJS = SOLVER_VARS_GLOBAL.o IPT_SC_MATS.o  IPT_SC_SOPT.o IPT_MATS.o IPT_SOPT.o  SOLVER_INTERFACE.o RDMFT_VARS_GLOBAL.o

ARGS= $(SFMODS) $(SFLIBS)
ARGS_DEB=$(SFMODS_DEB) $(SFLIBS_DEB)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
REV = $(shell git rev-parse HEAD)
VER='character(len=41),parameter :: revision = "$(REV)"' > revision.inc

#=================STANDARD COMPILATION====================================
all: FLAG=$(STD)
all: version $(OBJS)
	@echo " ........... compile: optimized ........... "
	@echo $(VER)
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

#================OPTIMIZED COMPILATION====================================
opt: FLAG=$(OPT)
opt: 	version $(OBJS)
	@echo " ........... compile: optimized ........... "
	@echo $(VER)
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)


#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB)
debug: 	version $(OBJS)
	@echo " ........... compile : debug   ........... "
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS_DEB)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)


pade: 	version $(OBJS)
	@echo " ........... compile: rdmft_pade ........... "
	$(FC) $(STD) $(OBJS) drivers/pade_matsubara_to_real.f90 -o $(DIREXE)/pade_matsubara_to_real $(LIBDMFT) $(SFLIBS) $(SFMODS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/pade_matsubara_to_real

data: 	version $(OBJS)
	@echo " ........... compile: get_data ........... "
	$(FC) $(STD) $(OBJS) $(DIR)/get_data_$(EXE).f90 -o $(DIREXE)/get_data_$(EXE)_$(BRANCH) $(LIBDMFT) $(SFLIBS) $(SFMODS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/get_data_$(EXE)


lib: 
	@make -C ./library_src


.f90.o:	
	$(FC) $(FLAG) -c $< $(SFMODS) 


clean: 
	@echo 'removing *.mod *.o *~'
	@rm -fv *.mod *.o *~ revision.inc

clean_all: clean
	@echo 'cleaning the library'
	@make -C ./library_src clean
version:
	@echo $(VER)
