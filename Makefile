#=========================================================================
include sfmake.inc
#=========================================================================
FC=$(SFMPI)/mpif90
DIREXE=$(HOME)/.bin
DIR=drivers
#EXE=ahm_real_trap
#EXE=ahm_real_disorder
EXE=ahm_matsubara_disorder

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

OBJS = SOLVER_VARS_GLOBAL.o IPT_SC_MATS.o  IPT_SC_SOPT.o IPT_MATS.o IPT_SOPT.o  SOLVER_INTERFACE.o RDMFT_VARS_GLOBAL.o RDMFT_FUNX.o RDMFT.o

ARGS= $(SFLIBS)
#ARGS=$(SFLIBS_DEB)
FLAG=$(STD)
#FLAG=$(OPT)
#FLAG=$(DEB)

#=================STANDARD COMPILATION====================================
compile: version $(OBJS)
	@echo " ........... compile: optimized ........... "
	@echo $(VER)
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)


.f90.o:	
	$(FC) $(FLAG) -c $< $(SFINCLUDE) 


clean: 
	@echo 'removing *.mod *.o *~'
	@rm -fv *.mod *.o *~ revision.inc

version:
	@echo $(VER)
