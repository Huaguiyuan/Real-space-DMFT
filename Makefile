#=========================================================================
include lib.mk
#=========================================================================
FC=$(SFMPI)/mpif90
EXE   =ahm_matsubara_disorder
#EXE   =ahm_matsubara_disorder_phase
DIR   =drivers/ahm_matsubara
DIREXE=$(HOME)/.bin

STD+=-static
OPT+=-static
DEB+=-static

OBJS=RDMFT_VARS_GLOBAL.o
OBJS_DEB=RDMFT_VARS_GLOBAL_DEB.o
all: 	version $(OBJS)
	@echo " ........... compile: standard ........... "
	$(FC) $(STD) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(SFLIBS) $(SFMODS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)

opt: 	version $(OBJS)
	@echo " ........... compile: optimized ........... "
	$(FC) $(OPT) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(SFLIBS) $(SFMODS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)


debug: 	version $(OBJS_DEB)
	@echo " ........... compile : debug   ........... "
	$(FC) $(DEB) $(OBJS_DEB) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT_DEB) $(SFLIBS_DEB) $(SFMODS_DEB) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)


pade: 	version $(OBJS)
	@echo " ........... compile: rdmft_pade ........... "
	$(FC) $(STD) $(OBJS) drivers/pade_matsubara_to_real.f90 -o $(DIREXE)/pade_matsubara_to_real $(LIBDMFT) $(SFLIBS) $(SFMODS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/pade_matsubara_to_real

# data: 	version $(OBJS)
# 	@echo " ........... compile: get_data ........... "
# 	$(FC) $(STD) $(OBJS) $(DIR)/get_data_$(EXE).f90 -o $(DIREXE)/get_data_$(EXE) $(LIBDMFT) $(SFLIBS) $(SFMODS) 
# 	@echo " ...................... done .............................. "
# 	@echo ""
# 	@echo "created" $(DIREXE)/get_data_$(EXE)


RDMFT_VARS_GLOBAL.o: RDMFT_VARS_GLOBAL.f90
	$(FC) $(STD) -c RDMFT_VARS_GLOBAL.f90 $(SFMODS) 
RDMFT_VARS_GLOBAL_DEB.o: RDMFT_VARS_GLOBAL.f90
	$(FC) $(DEB) -c RDMFT_VARS_GLOBAL.f90 $(SFMODS_DEB) -o  RDMFT_VARS_GLOBAL_DEB.o

clean: 
	@echo 'removing *.mod *.o *~'
	@rm -fv *.mod *.o *~ revision.inc


#########################################################################
include version.mk
#########################################################################
