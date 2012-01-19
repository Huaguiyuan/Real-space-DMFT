EXE	=hm_real_disorder
DIR	=drivers/hm_real
DIREXE	=$(HOME)/.bin
FC	=mpif90

#########################################################################
include $(HOME)/lib/lib.mk
include $(HOME)/lib/libdmft.mk
#########################################################################

OBJS=RDMFT_VARS_GLOBAL.o
OBJS_DEB=RDMFT_VARS_GLOBAL_DEB.o
all: 	$(OBJS)
	@echo " ........... compile: standard ........... "
	$(FC) $(STD) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(LIBS) $(MODS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)

opt: 	$(OBJS)
	@echo " ........... compile: optimized ........... "
	$(FC) $(OPT) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(LIBS) $(MODS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)


debug: 	$(OBJS_DEB)
	@echo " ........... compile : debug   ........... "
	$(FC) $(DEB) $(OBJS_DEB) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT_DEB) $(LIBS_DEB) $(MODS_DEB) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)


data: 	$(OBJS)
	@echo " ........... compile: get_data ........... "
	$(FC) $(STD) $(OBJS) $(DIR)/get_data_$(EXE).f90 -o $(DIREXE)/get_data_$(EXE) $(LIBDMFT) $(LIBS) $(MODS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/get_data_$(EXE)


RDMFT_VARS_GLOBAL.o: RDMFT_VARS_GLOBAL.f90
	$(FC) $(STD) -c RDMFT_VARS_GLOBAL.f90 $(MODS) 
RDMFT_VARS_GLOBAL_DEB.o: RDMFT_VARS_GLOBAL.f90
	$(FC) $(DEB) -c RDMFT_VARS_GLOBAL.f90 $(MODS_DEB) -o  RDMFT_VARS_GLOBAL_DEB.o

clean: 
	@echo 'removing *.mod *.o *~'
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
