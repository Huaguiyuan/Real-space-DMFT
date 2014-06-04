FC=$(SFMPI)/mpif90
#=========================================================================
include sfmake.inc
#=========================================================================

DIREXE=$(HOME)/.bin
DIR=drivers

#EXE=ahm_real_trap
#EXE=ahm_real_disorder
#EXE=ed_ahm_squareRDMFT
#EXE=ed_ahm_slabRDMFT
EXE=ed_ahm_disorder
#EXE=ed_hm_disorder

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

OBJS = RDMFT_INPUT_VARS.o RDMFT_VARS_GLOBAL.o RDMFT_AUX_FUNX.o RDMFT_WRAP_IPT.o RDMFT_WRAP_ED.o RDMFT.o


#=================STANDARD COMPILATION====================================
normal: FLAG=$(STD)
normal: ARGS= -I./IPT_SOLVER -I./ED_SOLVER -L./IPT_SOLVER -L./ED_SOLVER -lipt_rdmft -led_rdmft $(SFLIBS)
normal:compile

#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB)
debug: ARGS= -I./IPT_SOLVER -I./ED_SOLVER -L./IPT_SOLVER -L./ED_SOLVER -lipt_rdmft_deb -led_rdmft_deb $(SFLIBS_DEB)
debug: compile

all: lib normal
all_debug: lib_debug debug

lib: ipt_solver ed_solver
lib_debug: ipt_solver_debug ed_solver_debug

#=================STANDARD COMPILATION====================================
compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

ipt_solver:
	@make -C IPT_SOLVER/

ed_solver:
	@make -C ED_SOLVER/

ipt_solver_debug:
	@make debug -C IPT_SOLVER/

ed_solver_debug:
	@make debug -C ED_SOLVER/

.f90.o:	
	$(FC) $(FLAG) -c $< $(SFINCLUDE) -I./IPT_SOLVER -I./ED_SOLVER


clean: 
	@echo 'removing *.mod *.o *~'
	@rm -fv *.mod *.o *~ revision.inc

all_clean: clean
	@make -C IPT_SOLVER/ clean
	@make -C ED_SOLVER/ clean



version:
	@echo $(VER)
