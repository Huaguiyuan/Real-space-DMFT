#=========================================================================
include sfmake.inc
#=========================================================================
FC=ifort

OBJS = SOLVER_INPUT_VARS  IPT_VARS_GLOBAL.o IPT_MATS.o   IPT_SOPT.o IPT_SC_SOPT.o DMFT_IPT.o

#=================STANDARD COMPILATION====================================
all: FLAG=$(STD)
all: ARGS=$(SFLIBS)
all: libipt_rdmft.a


#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB)
debug: ARGS=$(SFLIBS_DEB)
debug:libipt_rdmft_deb.a


libipt_rdmft.a: $(OBJS)
	rm -fv $@
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@
	chmod +x $@
	rm -fv *.o
	@echo "LIBRARY $@ IS DONE........."
	@echo ' '


libipt_rdmft_deb.a:  $(OBJS)
	rm -fv $@
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@
	chmod +x $@
	rm -fv *.o
	@echo "LIBRARY $@ IS DONE........."
	@echo ' '


.f90.o:	
	$(FC) $(FLAG) -c $< $(SFINCLUDE) $(SFLIB) #$(MOPT)../ 
SOLVER_INPUT_VARS:
	$(FC) $(FLAG) -c ../$@.f90 $(SFINCLUDE) $(SFLIB) #$(MOPT)../ 
	$(FC) $(FLAG) -c ../$@.f90 $(SFINCLUDE) $(SFLIB) $(MOPT)../ 

clean: 
	@echo 'removing *.mod *.o *~'
	@rm -fv *.mod *.o *~ *.a

