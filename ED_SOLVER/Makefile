include sfmake.inc


#COMPILATION:
OBJS= RDMFT_INPUT_VARS MATRIX_SPARSE.o ED_EIGENSPACE.o ED_BATH_TYPE.o ED_INPUT_VARS.o ED_VARS_GLOBAL.o ARPACK_LANCZOS.o PLAIN_LANCZOS.o ED_AUX_FUNX.o ED_BATH.o ED_HAMILTONIAN.o ED_GREENS_FUNCTIONS.o ED_OBSERVABLES.o ED_CHI2FIT.o ED_DIAG.o DMFT_ED.o

#=================STANDARD COMPILATION====================================
all: FLAG=$(STD)
all: ARGS=$(SFLIBS)
all: libed_rdmft.a

#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB)
debug: ARGS=$(SFLIBS_DEB)
debug:libipt_rdmft_deb.a



libed_rdmft.a: $(OBJS)
	rm -fv $@
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@
	chmod +x $@
	rm -fv *.o
	@echo "LIBRARY $@ IS DONE........."
	@echo ' '


libed_rdmft_deb.a:  $(OBJS)
	rm -fv $@
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@
	chmod +x $@
	rm -fv *.o
	@echo "LIBRARY $@ IS DONE........."
	@echo ' '


.f90.o:	
	$(FC) $(FLAG) -c $< $(SFINCLUDE)
RDMFT_INPUT_VARS:
	$(FC) $(FLAG) -c ../$@.f90 $(SFINCLUDE) 
#$(FC) $(FLAG) -c ../$@.f90 $(SFINCLUDE) $(MOPT)../ 
clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc *.a

sync_ed:
	@rsync -avP --del --exclude=ED_INPUT_VARS.f90 $(HOME)/projects_src/ED/*.f90 .
