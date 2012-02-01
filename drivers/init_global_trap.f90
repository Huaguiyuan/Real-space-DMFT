  !START MPI:
  !=====================================================================
  call start_mpi()

  !READ INPUT FILES:
  !=====================================================================
  call read_input("inputIPT.in")
  call rdmft_read_input("inputRDMFT.in")

  !CHECK ODD NSIDE FOR TRAP:
  !=====================================================================
  if(mod(Nside,2)==0)then
     Nside=Nside-1
     if(mpiID==0)then 
        write(*,"(A)")bg_red("Nside has to be odd!")
        write(*,"(A,I,A)")bg_red("Using Nside="),Nside,bg_red(" instead")
     endif
  endif


  !ALLOCATE WORKING ARRAYS:
  !=====================================================================
  Ns    =Nside**2
  wmax  =wmax+Wdis
  allocate(etrap(Ns))
  allocate(H0(Ns,Ns))
  allocate(icol(Ns),irow(Ns))
  !definisco dei nuovi boundaries per il mapping da i a isite,jsite
  !NB Nside e' sempre dispari nel programma
  allocate(ij2site(-Nside/2:Nside/2,-Nside/2:Nside/2))  
  allocate(nii(Ns))
  allocate(dii(Ns))


  if(mpiID==0)write(*,"(A,I9)")"System size =",Ns


  !CREATE DATA_DIR FOR AND STORE THIS IDUM 
  name_dir="data_trap"
  if(mpiID==0)then
     call create_data_dir(trim(adjustl(trim(name_dir))))
  endif


  !BUILD THE LATTICE HAMILTONIAN:
  !=====================================================================
  call get_tb_hamiltonian(centered=.true.)


  !REDUCE THE PROBLEM BY TAKING INTO ACCOUTN TRAP SYMMETRIES:
  !=====================================================================
  if (symmflag) then
     Nindip = (Nside**2+4*Nside+3)/8
     if (mpiID==0)write(*,"(A,I4,A)")"Using Trap Symmetries: ",Nindip," independent lattice sites." 
     allocate(indipsites(Nindip))
     call get_indip_list
  endif


  !WORKING WITH FIXED TOTAL OCCUPATION:
  !=====================================================================
  if (N_wanted==0) then 
     densfixed=.false.
     if (mpiID==0) then
        write(*,*)"Working at fixed (global) chemical potential"
     endif
  else
     densfixed=.true.
     if (mpiID==0) then 
        write(*,"(A,I6)")"Working at fixed total particle number         =",N_wanted
        write(*,"(A,F12.9)")"Required tolerance over the number of particles=",N_tol
        write(*,"(A,F12.9)")"Starting value for the trap compressibility    =",chitrap
        write(*,"(A,F12.9)")"Initial step in mu                             =",ndelta
     endif
  endif


  !BUILD THE TRAP:
  !=====================================================================
  do is=1,Ns
     etrap(is)= 0.5d0*V0trap*trap_distance_square(is)
  enddo
