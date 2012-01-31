  !     START MPI:
  !=====================================================================
  call start_mpi()

  !     READ INPUT FILES:
  !=====================================================================
  call read_input("inputIPT.in")
  call rdmft_read_input("inputRDMFT.in")

  !     BUILD MATSUBARA/TAU MESH   
  !==================================================================== 
  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)

  if (mpiID==0) then
     write(*,*)"======================================================"
     write(*,"(A,I9,A)")"Using ",L," frequencies"
  endif

  wmax  =  wmax+Wdis  ! real frequency mesh for spectral functions (unused)

  !     ALLOCATE WORKING ARRAYS:
  !=====================================================================
  Ns    =Nside**2       !  per il momento alloco tutto... 
  n     =0.5d0          !  se ho problemi di memoria allochero' solo quello che mi serve
                        !  il numero di siti indipendenti per simmetria cubica in 2d
                        !   e' dato da

  if (mpiID==0) then 
     write(*,*)"======================================================"
     write(*,"(A,I9,A)")"System size =",Ns 
  endif
       
  allocate(etrap(Ns))  
  allocate(H0(Ns,Ns))                 ! tb Hamiltonian in principle I can treat every quadratic Hamiltonian ...
  allocate(icol(Ns),irow(Ns))         ! questo non va cambiato    
  allocate(ij2site(-Nside/2:Nside/2,-Nside/2:Nside/2))      ! definisco dei nuovi boundaries per il mapping da i a isite,jsite                                                                 ! NB Nside e' sempre dispari nel programma
  allocate(nii(Ns))
  allocate(dii(Ns))




  !     BUILD THE LATTICE HAMILTONIAN:
  !=====================================================================
  call get_tb_hamiltonian

  if (symmflag) then
   Nindip =(Nside**2+4*Nside+3)/8
   if (mpiID==0) then 
      write(*,"(A,I4,A)")"Using Trap Symmetries: only ",Nindip,"  independent lattice sites" 
      write(*,*)"======================================================"
   endif
   allocate(indipsites(Nindip))
   call get_indip_list
  endif

  
  if (N_wanted==0) then 
     densfixed=.false.
     if (mpiID==0) then
        write(*,*)"======================================================" 
        write(*,*)"Working at fixed (global) chemical potential"
        write(*,*)""
     endif
  else
     densfixed=.true.
     deltan=1.0d0                    ! to be read from input ?
     if (mpiID==0) then 
        write(*,*)"======================================================"
        write(*,*)"Working at fixed total particle number",N_wanted
        write(*,*)"Required tolerance over the number of particles",N_tol
        write(*,*)"Starting value for the trap compressibility",chitrap
        write(*,*)"Initial step in mu",deltan
        write(*,*)""
    endif
 endif

  !     BUILD THE TRAP:
  !=====================================================================
  do is=1,Ns
     etrap(is)= 0.5d0*V0trap*trap_distance_square(is) 
! + a0trap now it is defined outside so we can change it in the loops
  enddo

! if(mpiID==0)call splot("fort.53",etrap)


  !     CREATE DATA_DIR FOR TRAP DATA
  name_dir="data_trap"
  if(mpiID==0)then
     call create_data_dir(trim(adjustl(trim(name_dir))))
     write(*,*)"======================================================"
  endif

