!###############################################################
!     PROGRAM  : VARS_GLOBAL
!     TYPE     : Module
!     PURPOSE  : Contains global variables
!     AUTHORS  : Adriano Amaricci
!###############################################################
module RDMFT_VARS_GLOBAL
  !LOAD MODULE FROM LIBRARY:
  USE COMMON_VARS
  USE CHRONOBAR, ONLY:start_timer,stop_timer,eta
  USE IOTOOLS
  USE SLPLOT
  USE SLREAD
  USE MATRIX,    ONLY:mat_inversion_sym,mat_inversion_gj,mat_inversion
  USE RANDOM,    ONLY:nrand,drand,init_random_number
  USE INTEGRATE, ONLY:kronig
  USE TOOLS,     ONLY:fermi,check_convergence
  USE DMFT_IPT
  USE MPI
  USE OMP_LIB
  implicit none

  !Lattice size:
  !=========================================================
  integer,protected :: Nside
  integer           :: Ns,Nindip

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau
  real(8),dimension(:),allocatable :: wr,t

  !Large matrices for Lattice Hamiltonian/GF
  !=========================================================
  integer,dimension(:),allocatable   :: icol,irow
  integer,dimension(:,:),allocatable :: ij2site
  integer,dimension(:), allocatable  :: indipsites
  real(8),dimension(:,:),allocatable :: H0,Id
  real(8),dimension(:),allocatable   :: nii,dii !density profile and delta profile 

  !Global variables
  !=========================================================
  real(8) :: Wdis               !Disorder strength
  integer :: idum               !Disorder seed
  real(8) :: a0trap             !Trap bottom
  real(8) :: V0trap             !Trap curvature in x,y directions (assumed circular symmetry)
  integer :: N_wanted           !Required number of particles for canonical calculations [set 0 for fixmu]
  real(8) :: N_tol              !Tolerance over the number of particles
  real(8) :: chitrap            !Tentative value for the global trap compressibility dN/d(mu_{tot})
  real(8) :: gammatrap          !Trap_asymmetry in the third dimension (not implemented yet)
  integer :: dim                !Spatial dimension (1,2,3) not implemented yet     

  !Other variables:
  !=========================================================
  character(len=20)                  :: name_dir
  real(8)                            :: n,n0,xmu0
  logical                            :: pbcflag
  logical                            :: symmflag
  logical                            :: densfixed

  !Random/trap energies
  !=========================================================
  real(8),dimension(:),allocatable   :: erandom,etrap

  !Namelists:
  !=========================================================
  namelist/disorder/&
       Wdis,     &
       V0trap,   &
       a0trap,   &
       Nside,    &
       idum,     &
       pbcflag,  &
       symmflag, &
       N_wanted, &
       N_tol,    &
       chitrap,  &   
       gammatrap,&
       dim,      &
       omp_num_threads

contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine rdmft_read_input(inputFILE)
    character(len=*) :: inputFILE
    integer :: i
    !local variables: default values
    Wdis            =0.5d0
    Nside           =9
    a0trap          = 0.d0
    v0trap          = 0.1d0
    omp_num_threads = 1
    pbcflag         = .true.
    symmflag        =.false.
    N_wanted        = Nside**2/2
    N_tol            = 0.1d0
    chitrap         = 0.1d0 
    gammatrap       =1.d0
    dim             = 2
    idum            =1234567
  
    !SET SIZE THRESHOLD FOR FILE ZIPPING:
    store_size=1024

    open(10,file=adjustl(trim(inputFILE)))
    read(10,nml=disorder)
    close(10)
    allocate(help_buffer(26))
    help_buffer=([&
         'NAME',&
         '  xxx_disorder/trap',&
         '',&
         'DESCRIPTION',&
         '  This a layer interface to different codes solving the Real-space Dynamical Mean Field Theory.',&
         '  Each code must be implemented separetely, for convenience some driver routines are in drivers/ dir.',&
         '  The structure is very easy: an thight-binding hamiltonian is generated at fist and the DMFT problem  ',&
         '  is solved in the local wannier basis, solving the impurity problem at each lattice site. Simmetries ',&
         '  can be used to reduce the size of the problem. Lattice sites are related via self-consistency. ',&
         '  The code is MPI parallel. ',&
         '  ',&
         'OPTIONS',&
         ' wdis=[0.5]     -- degree of local disorder.',&
         ' Nside=[9]      -- linear size of the cluster to be solved.',&
         ' a0trap=[0.0]   -- bottom of the trap. here kept separated from mu.',&
         ' v0trap=[0.1]   -- fix the parabolic shape of the trap.',&
         ' n_wanted=[0]   -- Required number of particles in the trap. Fix to 0 (default) for mufixed',&
         ' n_tol=[0.1]    -- Tolerance over the total density',&
         ' chitrap=[0.1]  -- Tentative value of the global trap compressibility',&
         ' gammatrap=[1.0]-- Trap asymmetry in the z-direction (not yet implemented)',&
         ' dim=[2]        -- Dimension of space (not yet implemented)',&
         ' pbcflag=[T]    -- periodic boundary conditions.',&
         ' symmflag=[T]   -- Enforce trap cubic symmetry in the xy-plane.',&
         ' idum=[1234567]-- initial seed for the random variable sample.',&
         ' omp_num_threads=[1] -- fix the number of threads in OMP environment.',&
         '  '])
    call parse_cmd_help(help_buffer)
    call parse_cmd_variable(wdis,"WDIS")
    call parse_cmd_variable(v0trap,"V0TRAP")
    call parse_cmd_variable(a0trap,"A0TRAP")
    call parse_cmd_variable(Nside,"NSIDE")
    call parse_cmd_variable(n_wanted,"N_WANTED")
    call parse_cmd_variable(n_tol,"N_TOL")
    call parse_cmd_variable(chitrap,"CHITRAP")
    call parse_cmd_variable(gammatrap,"GAMMATRAP")
    call parse_cmd_variable(dim,"DIM")
    call parse_cmd_variable(pbcflag,"PBCFLAG")
    call parse_cmd_variable(symmflag,"SYMMFLAG")
    call parse_cmd_variable(omp_num_threads,"OMP_NUM_THREADS")
    call parse_cmd_variable(idum,"IDUM")
    
!    n_command_arg=command_argument_count()
!    open(10,file=adjustl(trim(inputFILE)))
!    read(10,nml=disorder)
!    close(10)

!    !Process command line variable change:
!    if(n_command_arg/=0)then
!       do i=1,n_command_arg
!          call get_command_argument(i,arg_buffer)
!          call cmd_help(arg_buffer,&
!               "exe nside=value[10] pbcflag=value[T] disorder[wdis=value[0.5] idum=value[1234567]] trap[v0trap=value[1] a0trap=value[1]]")
!          call cmd_var(arg_buffer)
!          select case(nml_name)
!          case("WDIS");read(nml_value,*)Wdis

!          case("V0TRAP");read(nml_value,*)v0trap
!          case("A0TRAP");read(nml_value,*)a0trap
!          case("NSIDE");read(nml_value,*)Nside
!          case("PBCFLAG");read(nml_value,*)pbcflag
!          case("SYMMFLAG");read(nml_value,*)symmflag
!          case("N_WANTED");read(nml_value,*)N_wanted
!          case("GAMMA");read(nml_value,*)gamma
!          case("DIM");read(nml_value,*)dim   
!          case("OMP_NUM_THREADS");read(nml_value,*)omp_num_threads
!          case("IDUM");read(nml_value,*)idum
!          case default
!             print*,"No corresponging variable in NML"
!          end select
!       enddo
!    endif


!!   it is better to have the center of the trap on a given site
!!   to avoid spurious degeneracies
  


    !SET OMP THREADS NUMBER
    call omp_set_num_threads(OMP_NUM_THREADS)

    !Print on the screen used vars
    if(mpiID==0)then
       write(*,nml=disorder)
       open(10,file="used."//adjustl(trim(inputFILE)))
       write(10,nml=disorder)
       close(10)
    endif

    if (mod(Nside,2)==0)then
       Nside=Nside-1 
       if (mpiID==0)then 
         write(*,*)"Nside has to be odd"
         write(*,*)"using instead Nside=",Nside
       endif
   endif

  end subroutine rdmft_read_input
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : Build tight-binding Hamiltonian
  !+----------------------------------------------------------------+
  subroutine get_tb_hamiltonian()
    integer :: i,jj,j,k,row,col,link(4),istate,jstate
    H0=0.d0

!    write(534,*)"Tabella dei siti"
!    write(534,*)"i,irow(i),icol(i)"
!    this should be changed to accomodate first the indipendent sites..TODO

    do row=0,Nside-1
       do col=0,Nside-1
          i=col+row*Nside+1                            ! questo lo lascio uguale

!          irow(i)=row+1
!          icol(i)=col+1                               ! vecchio
!          ij2site(row+1,col+1)=i
 
          irow(i)=-Nside/2+row                     ! cambio la tabella i -> isite,jsite
          icol(i)=-Nside/2+col                     ! per farla simmetrica.. aiutera' 
          ij2site(row-Nside/2,col-Nside/2)=i       ! a implementare le simmetrie

!          write(534,*)i,irow(i),icol(i)

          if(pbcflag)then
             !HOPPING w/ PERIODIC BOUNDARY CONDITIONS
             link(1)= row*Nside+1              + mod(col+1,Nside)  ;
             link(3)= row*Nside+1              + (col-1)           ; if((col-1)<0)link(3)=(Nside+(col-1))+row*Nside+1
             link(2)= mod(row+1,Nside)*Nside+1 + col               ; 
             link(4)= (row-1)*Nside+1          + col               ; if((row-1)<0)link(4)=col+(Nside+(row-1))*Nside+1
          else
             !without PBC = closed boundary conditions
             link(1)= row*Nside+1              + col+1   ; if((col+1)==Nside)link(1)=0
             link(3)= row*Nside+1              +(col-1)  ; if((col-1)<0)     link(3)=0
             link(2)= (row+1)*Nside+1 + col              ; if((row+1)==Nside)link(2)=0
             link(4)= (row-1)*Nside+1          + col     ; if((row-1)<0)     link(4)=0
          endif
          do jj=1,4
             if(link(jj)>0)H0(i,link(jj))=ts
          enddo
       enddo
    enddo

!    do istate=1,Ns
!       do jstate=1,Ns
!       write(44,*)istate,jstate,H0(istate,jstate)
!       enddo
!    enddo

  end subroutine get_tb_hamiltonian

  !******************************************************************
  !******************************************************************
  !******************************************************************

    !+----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : build state vector of indipendent sites 
  !+----------------------------------------------------------------+
  subroutine get_indip_list()
    integer :: i,row,col,istate,jstate

    i=0
    do col=0,Nside/2
       do row=0,col
          i= i+1
          indipsites(i)=ij2site(row,col)
!          write(424,*)i,row,col,indipsites(i)
       enddo
    enddo
  end subroutine get_indip_list

end module RDMFT_VARS_GLOBAL
