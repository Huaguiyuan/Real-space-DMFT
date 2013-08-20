!###############################################################
! PROGRAM  : RDMFT_VARS_GLOBAL
! AUTHORS  : Adriano Amaricci & Antonio Privitera
! NAME
!   xxx_disorder/trap  
! DESCRIPTION
!   This a layer interface to different codes solving the Real-space Dynamical Mean Field Theory.
!   Each code must be implemented separetely, for convenience some driver routines are in drivers/ dir.
!   The structure is very easy: a tight-binding hamiltonian is generated at first and the DMFT problem  
!   is solved in the local wannier basis, solving the impurity problem at each lattice site. Simmetries 
!   can be used to reduce the size of the problem. Lattice sites are related via self-consistency. 
!   The code is MPI parallel.     
! OPTIONS
!  wdis=[0.5]    -- degree of local disorder.
!  Nside=[10]    -- linear size of the cluster to be solved.
!  a0trap=[0]    -- bottom of the trap. here kept separated from mu.
!  v0trap=[0.1]  -- fix the parabolic shape of the trap.
!  nread=[0.0]   -- density value for chemical potential search.
!  ndelta=[0.1]  -- starting value for chemical potential shift.
!  nerror=[1.d-4]-- max error in adjusting chemical potential. 
!  symmflag=[T]  -- Enforce trap cubic symmetry in the xy-plane.
!  n_wanted=[0]  -- Required number of particles in the trap. Fix to 0 (default) for mufixed
!  n_tol=[0.1]   -- Tolerance over the total density
!  chitrap=[0.1] -- Tentative value of the global trap compressibility
!  pbcflag=[T]   -- periodic boundary conditions.
!  idum=[1234567]-- initial seed for the random variable sample.	
!###############################################################
module RDMFT_VARS_GLOBAL
  !Scientific library
  USE COMMON_VARS
  USE SCIFOR_VERSION
  USE TIMER, ONLY:start_timer,stop_timer,eta
  USE IOTOOLS
  USE MATRIX
  USE RANDOM,    ONLY:nrand,init_random_number
  USE STATISTICS
  USE INTEGRATE, ONLY:kronig
  USE FUNCTIONS, ONLY:fermi
  USE TOOLS,     ONLY:check_convergence
  !Impurity solver interface
  USE SOLVER_VARS_GLOBAL
  !parallel library
  USE MPI
  implicit none

  !Revision software:
  !=========================================================
  include "revision.inc"

  !Lattice size:
  !=========================================================
  integer   :: Nside,Ns,Nindip,iloop

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau
  real(8),dimension(:),allocatable :: wr,t


  !Large matrices for Lattice Hamiltonian/GF
  !=========================================================
  integer,dimension(:),allocatable   :: icol,irow
  integer,dimension(:,:),allocatable :: ij2site
  integer,dimension(:), allocatable  :: indipsites !***to be renamed***
  real(8),dimension(:,:),allocatable :: H0,Id


  !Local density and order parameter profiles:
  !=========================================================
  real(8),dimension(:),allocatable    :: nii,dii,gap_ii
  complex(8),dimension(:),allocatable :: cdii


  !Global variables
  !=========================================================
  real(8) :: Wdis               !Disorder strength
  integer :: idum               !disorder seed
  real(8) :: a0trap             !Trap bottom
  real(8) :: V0trap             !Trap curvature in x,y directions (assumed circular symmetry)
  !***To be renamed***
  integer :: N_wanted           !Required number of particles for canonical calculations [set 0 for fixmu]
  real(8) :: N_tol              !Tolerance over the number of particles
  real(8) :: chitrap            !Tentative value for the global trap compressibility dN/d(mu_{tot})


  !Other variables:
  !=========================================================
  logical           :: pbcflag
  logical           :: symmflag
  logical           :: densfixed
  real(8)           :: nread,nerror,ndelta
  integer           :: mix_type
  ! !Random energies
  ! !=========================================================
  ! real(8),dimension(:),allocatable   :: erandom,etrap


  !Restart files
  !=========================================================
  character(len=64)          :: fileSig,fileSelf



  !Namelist:
  !=========================================================
  namelist/disorder/&
       Wdis,     &
       V0trap,   &
       a0trap,   &
       Nside,    &
       idum,     &
       nread,    &
       nerror,   &
       ndelta,   &
       fileSig,  &
       fileSelf, &
       symmflag, &
       N_wanted, &
       N_tol,    &
       chitrap,  &
       mix_type, &   
       pbcflag



contains

  !+----------------------------------------------------------------+
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine rdmft_read_input(inputFILE)
    character(len=*) :: inputFILE
    integer          :: i
    logical          :: control
    !local variables: default values
    Wdis            = 0.5d0
    Nside           = 5
    a0trap          = 0.d0
    v0trap          = 0.1d0
    nread           = 0.d0
    nerror          = 1.d-4
    ndelta          = 0.1d0
    fileSig         = "LSigma.data"
    fileSelf        = "LSelf.data"
    symmflag        =.false.
    N_wanted        = Nside**2/2
    N_tol           = 0.1d0
    chitrap         = 0.1d0 
    pbcflag         = .true.
    mix_type        = 0
    idum            = 1234567


    !SET SIZE THRESHOLD FOR FILE ZIPPING:
    store_size=1024

    !Read input file (if any)
    inquire(file=adjustl(trim(inputFILE)),exist=control)
    if(control)then
       open(10,file=adjustl(trim(inputFILE)))
       read(10,nml=disorder)
       close(10)
    else
       open(10,file="default."//adjustl(trim(inputFILE)))
       write(10,nml=disorder)
       close(10)
       call abort("can not open INPUT file, dumping a default version in +default."//adjustl(trim(inputFILE)))
    endif

    !Parse cmd.line arguments (if any)
    call parse_cmd_variable(wdis,"WDIS")
    call parse_cmd_variable(v0trap,"V0TRAP")
    call parse_cmd_variable(a0trap,"A0TRAP")
    call parse_cmd_variable(Nside,"NSIDE")
    call parse_cmd_variable(nread,"NREAD")
    call parse_cmd_variable(nerror,"NERROR")
    call parse_cmd_variable(ndelta,"NDELTA")
    call parse_cmd_variable(n_wanted,"NWANTED")
    call parse_cmd_variable(n_tol,"NTOL")
    call parse_cmd_variable(fileSig,"FILESIG")
    call parse_cmd_variable(fileSelf,"FILESELF")
    call parse_cmd_variable(chitrap,"CHITRAP")
    call parse_cmd_variable(symmflag,"SYMMFLAG")
    call parse_cmd_variable(pbcflag,"PBCFLAG")
    call parse_cmd_variable(mix_type,"MIX_TYPE")
    call parse_cmd_variable(idum,"IDUM")

    !Print on the screen used vars
    if(mpiID==0)then
       write(*,nml=disorder)
       open(10,file="used."//adjustl(trim(inputFILE)))
       write(10,nml=disorder)
       close(10)
    endif

    !------SET NUMBER OF LATTICE SITES------!
    Ns    =Nside**2
    !---------------------------------------!

    allocate(nii(Ns))
    allocate(dii(Ns))

    !STORE THIS IDUM 
    !=====================================================================
    if(mpiID==0)then
       open(10,file="list_idum",access='append')
       write(10,*)idum
       close(10)
    endif

    call version(revision)
  end subroutine rdmft_read_input



end module RDMFT_VARS_GLOBAL
