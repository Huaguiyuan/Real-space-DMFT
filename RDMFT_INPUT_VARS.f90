!###############################################################
! PROGRAM  : RDMFT_VARS_GLOBAL
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
module RDMFT_INPUT_VARS
  USE SCIFOR_VERSION
  USE COMMON_VARS
  !Impurity solver interface
  USE SOLVER_INPUT_VARS
  !parallel library
  USE MPI
  implicit none

  !Revision software:
  !=========================================================
  include "revision.inc"

  !Lattice size:
  !=========================================================
  integer   :: Nside


  !Global variables
  !=========================================================
  real(8) :: Wdis               !Disorder strength
  integer :: idum               !disorder seed
  real(8) :: a0trap             !Trap bottom
  real(8) :: V0trap             !Trap curvature in x,y directions (assumed circular symmetry)
  real(8) :: chitrap            !Tentative value for the global trap compressibility dN/d(mu_{tot})
  integer :: N_wanted           !Required number of particles for canonical calculations [set 0 for fixmu]
  real(8) :: N_tol              !Tolerance over the number of particles


  !Other variables:
  !=========================================================
  logical           :: pbcflag
  logical           :: symmflag
  real(8)           :: rdmft_nread,rdmft_nerror,rdmft_ndelta
  integer           :: mix_type


  !Restart files
  !=========================================================
  character(len=64)          :: fileSig,fileSelf



  !Namelist:
  !=========================================================
  namelist/rdmft_vars/&
       Wdis,     &
       V0trap,   &
       a0trap,   &
       Nside,    &
       idum,     &
       rdmft_nread,    &
       rdmft_nerror,   &
       rdmft_ndelta,   &
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

    call version(revision)

    !local variables: default values
    Wdis            = 0.5d0
    Nside           = 5
    a0trap          = 0.d0
    v0trap          = 0.1d0
    rdmft_nread           = 0.d0
    rdmft_nerror          = 1.d-4
    rdmft_ndelta          = 0.1d0
    fileSig         = "LSigma.data"
    fileSelf        = "LSelf.data"
    symmflag        =.false.
    N_wanted        = Nside**2/2
    N_tol           = 0.1d0
    chitrap         = 0.1d0 
    pbcflag         = .true.
    mix_type        = 0
    idum            = 1234567

    !Read input file (if any)
    inquire(file=inputFILE,exist=control)
    if(control)then
       open(10,file=inputFILE)
       read(10,nml=rdmft_vars)
       close(10)
    else
       if(mpiID==0)then
          open(10,file="default."//inputFILE)
          write(10,nml=rdmft_vars)
          close(10)
          write(*,*) "can not open INPUT file, dumping a default version in +default."//inputFILE
       endif
       stop
    endif

    !Parse cmd.line arguments (if any)
    call parse_cmd_variable(wdis,"WDIS")
    call parse_cmd_variable(v0trap,"V0TRAP")
    call parse_cmd_variable(a0trap,"A0TRAP")
    call parse_cmd_variable(Nside,"NSIDE")
    call parse_cmd_variable(rdmft_nread,"RDMFT_NREAD")
    call parse_cmd_variable(rdmft_nerror,"RDMFT_NERROR")
    call parse_cmd_variable(rdmft_ndelta,"RDMFT_NDELTA")
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
       write(*,nml=rdmft_vars)
       open(10,file="used."//inputFILE)
       write(10,nml=rdmft_vars)
       close(10)
    endif

    
  end subroutine rdmft_read_input



end module RDMFT_INPUT_VARS
