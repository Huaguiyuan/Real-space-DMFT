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
module RDMFT_VARS_GLOBAL
  !Scientific library
  USE COMMON_VARS
  
  ! !Impurity solver interface
  ! USE SOLVER_INPUT_VARS
  !parallel library
  USE MPI
  implicit none


  !Lattice size:
  !=========================================================
  integer   :: Nlat,Nindip

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


  !Local density and order parameter profiles:
  !=========================================================
  real(8),dimension(:),allocatable    :: nii,dii,gap_ii
  complex(8),dimension(:),allocatable :: cdii
  logical                             :: densfixed


end module RDMFT_VARS_GLOBAL
