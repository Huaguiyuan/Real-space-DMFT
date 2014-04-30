module RDMFT_VARS_GLOBAL
  USE CONSTANTS
  USE MPI_VARS
  USE MPI
  implicit none


  !Revision software:
  !=========================================================
  include "revision.inc"


  !Lattice size:
  !=========================================================
  integer                             :: Nlat,Nindip

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable    :: wm,tau
  real(8),dimension(:),allocatable    :: wr,t


  !Large matrices for Lattice Hamiltonian/GF
  !=========================================================
  integer,dimension(:),allocatable    :: icol,irow
  integer,dimension(:,:),allocatable  :: ij2site
  integer,dimension(:), allocatable   :: indipsites
  real(8),dimension(:,:),allocatable  :: H0,Id


  !Local density and order parameter profiles:
  !=========================================================
  real(8),dimension(:),allocatable    :: nii,dii,pii,gap_ii
  complex(8),dimension(:),allocatable :: cdii
  logical                             :: densfixed


end module RDMFT_VARS_GLOBAL
