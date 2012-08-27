!########################################################
!     Program  : DHMMPT_MATSUBARA
!     TYPE     : Main program
!     PURPOSE  :solve the disordered (D) Hubbard model (HM) 
!using Modified Perturbation THeory (MPT), w/ DMFT
!     COMMENTS :
!DISORDER REALIZATION DEPENDS ON THE PARAMETER int idum:
!SO THAT DIFFERENT REALIZATIONS (STATISTICS) ARE PERFORMED 
!CALLING THIS PROGRAM MANY TIMES WHILE PROVIDING A *DIFFERENT* SEED.
!THE RESULT OF EACH CALCULATION IS STORED IN A DIFFERENT DIR
!     AUTHORS  : A.Amaricci
!########################################################
module COMMON_BROYDN
  implicit none
  integer                :: isite
  complex(8),allocatable :: fg(:,:),sigma(:,:)
  complex(8),allocatable :: fg0(:),gamma(:)
  real(8),allocatable    :: fgt(:,:),fg0t(:)
end module COMMON_BROYDN

function funcv(x)
  USE RDMFT_VARS_GLOBAL
  USE COMMON_BROYDN
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  xmu0=x(1)
  fg0 = one/(one/gamma +erandom(isite)+xmu-xmu0-U*(n-0.5d0))
  call fftgf_iw2tau(fg0,fg0t,beta)
  n0=-real(fg0t(L))
  funcv(1)=n-n0
  write(101+mpiID,"(3(f13.9))")n,n0,xmu0
end function funcv

program hmmpt_matsubara
  !LOCAL:
  USE RDMFT_VARS_GLOBAL
  USE COMMON_BROYDN
  implicit none
  integer    :: is
  real(8)    :: x(1),r
  logical    :: check,converged  


  !GLOBAL INITIALIZATION:
  !=====================================================================
  include "init_global.f90"


  !ALLOCATE WORKING ARRAYS:
  !=====================================================================
  allocate(fg(Ns,L),sigma(Ns,L))
  allocate(fg0(L),gamma(L))
  allocate(fgt(Ns,0:L),fg0t(0:L))


  !START DMFT LOOP SEQUENCE:SOLVE FOR \SIGMA_II(W), G_II(W)=FG
  !=====================================================================
  call setup_initial_sigma()
  iloop=0 ; converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     if(mpiID==0)write(*,"(A,I5,A,L6)")"DMFT-loop",iloop," convergence:",converged
     call get_gloc_mpi()       !SOLVE G_II (GLOCAL) \FORALL FREQUENCY W
     call solve_impurity_mpi() !SOLVE IMPURITY MODEL, \FORALL LATTICE SITES
     converged=check_convergence(sigma,eps_error,Nsuccess,nloop,id=0)
     call print_out(converged)
     call msg("============================================",lines=2)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
  enddo
  if(mpiID==0)call system("mv -vf *.err "//trim(adjustl(trim(name_dir)))//"/")
  call close_mpi()


contains

  include "solve_matsubara_routines.f90"

end program hmmpt_matsubara






