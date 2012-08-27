
MODULE COMMON_BROYDN
  implicit none
  integer                :: isite    
  complex(8),allocatable :: fg(:,:),sigma(:,:)
  complex(8),allocatable :: fg0(:),gamma(:)
END MODULE COMMON_BROYDN

FUNCTION FUNCV(x)
  USE RDMFT_VARS_GLOBAL
  USE COMMON_BROYDN
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  xmu0=x(1)
  fg0 = one/(one/gamma +erandom(isite)+xmu-xmu0-U*(n-0.5d0))
  n0=sum(fermi(wr,beta)*aimag(fg0))/sum(aimag(fg0))
  funcv(1)=n-n0
  write(101+mpiID,"(3(f13.9))")n,n0,xmu0
END FUNCTION FUNCV
