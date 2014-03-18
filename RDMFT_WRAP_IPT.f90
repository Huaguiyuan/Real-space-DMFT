!###############################################################
! PROGRAM  : RDMFT_WRAP_IPT
! PURPOSE  : Contains the main function performin RDMFT calculation
! for the IPT solver
!###############################################################
module RDMFT_WRAP_IPT
  USE RDMFT_INPUT_VARS
  USE RDMFT_VARS_GLOBAL
  USE TIMER
  USE FFTGF
  USE IOTOOLS,   only:reg,sread
  USE RANDOM,    only:nrand,init_random_number
  USE STATISTICS
  USE FUNCTIONS, only:fermi
  !Impurity solver interface
  USE DMFT_IPT
  implicit none
  private


  public :: ipt_solve_sc_impurity_mats
  public :: ipt_solve_sc_impurity_real


contains


  !+----------------------------------------------------------------+
  !PURPOSE  : parallel solution of impurity problem SC-version
  !+----------------------------------------------------------------+
  subroutine ipt_solve_sc_impurity_mats(fg,sigma)
    integer    :: is,i
    complex(8) :: fg(2,Nlat,Lmats),sigma(2,Nlat,Lmats)
    real(8)    :: nii_tmp(Nlat),dii_tmp(Nlat)
    complex(8) :: sigma_tmp(2,Nlat,Lmats)
    if(mpiID==0)write(LOGfile,*)"Solve impurity:"
    call start_timer
    sigma_tmp=zero
    nii_tmp  =0.d0
    dii_tmp  =0.d0
    do is=1+mpiID,Nlat,mpiSIZE
       call ipt_solve_per_site_mats(is,fg,sigma,sigma_tmp(:,is,:),nii_tmp(is),dii_tmp(is))
       call eta(is,Nlat,file="Impurity.eta")
    enddo
    call stop_timer
    call MPI_ALLREDUCE(sigma_tmp,sigma,2*Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ipt_solve_sc_impurity_mats

  subroutine ipt_solve_sc_impurity_real(fg,sigma)
    integer    :: is,i
    complex(8) :: fg(2,Nlat,Lreal),sigma(2,Nlat,Lreal)
    real(8)    :: nii_tmp(Nlat),dii_tmp(Nlat)
    complex(8) :: sigma_tmp(2,Nlat,Lreal)
    if(mpiID==0)write(LOGfile,*)"Solve impurity:"
    call start_timer
    sigma_tmp=zero
    nii_tmp  =0.d0
    dii_tmp  =0.d0
    do is=1+mpiID,Nlat,mpiSIZE
       call ipt_solve_per_site_real(is,fg,sigma,sigma_tmp(:,is,:),nii_tmp(is),dii_tmp(is))
       call eta(is,Nlat,file="Impurity.eta")
    enddo
    call stop_timer
    call MPI_ALLREDUCE(sigma_tmp,sigma,2*Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ipt_solve_sc_impurity_real







  !+----------------------------------------------------------------+
  !PURPOSE  : IPT solution of the is^th-impurity problem.
  !+----------------------------------------------------------------+
  subroutine ipt_solve_per_site_mats(is,fg,sigma,sig_tmp,ntmp,dtmp)
    integer,intent(in)                            :: is
    complex(8),dimension(2,Nlat,Lmats),intent(in) :: fg,sigma !Local GF,Self-energy
    complex(8),dimension(2,Lmats),intent(out)     :: sig_tmp
    real(8),intent(out)                           :: ntmp,dtmp
    complex(8),dimension(Lmats)                   :: det
    complex(8),dimension(:,:,:),allocatable,save  :: Wold
    complex(8),dimension(2,Lmats)                 :: calG,fg0
    real(8),dimension(2,0:Lmats)                  :: fgt,fg0t
    real(8)                                       :: n,n0,delta,delta0
    if(.not.allocated(Wold))allocate(Wold(2,Nlat,Lmats))
    if(mix_type==1)Wold(:,is,:) = sigma(:,is,:)
    !
    call fftgf_iw2tau(fg(1,is,:),fgt(1,0:Lmats),beta)
    call fftgf_iw2tau(fg(2,is,:),fgt(2,0:Lmats),beta,notail=.true.)
    n = -fgt(1,Lmats) ; delta = -u*fgt(2,Lmats)
    !
    ntmp=2.d0*n; dtmp=delta
    !
    fg0=zero ; calG=zero
    det       = abs(fg(1,is,:))**2    + fg(2,is,:)**2
    fg0(1,:)  = conjg(fg(1,is,:))/det + sigma(1,is,:) + U*(n-0.5d0)
    fg0(2,:)  = fg(2,is,:)/det        + sigma(2,is,:) + delta
    det       = abs(fg0(1,:))**2      + fg0(2,:)**2
    calG(1,:) = conjg(fg0(1,:))/det
    calG(2,:) = fg0(2,:)/det
    !
    ! if(mix_type==0)then
    !    if(iloop>1)calG(:,:) =  weight*calG(:,:) + (1.d0-weight)*Wold(:,is,:)
    !    Wold(:,is,:)  =  calG(:,:)
    ! endif
    !
    call fftgf_iw2tau(calG(1,:),fg0t(1,:),beta)
    call fftgf_iw2tau(calG(2,:),fg0t(2,:),beta,notail=.true.)
    n0=-fg0t(1,Lmats) ; delta0= -u*fg0t(2,Lmats)
    write(750,"(2I4,4(f16.12))",advance="yes")mpiID,is,n,n0,delta,delta0
    !
    sig_tmp(:,:) =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)

    sig_tmp(2,:) = dreal(sig_tmp(2,:))+zero    !!##ACTHUNG!!!
  end subroutine ipt_solve_per_site_mats

  subroutine ipt_solve_per_site_real(is,fg,sigma,sig_tmp,ntmp,dtmp)
    integer,intent(in)                            :: is
    complex(8),dimension(2,Nlat,Lreal),intent(in) :: fg,sigma !Local GF,Self-energy
    complex(8),dimension(2,Lreal),intent(out)     :: sig_tmp
    real(8),intent(out)                           :: ntmp,dtmp
    complex(8)                                    :: det
    complex(8),dimension(:,:,:),allocatable,save  :: Wold
    complex(8),dimension(2,1:Lreal)               :: calG,fg0
    real(8)                                       :: n,n0,delta,delta0
    integer                                       :: i
    if(.not.allocated(Wold))allocate(Wold(2,Nlat,1:Lreal))
    if(mix_type==1)Wold(:,is,:) = sigma(:,is,:)
    !
    n    = -sum(dimag(fg(1,is,:))*fermi(wr,beta))*fmesh/pi ! densita' per spin
    delta= -u*sum(dimag(fg(2,is,:))*fermi(wr,beta))*fmesh/pi
    !
    ntmp=2.d0*n; dtmp=delta
    !
    fg0=zero ; calG=zero 
    do i=1,Lreal
       det     = fg(1,is,i)*conjg(fg(1,is,Lreal+1-i)) + conjg(fg(2,is,Lreal+1-i))*fg(2,is,i)
       fg0(1,i)= conjg(fg(1,is,Lreal+1-i))/det  + sigma(1,is,i)         + u*(n-0.5d0)
       fg0(2,i)= conjg(fg(2,is,Lreal+1-i))/det  + sigma(2,is,Lreal+1-i) + delta
    end do
    do i=1,Lreal
       det     =  fg0(1,i)*conjg(fg0(1,Lreal+1-i)) + conjg(fg0(2,Lreal+1-i))*fg0(2,i)
       calG(1,i)=  conjg(fg0(1,Lreal+1-i))/det
       calG(2,i)=  conjg(fg0(2,Lreal+1-i))/det
    end do
    !
    ! if(mix_type==0)then
    !    if(iloop>1)calG(:,:) =  weight*calG(:,:) + (1.d0-weight)*Wold(:,is,:)
    !    Wold(:,is,:)  =  calG(:,:)
    ! endif
    !
    n0    = -sum(dimag(calG(1,:))*fermi(wr,beta))*fmesh/pi
    delta0= -u*sum(dimag(calG(2,:))*fermi(wr,beta))*fmesh/pi
    write(750,"(2I4,4(f16.12))",advance="yes")mpiID,is,n,n0,delta,delta0
    !
    sig_tmp(:,:)  =  solve_mpt_sc_sopt(calG,wr,n,n0,delta,delta0,Lreal)
    if(mix_type==1)sig_tmp(:,:) =  weight*sig_tmp(:,:) + (1.d0-weight)*Wold(:,is,:)
    !
  end subroutine ipt_solve_per_site_real






end module RDMFT_WRAP_IPT
