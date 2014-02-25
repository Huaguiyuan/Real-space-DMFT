!###############################################################
!     PURPOSE  : SOLVE DMFT-IPT REPULSIVE IN MATSUBARA FREQUENCY
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_MATS
  USE IPT_VARS_GLOBAL
  USE FFTGF
  implicit none
  private

  interface solve_ipt_sc_matsubara
     module procedure solve_ipt_sc_matsubara_r,solve_ipt_sc_matsubara_c
  end interface solve_ipt_sc_matsubara

  interface solve_mpt_sc_matsubara
     module procedure solve_mpt_sc_matsubara_r,solve_mpt_sc_matsubara_c
  end interface solve_mpt_sc_matsubara

  public :: solve_ipt_matsubara !half-filling
  public :: solve_ipt_sc_matsubara
  public :: solve_mpt_matsubara !away from h-f
  public :: solve_mpt_sc_matsubara

contains


  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order perturbation theory in Matsubara for 
  ! the repulsive model at half-filling
  !+-------------------------------------------------------------------+
  function solve_ipt_matsubara(fg0_iw) result(sigma_iw)
    complex(8),dimension(:)            :: fg0_iw
    complex(8),dimension(size(fg0_iw)) :: sigma_iw
    real(8),dimension(0:size(fg0_iw))  :: fg0_tau,sigma_tau
    integer                            :: i,Lf
    Lf=size(fg0_iw)
    call fftgf_iw2tau(fg0_iw,fg0_tau(0:),beta)
    forall(i=0:Lf)sigma_tau(i)=U**2*(fg0_tau(i))**2*fg0_tau(Lf-i)
    call fftgf_tau2iw(sigma_tau(0:),sigma_iw,beta)
    open(100,file="Sigma_tau.ipt")
    do i=0,Lf
       write(100,*)i*beta/dble(Lf),sigma_tau(i)
    enddo
    close(100)
  end function solve_ipt_matsubara





  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order modified perturbation theory in Matsubara for 
  ! the repulsive model away
  !+-------------------------------------------------------------------+
  function solve_mpt_matsubara(fg0_iw,n,n0,xmu0) result(sigma_iw)
    complex(8),dimension(:)            :: fg0_iw
    complex(8),dimension(size(fg0_iw)) :: sigma_iw
    real(8),dimension(0:size(fg0_iw))  :: fg0_tau,sigma_tau
    real(8)                            :: n,n0,xmu0
    real(8)                            :: A,B
    integer                            :: i,Lf
    Lf=size(fg0_iw)
    call fftgf_iw2tau(fg0_iw,fg0_tau(0:),beta)
    forall(i=0:Lf)sigma_tau(i)=U**2*(fg0_tau(i))**2*fg0_tau(Lf-i)
    call fftgf_tau2iw(sigma_tau(0:),sigma_iw,beta)
    call get_A
    call get_B
    sigma_iw = U*(n-0.5d0) + A*sigma_iw/(1.d0-B*sigma_iw)
    open(100,file="Sigma_tau.ipt")
    do i=0,Lf
       write(100,*)i*beta/dble(Lf),sigma_tau(i)
    enddo
    close(100)
  contains
    subroutine get_A
      real(8)                          :: A1,A2
      A1= n*(1.d0-n)
      A2= n0*(1.d0-n0)
      A = A1/A2
    end subroutine get_A
    subroutine get_B
      real(8)                          :: B1,B2
      B1 = (xmu0-xmu) + U*(1.d0-2.d0*n)
      B2 = n0*(1.d0-n0)*U**2
      B  = B1/B2
    end subroutine get_B
  end function solve_mpt_matsubara






  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order perturbation theory in Matsubara for 
  ! the attractive model at half-filling, REAL order parameter case
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_matsubara_r(fg0_iw,delta) result(sigma_iw)
    complex(8),dimension(:,:)                           :: fg0_iw
    complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
    real(8)                                             :: delta
    complex(8),dimension(:),allocatable                 :: calG11,calG22,calF
    real(8),dimension(:),allocatable                    :: calG11t,calG22t,calFt
    real(8),dimension(:),allocatable                    :: sigmat,selft
    integer                                             :: i,LM
    LM=size(fg0_iw,2)
    if(size(fg0_iw,1)/=2)stop "solve_ipt_sc_matsubara_r: size(input,1)!= 2"
    allocate(calG11(LM),calG11t(0:LM))
    allocate(calG22(LM),calG22t(0:LM))
    allocate(calF(LM),calFt(0:LM))
    allocate(sigmat(0:LM),selft(0:LM))
    !
    !Get the HF-corrected Weiss-Fields:
    calG11 =  fg0_iw(1,:)
    calG22 = -conjg(fg0_iw(1,:))
    calF   =  fg0_iw(2,:)
    call fftgf_iw2tau(calG11,calG11t(0:),beta)
    call fftgf_iw2tau(calG22,calG22t(0:),beta)
    call fftgf_iw2tau(calF,calFt(0:),beta,notail=.true.)
    !Get the 2nd-order Sigma:
    forall(i=0:LM)
       sigmat(i)=  U**2*(calG11t(i)*calG22t(i) - calFt(i)**2)*calG22t(LM-i)
       selft(i)= -U**2*(calFt(i)**2 - calG11t(i)*calG22t(i))*calFt(i)
    end forall
    call fftgf_tau2iw(sigmat(0:),sigma_iw(1,:),beta)
    call fftgf_tau2iw(selft(0:),sigma_iw(2,:),beta)
    !sigma_iw(1,:)=sigma_iw(1,:)
    sigma_iw(2,:)=sigma_iw(2,:) - delta
    !
    open(11,file="Sigma_tau.ipt",access="append")
    open(12,file="Self_tau.ipt",access="append")
    do i=0,LM
       write(11,*)i*beta/dble(LM),sigmat(i)
       write(12,*)i*beta/dble(LM),selft(i)
    enddo
    close(11);close(12)
    !
    deallocate(calG11,calG11t)
    deallocate(calG22,calG22t)
    deallocate(calF,calFt)
    deallocate(sigmat,selft)
  end function solve_ipt_sc_matsubara_r





  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order modified perturbation theory in Matsubara for 
  ! the attractice model away, REAL order parameter case
  !+-------------------------------------------------------------------+
  function solve_mpt_sc_matsubara_r(fg0_iw,n,n0,delta,delta0) result(sigma_iw)
    complex(8),dimension(:,:)                           :: fg0_iw
    complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
    complex(8),dimension(:),allocatable                 :: calG11,calG22,calF
    real(8),dimension(:),allocatable                    :: calG11t,calG22t,calFt
    real(8),dimension(:),allocatable                    :: sigmat,selft
    integer                                             :: i,LM
    real(8)                                             :: n,n0
    real(8)                                             :: delta,delta0
    real(8)                                             :: A,B
    LM=size(fg0_iw,2)
    if(size(fg0_iw,1)/=2)stop "solve_mipt_sc_matsubara_r: size(input,1)!= 2"
    allocate(calG11(LM),calG11t(0:LM))
    allocate(calG22(LM),calG22t(0:LM))
    allocate(calF(LM),calFt(0:LM))
    allocate(sigmat(0:LM),selft(0:LM))
    !
    !Get the HF-corrected Weiss-Fields:
    calG11 =  fg0_iw(1,:)
    calG22 = -conjg(fg0_iw(1,:))
    calF   =  fg0_iw(2,:)
    call fftgf_iw2tau(calG11,calG11t(0:),beta)
    call fftgf_iw2tau(calG22,calG22t(0:),beta)
    call fftgf_iw2tau(calF,calFt(0:),beta,notail=.true.)
    !Get the 2nd-order Sigma:
    forall(i=0:LM)
       sigmat(i)=  U**2*(calG11t(i)*calG22t(i) - calFt(i)**2)*calG22t(LM-i)
       selft(i)= -U**2*(calFt(i)**2 - calG11t(i)*calG22t(i))*calFt(i)
    end forall
    call fftgf_tau2iw(sigmat(0:),sigma_iw(1,:),beta)
    call fftgf_tau2iw(selft(0:),sigma_iw(2,:),beta)
    !
    A=U**2*n*(1.d0-n)-delta**2
    B=U**2*n0*(1.d0-n0)-delta0**2
    sigma_iw(1,:) =-U*(n-0.5d0) + sigma_iw(1,:)*A/B
    sigma_iw(2,:) =-delta       + sigma_iw(2,:)*A/B
    !
    open(11,file="Sigma_tau.ipt",access="append")
    open(12,file="Self_tau.ipt",access="append")
    do i=0,LM
       write(11,*)i*beta/dble(LM),sigmat(i)
       write(12,*)i*beta/dble(LM),selft(i)
    enddo
    close(11);close(12)
    !
    deallocate(calG11,calG11t)
    deallocate(calG22,calG22t)
    deallocate(calF,calFt)
    deallocate(sigmat,selft)
  end function solve_mpt_sc_matsubara_r







  !*******************************************************************
  !*******************************************************************
  !*******************************************************************






  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order perturbation theory in Matsubara for 
  ! the attractive model at half-filling, REAL order parameter case
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_matsubara_c(fg0_iw,delta) result(sigma_iw)
    complex(8),dimension(:,:)                           :: fg0_iw
    complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
    complex(8)                                          :: delta
    complex(8),dimension(:),allocatable                 :: calG11,calG22,calF
    real(8),dimension(:),allocatable                    :: calG11t,calG22t
    real(8),dimension(:),allocatable                    :: sigmat
    complex(8),dimension(:),allocatable                 :: calFt,selft
    integer                                             :: i,LM
    LM=size(fg0_iw,2)
    if(size(fg0_iw,1)/=2)stop "solve_ipt_sc_matsubara_c: size(input,1)!= 2"
    allocate(calG11(LM),calG11t(0:LM))
    allocate(calG22(LM),calG22t(0:LM))
    allocate(calF(LM),calFt(0:LM))
    allocate(sigmat(0:LM),selft(0:LM))
    !
    !Get the HF-corrected Weiss-Fields:
    calG11 =  fg0_iw(1,:)
    calG22 = -conjg(fg0_iw(1,:))
    calF   =  fg0_iw(2,:)
    call fftgf_iw2tau(calG11,calG11t(0:),beta)
    call fftgf_iw2tau(calG22,calG22t(0:),beta)
    call fftff_iw2tau(calF,calFt(0:),beta)
    !Get the 2nd-order Sigma:
    forall(i=0:LM)
       sigmat(i)=  U**2*(calG11t(i)*calG22t(i) - abs(calFt(i))**2)*calG22t(LM-i)
       selft(i) = -U**2*(abs(calFt(i))**2 - calG11t(i)*calG22t(i))*calFt(i) !ACTHUNG HERE: time inversion simmetry
    end forall
    call fftgf_tau2iw(sigmat(0:),sigma_iw(1,:),beta)
    call fftff_tau2iw(selft(0:),sigma_iw(2,:),beta)
    !
    sigma_iw(2,:)=sigma_iw(2,:) - delta
    !
    open(11,file="Sigma_tau.ipt",access="append")
    open(12,file="Self_tau.ipt",access="append")
    do i=0,LM
       write(11,*)i*beta/dble(LM),sigmat(i)
       write(12,*)i*beta/dble(LM),dreal(selft(i)),dimag(selft(i))
    enddo
    close(11);close(12)
    !
    deallocate(calG11,calG11t)
    deallocate(calG22,calG22t)
    deallocate(calF,calFt)
    deallocate(sigmat,selft)
  end function solve_ipt_sc_matsubara_c



  !+-------------------------------------------------------------------+
  !PURPOSE: Solve 2nd order modified perturbation theory in Matsubara for 
  ! the attractice model away, REAL order parameter case
  !+-------------------------------------------------------------------+
  function solve_mpt_sc_matsubara_c(fg0_iw,n,n0,delta,delta0) result(sigma_iw)
    complex(8),dimension(:,:)                           :: fg0_iw
    complex(8),dimension(size(fg0_iw,1),size(fg0_iw,2)) :: sigma_iw
    complex(8),dimension(:),allocatable                 :: calG11,calG22,calF
    real(8),dimension(:),allocatable                    :: calG11t,calG22t
    real(8),dimension(:),allocatable                    :: sigmat
    complex(8),dimension(:),allocatable                 :: calFt,selft
    integer                                             :: i,LM
    real(8)                                             :: n,n0
    complex(8)                                          :: delta,delta0
    real(8)                                             :: A,B
    LM=size(fg0_iw,2)
    if(size(fg0_iw,1)/=2)stop "solve_ipt_sc_matsubara_c: size(input,1)!= 2"
    allocate(calG11(LM),calG11t(0:LM))
    allocate(calG22(LM),calG22t(0:LM))
    allocate(calF(LM),calFt(0:LM))
    allocate(sigmat(0:LM),selft(0:LM))
    !
    !Get the HF-corrected Weiss-Fields:
    calG11 =  fg0_iw(1,:)
    calG22 = -conjg(fg0_iw(1,:))
    calF   =  fg0_iw(2,:)
    call fftgf_iw2tau(calG11,calG11t(0:),beta)
    call fftgf_iw2tau(calG22,calG22t(0:),beta)
    call fftff_iw2tau(calF,calFt(0:),beta)
    !Get the 2nd-order Sigma:
    forall(i=0:LM)
       sigmat(i)=  U**2*(calG11t(i)*calG22t(i) - abs(calFt(i))**2)*calG22t(LM-i)
       selft(i) = -U**2*(abs(calFt(i))**2 - calG11t(i)*calG22t(i))*calFt(i) !ACTHUNG HERE: time inversion simmetry
    end forall
    call fftgf_tau2iw(sigmat(0:),sigma_iw(1,:),beta)
    call fftff_tau2iw(selft(0:),sigma_iw(2,:),beta)
    !
    !This is not obvious, just a guess now but I need to check!!
    A=U**2*n*(1.d0-n)-abs(delta)**2
    B=U**2*n0*(1.d0-n0)-abs(delta0)**2
    sigma_iw(1,:) =-U*(n-0.5d0) + sigma_iw(1,:)*A/B
    sigma_iw(2,:) =-delta       + sigma_iw(2,:)*A/B
    !
    open(11,file="Sigma_tau.ipt",access="append")
    open(12,file="Self_tau.ipt",access="append")
    do i=0,LM
       write(11,*)i*beta/dble(LM),sigmat(i)
       write(12,*)i*beta/dble(LM),dreal(selft(i)),dimag(selft(i))
    enddo
    close(11);close(12)
    !
    deallocate(calG11,calG11t)
    deallocate(calG22,calG22t)
    deallocate(calF,calFt)
    deallocate(sigmat,selft)
  end function solve_mpt_sc_matsubara_c







end module IPT_MATS
