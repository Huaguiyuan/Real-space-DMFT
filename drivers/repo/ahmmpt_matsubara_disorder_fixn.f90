!########################################################
!     PURPOSE  :solve the attractive (A) disordered (D) Hubbard
! model (HM) using Modified Perturbation THeory (MPT), w/ DMFT
!     COMMENTS : DISORDER REALIZATION DEPENDS ON THE PARAMETER int idum:
!SO THAT DIFFERENT REALIZATIONS (STATISTICS) ARE PERFORMED 
!CALLING THIS PROGRAM MANY TIMES WHILE PROVIDING A *DIFFERENT* SEED.
!THE RESULT OF EACH CALCULATION IS STORED IN A DIFFERENT DIR
!     AUTHORS  : A.Amaricci
!########################################################
program adhmipt_matsubara
  USE RDMFT_VARS_GLOBAL
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma,gf_tmp
  complex(8),allocatable,dimension(:,:)   :: fg0
  logical                                 :: converged
  real(8)                                 :: r,delta,delta0,naverage
  integer                                 :: i,is,esp,lm
  ! real(8) :: ndelta,ndelta1
  ! integer :: nindex,nindex1

  !GLOBAL INITIALIZATION:
  !===============================================================
  include "init_global_disorder.f90"

  !ALLOCATE WORKING ARRAYS
  !===============================================================
  allocate(fg(2,Ns,L))
  allocate(fg0(2,L))
  allocate(sigma(2,Ns,L))
  allocate(gf_tmp(2,Ns,L))

  nread=0.875d0
  ndelta=0.1d0
  nerror=1.d-4

  !START DMFT LOOP SEQUENCE:
  !==============================================================
  call setup_initial_sc_sigma()
  iloop=0 ; converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     call get_sc_gloc_mpi()      !SOLVE G_II (GLOCAL) \FORALL FREQUENCIES
     call solve_sc_impurity_mpi()!SOLVE IMPURITY MODEL, \FORALL LATTICE SITES:
     converged=check_convergence(sigma(1,:,:)+sigma(2,:,:),eps_error,Nsuccess,nloop,id=0)
     naverage=get_naverage()
     if(mpiID==0)then
        nindex1=nindex
        ndelta1=ndelta
        if((naverage.ge.nread+nerror))then
           nindex=-1
        elseif(naverage.le.nread-nerror)then
           nindex=1
        else
           nindex=0
        endif
        if(nindex1+nindex.eq.0)then !avoid loop forth and back
           ndelta=ndelta1/2.d0 !decreasing the step
           xmu=xmu+real(nindex,8)*ndelta
        else
           ndelta=ndelta1
           xmu=xmu+dble(nindex)*ndelta
        endif
        write(*,"(A,f15.12,A,f15.12,A,f15.12)")" n=",naverage,"/",nread,"| ",ndelta
        print*,""
        print*,abs(naverage-nread),xmu
        print*,""        
        if(abs(naverage-nread)>1.d-5)converged=.false.
        call splot(trim(adjustl(trim(name_dir)))//"/muVSiter.ipt",iloop,xmu,abs(naverage-nread),append=.true.)
     endif
     call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_sc_out(converged)
     call end_loop()
  enddo
  if(mpiID==0)call system("mv -vf *.err "//trim(adjustl(trim(name_dir)))//"/")
  call close_mpi()



contains


  subroutine setup_initial_sc_sigma()
    logical :: check1,check2,check
    inquire(file="LSigma.ipt",exist=check1)
    inquire(file="LSelf.ipt",exist=check2)
    check=check1*check2
    if(check)then
       if(mpiID==0)then
          write(*,*)"Reading Sigma in input:"
          call sread("LSigma.ipt",sigma(1,1:Ns,1:L))
          call sread("LSelf.ipt",sigma(2,1:Ns,1:L))
       endif
       call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    else
       n=0.5d0 ; delta=deltasc
       sigma(1,:,:)=zero ; sigma(2,:,:)=-delta
    endif
  end subroutine setup_initial_sc_sigma


  !******************************************************************
  !******************************************************************


  subroutine get_sc_gloc_mpi()
    complex(8) :: Gloc(2*Ns,2*Ns)
    integer    :: i,is
    call msg("Get local GF: (ETA --> fort.999)",id=0)
    call start_timer
    fg=zero ; gf_tmp=zero
    do i=1+mpiID,L,mpiSIZE
       Gloc=zero
       Gloc(1:Ns,1:Ns)          = -H0
       Gloc(Ns+1:2*Ns,Ns+1:2*Ns)=  H0
       do is=1,Ns
          Gloc(is,is)      =  xi*wm(i)+xmu-sigma(1,is,i)-erandom(is)
          Gloc(Ns+is,Ns+is)= -conjg(Gloc(is,is))
          Gloc(is,Ns+is)   = -sigma(2,is,i)
          Gloc(Ns+is,is)   = -sigma(2,is,-i)
       enddo
       call mat_inversion_sym(Gloc,2*Ns)
       forall(is=1:Ns)
          gf_tmp(1,is,i) = Gloc(is,is)
          gf_tmp(2,is,i) = Gloc(is,Ns+is)
       end forall
       call eta(i,L,999)
    enddo
    call stop_timer
    call MPI_REDUCE(gf_tmp,fg,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(fg,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_sc_gloc_mpi


  !******************************************************************
  !******************************************************************



  subroutine solve_sc_impurity_mpi()
    integer    :: is
    complex(8) :: zsigma(2,Ns,L)
    call msg("Solve impurity: (ETA --> fort.998)")
    if(Wdis/=0.d0)then
       call start_timer
       zsigma=zero ; sigma=zero
       do is=1+mpiID,Ns,mpiSIZE
          call solve_per_site(is)
          call eta(is,Ns,998)
       enddo
       call stop_timer
       call MPI_REDUCE(sigma,zsigma,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(zsigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       sigma=zsigma
    else
       call solve_per_site(is=1)
       forall(is=2:Ns)sigma(:,is,:)=sigma(:,1,:)
    end if
  end subroutine solve_sc_impurity_mpi



  !******************************************************************
  !******************************************************************



  subroutine solve_per_site(is)
    integer                                      :: is
    complex(8)                                   :: det(L)
    complex(8),dimension(:,:,:),allocatable,save :: sold
    complex(8),dimension(2,L)                    :: calG
    real(8),dimension(2,0:L)                     :: fgt,fg0t

    if(.not.allocated(sold))allocate(sold(2,Ns,L))
    sold(:,is,:)  =  sigma(:,is,:)

    call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
    call fftgf_iw2tau(fg(2,is,:),fgt(2,0:L),beta,notail=.true.)
    n    = -real(fgt(1,L),8) ; delta= -u*fgt(2,0)

    calG=zero ; fg0=zero
    det       = abs(fg(1,is,:))**2    + fg(2,is,:)**2
    fg0(1,:)  = conjg(fg(1,is,:))/det + sigma(1,is,:) - U*(n-0.5d0)
    fg0(2,:)  = fg(2,is,:)/det        + sigma(2,is,:) + delta

    det       =  abs(fg0(1,:))**2 + fg0(2,:)**2
    calG(1,:) =  conjg(fg0(1,:))/det
    calG(2,:) =  fg0(2,:)/det

    call fftgf_iw2tau(calG(1,:),fg0t(1,:),beta)
    call fftgf_iw2tau(calG(2,:),fg0t(2,:),beta,notail=.true.)
    n0=-real(fg0t(1,L)) ; delta0= -u*fg0t(2,0)
    write(750,"(I4,4(f16.12))",advance="yes")is,n,n0,delta,delta0
    sigma(:,is,:) =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)
    sigma(:,is,:) =  weigth*sigma(:,is,:) + (1.d0-weigth)*sold(:,is,:)
  end subroutine solve_per_site



  !******************************************************************
  !******************************************************************


  function get_naverage result(naverage)
    integer :: is
    real(8) :: naverage
    real(8) :: nii(Ns)
    real(8),dimension(2,0:L) :: fgt
    if(mpiID==0)then
       naverage=0.d0
       do is=1,Ns
          call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
          nii(is) = -2.d0*real(fgt(1,L))
       enddo
       naverage = sum(nii)/dble(Ns)
    endif
    call MPI_BCAST(naverage,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end function get_naverage


  !******************************************************************
  !******************************************************************


  subroutine print_sc_out(converged)
    integer                  :: is
    real(8)                  :: nimp,delta
    complex(8)               :: afg(2,L),asig(2,L)
    real(8)                  :: nii(Ns),dii(Ns)
    logical                  :: converged
    character(len=4)         :: loop
    real(8),dimension(2,0:L) :: fgt

    if(mpiID==0)then
       nimp=0.d0 ; delta=0.d0
       do is=1,Ns
          call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
          call fftgf_iw2tau(fg(2,is,:),fgt(2,0:L),beta,notail=.true.)
          nii(is) = -2.d0*real(fgt(1,L))
          dii(is) = -u*fgt(2,0)
       enddo
       nimp = sum(nii)/dble(Ns)
       delta= sum(dii)/dble(Ns)
       print*,"nimp  =",nimp
       print*,"delta =",delta
       call splot(trim(adjustl(trim(name_dir)))//"/navVSiloop.ipt",iloop,nimp,append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/davVSiloop.ipt",iloop,delta,append=TT)

       ! write(loop,"(I4)")iloop
       ! do i=1,Ns
       !    call splot("Gloc_iw_site."//trim(adjustl(trim(loop)))//".ipt",wm,fg(1,i,1:L),append=TT)
       !    call splot("Floc_iw_site."//trim(adjustl(trim(loop)))//".ipt",wm,fg(2,i,1:L),append=TT)       
       !    call splot("Sigma_iw_site."//trim(adjustl(trim(loop)))//".ipt",wm,sigma(1,i,1:L),append=TT)
       !    call splot("Self_iw_site."//trim(adjustl(trim(loop)))//".ipt",wm,sigma(2,i,1:L),append=TT)
       ! enddo


       if(converged)then
          call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.ipt",nii)
          call splot(trim(adjustl(trim(name_dir)))//"/deltaVSisite.ipt",dii)
          call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.ipt",erandom)
          call splot(trim(adjustl(trim(name_dir)))//"/LSigma.ipt",sigma(1,1:Ns,1:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LSelf.ipt",sigma(2,1:Ns,1:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LG.ipt",fg(1,1:Ns,1:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LF.ipt",fg(2,1:Ns,1:L))
       endif
    end if
    return
  end subroutine print_sc_out





end program adhmipt_matsubara
