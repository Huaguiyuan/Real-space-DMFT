!########################################################
!     PURPOSE  :solve the attractive (A) Trapped (T) Hubbard
! model (HM) using Modified Perturbation THeory (MPT), w/ DMFT
!     COMMENTS : 
!     AUTHORS  : A.Amaricci
!########################################################
program adhmipt_matsubara_trap
  USE RDMFT_VARS_GLOBAL
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma,gf_tmp
  complex(8),allocatable,dimension(:,:)   :: fg0
  logical                                 :: converged
  real(8)                                 :: r,delta,delta0
  integer                                 :: is,esp,lm


  !GLOBAL INITIALIZATION:
  !===============================================================
  include "init_global_trap.f90"


  !ALLOCATE WORKING ARRAYS
  !===============================================================
  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)
  write(*,"(A,I9,A)")"Using ",L," frequencies"

  allocate(fg(2,Ns,L))
  allocate(fg0(2,L))
  allocate(sigma(2,Ns,L))
  allocate(gf_tmp(2,Ns,L))


  !START DMFT LOOP SEQUENCE:
  !==============================================================
  call setup_initial_sc_sigma()
  iloop=0 ; converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     if(mpiID==0)write(*,"(A,I5,A,L6)")"DMFT-loop",iloop," convergence:",converged
     call get_sc_gloc_mpi()      !SOLVE G_II (GLOCAL) \FORALL FREQUENCIES:
     call solve_sc_impurity_mpi()!SOLVE IMPURITY MODEL,\FORALL LATTICE SITES:
     converged=check_convergence(sigma(1,:,:)+sigma(2,:,:),eps_error,Nsuccess,nloop,id=0)
     call print_sc_out(converged)
     call msg("============================================",lines=2)
  enddo
  if(mpiID==0)call system("mv -vf *.err "//trim(adjustl(trim(name_dir)))//"/")
  call close_mpi()

contains


  pure function trap_distance_square(is) result(dist)
    integer,intent(in) :: is
    integer            :: center
    real(8)            :: dist
    center=(Ns+1)/2
    dist=dble(irow(is)-irow(center))**2 + dble(icol(is)-icol(center))**2
  end function trap_distance_square



  !******************************************************************
  !******************************************************************



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
          Gloc(is,is)      =  xi*wm(i)+xmu-sigma(1,is,i)-etrap(is)
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
    call fft_iw2tau(fg(2,is,:),fgt(2,0:L),beta,L)
    n    = -real(fgt(1,L),8) ; delta= -u*fgt(2,0)

    calG=zero ; fg0=zero
    det       = abs(fg(1,is,:))**2    + fg(2,is,:)**2
    fg0(1,:)  = conjg(fg(1,is,:))/det + sigma(1,is,:) - U*(n-0.5d0)
    fg0(2,:)  = fg(2,is,:)/det        + sigma(2,is,:) + delta

    det       =  abs(fg0(1,:))**2 + fg0(2,:)**2
    calG(1,:) =  conjg(fg0(1,:))/det
    calG(2,:) =  fg0(2,:)/det

    call fftgf_iw2tau(calG(1,:),fg0t(1,:),beta)
    call fft_iw2tau(calG(2,:),fg0t(2,:),beta,L)
    n0=-real(fg0t(1,L)) ; delta0= -u*fg0t(2,0)
    write(750,"(I4,4(f16.12))",advance="yes")is,n,n0,delta,delta0
    sigma(:,is,:) =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)
    sigma(:,is,:) =  weigth*sigma(:,is,:) + (1.d0-weigth)*sold(:,is,:)
  end subroutine solve_per_site



  !******************************************************************
  !******************************************************************



  subroutine print_sc_out(converged)
    integer                  :: is,row,col
    real(8)                  :: nimp,delta
    complex(8)               :: afg(2,L),asig(2,L)
    real(8)                  :: nii(Ns),dii(Ns)
    logical                  :: converged
    character(len=4)         :: loop
    real(8),dimension(2,0:L) :: fgt
    real(8),allocatable                     :: grid_x(:),grid_y(:)
    real(8),allocatable                     :: dij(:,:),nij(:,:),eij(:,:)

    if(mpiID==0)then
       nimp=0.d0 ; delta=0.d0
       do is=1,Ns
          call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
          call fft_iw2tau(fg(2,is,:),fgt(2,0:L),beta,L)
          nii(is) = -2.d0*real(fgt(1,L))
          dii(is) = -u*fgt(2,0)
       enddo
       nimp = sum(nii)/dble(Ns)
       delta= sum(dii)/dble(Ns)
       print*,"nimp  =",nimp
       print*,"delta =",delta
       call splot(trim(adjustl(trim(name_dir)))//"/navVSiloop.ipt",iloop,nimp,append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/davVSiloop.ipt",iloop,delta,append=TT)

       call system("rm -fv Gloc_iw*ipt Floc_iw*ipt Sigma_iw*ipt Self_iw*ipt")
       do i=1,Ns
          call splot("Gloc_iw_site.ipt",wm,fg(1,i,1:L),append=TT)
          call splot("Floc_iw_site.ipt",wm,fg(2,i,1:L),append=TT)
          call splot("Sigma_iw_site.ipt",wm,sigma(1,i,1:L),append=TT)
          call splot("Self_iw_site.ipt",wm,sigma(2,i,1:L),append=TT)
       enddo


       if(converged)then
          call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.ipt",nii)
          call splot(trim(adjustl(trim(name_dir)))//"/deltaVSisite.ipt",dii)
          call splot(trim(adjustl(trim(name_dir)))//"/etrapVSisite.ipt",etrap)
          call splot(trim(adjustl(trim(name_dir)))//"/LSigma.ipt",sigma(1,1:Ns,1:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LSelf.ipt",sigma(2,1:Ns,1:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LG.ipt",fg(1,1:Ns,1:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LF.ipt",fg(2,1:Ns,1:L))
          ! afg(:,:)   =sum(fg(1:2,1:Ns,1:L),dim=2)/dble(Ns)
          ! asig(:,:)  =sum(sigma(1:2,1:Ns,1:L),dim=2)/dble(Ns)
          ! call splot(trim(adjustl(trim(name_dir)))//"/Gave_iw.ipt",wm,afg(1,1:L))
          ! call splot(trim(adjustl(trim(name_dir)))//"/Sigmaave_iw.ipt",wm,asig(1,1:L))


          !BUILD A GRID FOR  LATTICE PLOTS:
          allocate(grid_x(Nside),grid_y(Nside))
          allocate(nij(Nside,Nside),dij(Nside,Nside),eij(Nside,Nside))
          do row=1,Nside
             grid_x(row)=dble(row)
             grid_y(row)=dble(row)
          enddo

          do row=1,Nside
             do col=1,Nside
                i            = ij2site(row,col)
                nij(row,col) = nii(i)
                dij(row,col) = dii(i)
                eij(row,col) = etrap(i)
             enddo
          enddo
          call splot("3d_nVSij.data",grid_x,grid_y,nij)
          call splot("3d_deltaVSij.data",grid_x,grid_y,dij)
          call splot("3d_etrapVSij.data",grid_x,grid_y,eij)
          call system("mv -vf 3d_*.data plot_3d* "//trim(adjustl(trim(name_dir)))//"/ ")
       endif
    end if
    return
  end subroutine print_sc_out



end program adhmipt_matsubara_trap
