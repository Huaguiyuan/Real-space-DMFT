!########################################################
!     Program  : ADHMMPT
!     TYPE     : Main program
!     PURPOSE  :solve the attractive (A) disordered (D) Hubbard
! model (HM) using Modified Perturbation THeory (MPT), w/ DMFT
!     COMMENTS :
!DISORDER REALIZATION DEPENDS ON THE PARAMETER int idum:
!SO THAT DIFFERENT REALIZATIONS (STATISTICS) ARE PERFORMED 
!CALLING THIS PROGRAM MANY TIMES WHILE PROVIDING A *DIFFERENT* SEED.
!THE RESULT OF EACH CALCULATION IS STORED IN A DIFFERENT DIR
!     AUTHORS  : A.Amaricci
!########################################################
program adhmipt
  !LOCAL: 
  USE RDMFT_VARS_GLOBAL
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: fg,fg0,sigma,gf_tmp
  real(8)                                 :: r
  logical                                 :: converged
  real(8)                                 :: delta,delta0
  integer :: is
  complex(8) :: zeta

  print*," ###############################################################"
  print*," This code got a problem:"
  print*," it looks like it works, but check out the self-energies!"
  print*," these have a weird behavior, I suspect it comes from the"
  print*," new version of kramers-kronig integral. but I can not exclude"
  print*," it is an effect of a wrong MPI copy."
  print*," ###############################################################"
  print*," ...continue..."
  call sleep(5)

  !GLOBAL INITIALIZATION:
  !===============================================================
  include "init_global.h"


  !ALLOCATE WORKING ARRAYS
  !===============================================================
  allocate(fg(2,Ns,-L:L))
  allocate(fg0(2,Ns,-L:L))
  allocate(sigma(2,Ns,-L:L))
  allocate(gf_tmp(2,Ns,-L:L))


  !START DMFT LOOP SEQUENCE:
  !===============================================================================
  call setup_initial_sc_sigma()
  iloop=0 ; converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     if(mpiID==0)write(*,"(A,I5,A,L6)")"DMFT-loop",iloop," convergence:",converged
     call get_sc_gloc_mpi()      !SOLVE G_II (GLOCAL) \FORALL FREQUENCIES
     call solve_sc_impurity_mpi()!SOLVE IMPURITY MODEL, \FORALL LATTICE SITES:
     converged=check_convergence(sigma(1,:,:)+sigma(2,:,:),eps_error,Nsuccess,nloop,id=0)
     call print_sc_out(converged)
     call msg("============================================",lines=2)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
  enddo
  if(mpiID==0)call system("mv -vf *.err "//trim(adjustl(trim(name_dir)))//"/")
  call close_mpi()


contains


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
          call sread("LSigma.ipt",sigma(1,1:Ns,-L:L))
          call sread("LSelf.ipt",sigma(2,1:Ns,-L:L))
       endif
       call MPI_BCAST(sigma,2*Ns*(2*L+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    else
       n=0.5d0 ; delta=deltasc
       sigma(1,:,:)=zero ; sigma(2,:,:)=-delta
    endif
  end subroutine setup_initial_sc_sigma


  !******************************************************************
  !******************************************************************


  subroutine get_sc_gloc_mpi()
    complex(8) :: Gloc(2*Ns,2*Ns),zeta1,zeta2
    integer    :: i,is
    call msg("Get local GF:",id=0)
    call start_timer
    fg=zero ; gf_tmp=zero

    !(w + ie + mu%ran - H) - Sigma(w)   :: -(-w + ie + mu%ran - H - Sigma(-w))*
    ! -(-w -ie + mu%ran -H -Sigma*(-w)) ::    w + ie - mu%ran + H + Sigma*(-w)
    do i=-L+mpiID,L,mpiSIZE
       Gloc=zero
       Gloc(1:Ns,1:Ns)          = -H0
       Gloc(Ns+1:2*Ns,Ns+1:2*Ns)=  H0
       do is=1,Ns
          zeta1=cmplx(wr(i),eps,8)  + xmu - sigma(1,is,i)         - erandom(is)
          zeta2=cmplx(wr(i),eps,8)  + xmu + conjg(sigma(1,is,-i)) + erandom(is)
          Gloc(is,is)      = Gloc(is,is)       + zeta1
          Gloc(Ns+is,Ns+is)= Gloc(Ns+is,Ns+is) + zeta2
          Gloc(is,Ns+is)   = -sigma(2,is,i)
          Gloc(Ns+is,is)   =  sigma(2,is,-i)
       enddo
       call mat_inversion_sym(Gloc,2*Ns)
       forall(is=1:Ns)
          gf_tmp(1,is,i) = Gloc(is,is)
          gf_tmp(2,is,i) = Gloc(is,Ns+is)
       end forall
       call eta(L+i+1,2*L+1,999)
    enddo
    call stop_timer
    call MPI_REDUCE(gf_tmp,fg,2*Ns*(2*L+1),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(fg,2*Ns*(2*L+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_sc_gloc_mpi


  !******************************************************************
  !******************************************************************



  subroutine solve_sc_impurity_mpi()
    integer    :: is
    complex(8) :: zsigma(2,Ns,-L:L)
    call msg("Solve impurity:")
    if(Wdis/=0.d0)then
       call start_timer
       zsigma=zero ; sigma=zero
       do is=1+mpiID,Ns,mpiSIZE
          call solve_per_site(is)
          call eta(is,Ns,998)
       enddo
       call stop_timer
       call MPI_REDUCE(sigma,zsigma,2*Ns*(2*L+1),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(zsigma,2*Ns*(2*L+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
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
    complex(8)                                   :: det
    complex(8),dimension(:,:,:),allocatable,save :: sold
    complex(8),dimension(2,-L:L)                 :: calG
    if(.not.allocated(sold))then
       allocate(sold(2,Ns,-L:L))
       sold=sigma
    endif
    n    = -sum(fermi(wr,beta)*dimag(fg(1,is,-L:L)))*fmesh/pi!/sum(dimag(fg(1,is,-L:L)))
    delta= -u*sum(dimag(fg(2,is,:))*fermi(wr,beta))*fmesh/pi

    calG=zero ; fg0(:,is,:)=zero
    do i=-L,L
       det       = fg(1,is,i)*conjg(fg(1,is,-i)) - fg(2,is,i)*fg(2,is,-i)
       calG(1,i) = conjg(fg(1,is,-i))/det + sigma(1,is,i) - U*(n-0.5d0)
       calG(2,i) = -fg(2,is,-i)/det       + sigma(2,is,i) + delta
    end do
    do i=-L,L
       det         = calG(1,i)*conjg(calG(1,-i)) - calG(2,i)*calG(2,-i)
       fg0(1,is,i) = conjg(calG(1,-i))/det
       fg0(2,is,i) = -calG(2,-i)/det
    end do
    n0    = -sum(dimag(fg0(1,is,:))*fermi(wr,beta))*fmesh/pi!/sum(dimag(fg(1,is,:)))
    delta0= -u*sum(dimag(fg0(2,is,:))*fermi(wr,beta))*fmesh/pi
    write(*,"(I4,4(f16.12))",advance="yes")is,n,n0,delta,delta0

    sigma(1:2,is,-L:L)             =  solve_mpt_sc_sopt(fg0(1:2,is,-L:L),n,n0,delta,delta0)
    sigma(1:2,is,-L:L)             =  weigth*sigma(1:2,is,-L:L) + (1.d0-weigth)*sold(1:2,is,-L:L)
    sold(:,is,:)                   =  sigma(:,is,:)
  end subroutine solve_per_site



  !******************************************************************
  !******************************************************************



  subroutine print_sc_out(converged)
    integer   :: is
    real(8)   :: nimp,delta
    complex(8):: afg(2,-L:L),asig(2,-L:L)
    real(8)   :: nii(Ns),dii(Ns)
    logical   :: converged
    character(len=4) :: loop

    if(mpiID==0)then
       nimp=0.d0
       delta=0.d0
       do is=1,Ns
          nii(is) = 2.d0*sum(fermi(wr,beta)*dimag(fg(1,is,-L:L)))/sum(dimag(fg(1,is,-L:L)))
          dii(is) = -u*sum(dimag(fg(2,is,-L:L))*fermi(wr(-L:L),beta))*fmesh/pi
       enddo
       nimp = sum(nii)/dble(Ns)
       delta= sum(dii)/dble(Ns)
       print*,"nimp  =",nimp
       print*,"delta =",delta
       call splot(trim(adjustl(trim(name_dir)))//"/navVSiloop.ipt",iloop,nimp,append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/davVSiloop.ipt",iloop,delta,append=TT)

       write(loop,"(I4)")iloop
       do i=1,Ns
          call splot("Gloc_site."//trim(adjustl(trim(loop)))//".ipt",wr,fg(1,i,-L:L),append=TT)
          call splot("Floc_site."//trim(adjustl(trim(loop)))//".ipt",wr,fg(2,i,-L:L),append=TT)       
          call splot("Sigma_site."//trim(adjustl(trim(loop)))//".ipt",wr,sigma(1,i,-L:L),append=TT)
          !call splot("G0_site.ipt",wr,fg0(1,i,-L:L),append=TT)
          call splot("Self_site."//trim(adjustl(trim(loop)))//".ipt",wr,sigma(2,i,-L:L),append=TT)
          !call splot("F0_site.ipt",wr,fg0(2,i,-L:L),append=TT)
       enddo


       if(converged)then
          call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.ipt",nii)
          call splot(trim(adjustl(trim(name_dir)))//"/deltaVSisite.ipt",dii)
          call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.ipt",erandom)
          call splot(trim(adjustl(trim(name_dir)))//"/LSigma.ipt",sigma(1,1:Ns,-L:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LSelf.ipt",sigma(2,1:Ns,-L:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LG.ipt",fg(1,1:Ns,-L:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LF.ipt",fg(2,1:Ns,-L:L))
          afg(:,:)   =sum(fg(1:2,1:Ns,-L:L),dim=2)/dble(Ns)
          call splot(trim(adjustl(trim(name_dir)))//"/DOS.disorder.ipt",wr,-dimag(afg(1,-L:L))/pi,append=TT)
       endif
    end if
    return
  end subroutine print_sc_out


  !******************************************************************
  !******************************************************************

end program adhmipt
