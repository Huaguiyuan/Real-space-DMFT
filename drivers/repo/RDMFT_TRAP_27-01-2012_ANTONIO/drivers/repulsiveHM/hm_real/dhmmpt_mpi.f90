  !########################################################
  !     Program  : DHMMPT
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
  include "broydn_common.h"
  program hmmpt
    !LOCAL:
    USE RDMFT_VARS_GLOBAL
    USE COMMON_BROYDN
    implicit none
    integer    :: is
    real(8)    :: x(1),r
    logical    :: check,converged  


    !GLOBAL INITIALIZATION:
    !=====================================================================
    include "init_global.h"


    !ALLOCATE WORKING ARRAYS:
    !=====================================================================
    allocate(fg(Ns,-L:L),sigma(Ns,-L:L))
    allocate(fg0(-L:L),gamma(-L:L))


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


    !******************************************************************
    !******************************************************************


    subroutine setup_initial_sigma()
      logical :: check1,check2,check
      inquire(file="LSigma.ipt",exist=check)
      if(check)then
         if(mpiID==0)then
            write(*,*)"Reading Sigma in input:"
            call sread("LSigma.ipt",sigma(1:Ns,-L:L))
         endif
         call MPI_BCAST(sigma,Ns*(2*L+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
      else
         sigma=u*(n-0.5d0)
      endif
    end subroutine setup_initial_sigma



    !******************************************************************
    !******************************************************************



    subroutine get_gloc_mpi() 
      complex(8) :: zeta,Gloc(Ns,Ns),gf_tmp(Ns,-L:L)
      integer    :: i,is
      call msg("Get local GF:",id=0)
      call start_timer
      gf_tmp=zero ; fg=zero
      do i=-L+mpiID,L,mpiSIZE
         Gloc=zero
         zeta  = cmplx(wr(i),eps,8) + xmu
         Gloc  = -H0
         do is=1,Ns
            Gloc(is,is) = zeta - sigma(is,i) - erandom(is)
         enddo
         call mat_inversion_sym(Gloc,Ns)
         do is=1,Ns
            gf_tmp(is,i) = Gloc(is,is)
         enddo
         call eta(L+i+1,(2*L+1),999)
      enddo
      call stop_timer
      call MPI_REDUCE(gf_tmp,fg,Ns*(2*L+1),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
      call MPI_BCAST(fg,Ns*(2*L+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
      call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
    end subroutine get_gloc_mpi


    !******************************************************************
    !******************************************************************


    subroutine solve_impurity_mpi()
      integer    :: is
      complex(8) :: zsigma(Ns,-L:L)
      call msg("Solve impurity:")
      if(Wdis/=0.d0)then
         call start_timer
         zsigma=zero ; sigma=zero !sigma must be set to zero or reduction will ends up in a sum w/ old sigmas
         do is=1+mpiID,Ns,mpiSIZE
            call solve_per_site(is)
            call eta(is,Ns,998)
         enddo
         call stop_timer
         call MPI_REDUCE(sigma,zsigma,Ns*(2*L+1),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
         call MPI_BCAST(zsigma,Ns*(2*L+1),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
         sigma=zsigma
      else
         call solve_per_site(is=1)
         forall(is=2:Ns)sigma(is,:)=sigma(1,:)
      endif
    end subroutine solve_impurity_mpi


    !******************************************************************
    !******************************************************************


    subroutine solve_per_site(is)
      integer                     :: is
      complex(8),allocatable,save :: Sold(:,:)
      isite=is
      if(.not.allocated(Sold))then
         allocate(Sold(Ns,-L:L))
         Sold=Sigma
      endif
      n=sum(fermi(wr,beta)*dimag(fg(is,:)))/sum(dimag(fg(is,:)))
      gamma = one/(one/fg(is,:) + sigma(is,:))
      xmu0=0.d0
      x(1)=xmu
      call broydn(x,check)
      xmu0=x(1)
      sigma(is,:) =  solve_mpt_sopt(fg0,n,n0,xmu0)
      sigma(is,:) =  weigth*sigma(is,:) + (1.d0-weigth)*sold(is,:)
      sold(is,:)  =  sigma(is,:)
    end subroutine solve_per_site


    !******************************************************************
    !******************************************************************


    subroutine print_out(success)
      real(8)   :: nimp,nii(Ns)
      complex(8):: afg(-L:L),asig(-L:L)
      logical   :: success
      integer   :: i,is
      character(len=4) :: loop
      if(mpiID==0)then
         nimp=0.d0
         do is=1,Ns
            nimp = nimp + 2.d0*sum(fermi(wr,beta)*dimag(fg(is,-L:L)))/sum(dimag(fg(is,-L:L)))/dble(Ns)
         enddo

         print*,"nimp  =",nimp
         call splot(trim(adjustl(trim(name_dir)))//"/navVSiloop.ipt",iloop,nimp,append=TT)


         write(loop,"(I4)")iloop
         do i=1,Ns
            call splot("Gloc_site."//trim(adjustl(trim(loop)))//".ipt",wr,fg(i,-L:L),append=TT)
            call splot("Sigma_site."//trim(adjustl(trim(loop)))//".ipt",wr,sigma(i,-L:L),append=TT)
         enddo

         if(success)then
            do i=1,Ns
               nii(i) = 2.d0*sum(fermi(wr,beta)*dimag(fg(i,-L:L)))/sum(dimag(fg(i,-L:L)))
            enddo
            call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.ipt",nii)
            call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.ipt",erandom)
            call splot(trim(adjustl(trim(name_dir)))//"/LSigma.ipt",sigma(1:Ns,-L:L))
            call splot(trim(adjustl(trim(name_dir)))//"/LG.ipt",fg(1:Ns,-L:L))
            afg(:)   =sum(fg(1:Ns,-L:L),dim=1)/dble(Ns)
            call splot(trim(adjustl(trim(name_dir)))//"/DOS.disorder.ipt",wr,-dimag(afg(-L:L))/pi,append=TT)
         endif
      endif
      return
    end subroutine print_out



  end program hmmpt






