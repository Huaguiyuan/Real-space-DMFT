
!PURPOSE  :solve the disordered (D) Hubbard model (HM) 
! using Modified Perturbation THeory (MPT), w/ DMFT
!COMMENTS :
! DISORDER REALIZATION DEPENDS ON THE PARAMETER int idum:
! SO THAT DIFFERENT REALIZATIONS (STATISTICS) ARE PERFORMED 
! CALLING THIS PROGRAM MANY TIMES WHILE PROVIDING A *DIFFERENT* SEED.
! THE RESULT OF EACH CALCULATION IS STORED IN A DIFFERENT DIR
!AUTHORS  : A.Amaricci


MODULE COMMON_BROYDN
  USE BROYDEN
  implicit none
  integer                                 :: siteid,spinid
  real(8),dimension(2)                    :: n,n0,xmu0
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma
  complex(8),allocatable,dimension(:,:)   :: fg0,gamma
END MODULE COMMON_BROYDN

FUNCTION FUNCV(x)
  USE RDMFT_VARS_GLOBAL
  USE COMMON_BROYDN
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  xmu0(spinId)=x(1)
  fg0(spinId,:) = one/(one/gamma(spinId,:) +erandom(siteId)+xmu-xmu0(spinId)-U*(n(3-spinId)-0.5d0))
  !fg0(spinId,:) = one/(one/gamma(spinId,:) +xmu-xmu0-U*(n(3-spinId)-0.5d0))
  n0(spinId)=sum(fermi(wr,beta)*aimag(fg0(spinId,:)))/sum(aimag(fg0(spinId,:)))
  funcv(1)=n(spinId)-n0(spinId)
  write(101+mpiID,"(3(f13.9))")n(spinId),n0(spinId),xmu0(spinId)
END FUNCTION FUNCV

program hmmpt_spin
  USE RDMFT_VARS_GLOBAL
  USE COMMON_BROYDN
  implicit none
  integer    :: i,is
  real(8)    :: x(1),r
  logical    :: check,converged  


  !GLOBAL INITIALIZATION:
  !=====================================================================
  include "init_global_disorder.f90"


  !ALLOCATE WORKING ARRAYS:
  !=====================================================================
  allocate(wr(L))
  allocate(fg(2,Ns,L),sigma(2,Ns,L))
  allocate(fg0(2,L),gamma(2,L))
  wr=linspace(-wmax,wmax,L,mesh=fmesh)


  !START DMFT LOOP SEQUENCE:SOLVE FOR \SIGMA_II(W), G_II(W)=FG
  !=====================================================================
  call setup_initial_sigma()
  iloop=0 ; converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     call get_gloc_mpi()       !SOLVE G_II (GLOCAL) \FORALL FREQUENCY W
     call solve_impurity_mpi() !SOLVE IMPURITY MODEL, \FORALL LATTICE SITES
     converged=check_convergence(sigma(1,:,:),eps_error,Nsuccess,nloop,id=0)
     if(nread/=0.d0)call search_mu(converged)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_out(converged)
     call end_loop()
  enddo
  if(mpiID==0)call system("mv -vf *.err "//trim(adjustl(trim(name_dir)))//"/")
  call close_mpi()


contains


  !******************************************************************
  !******************************************************************


  subroutine setup_initial_sigma()
    logical :: check1,check2,check
    integer :: i
    inquire(file="LSigma.ipt",exist=check)
    if(check)then
       if(mpiID==0)then
          write(*,*)"Reading Sigma in input:"
          call sread("LSigma.ipt",sigma(1,1:Ns,1:L))
       endif
       sigma(2,:,:)=sigma(1,:,:)
       call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    else
       forall(i=1:2)sigma(i,:,:)=u*(n(i)-0.5d0)!+(3.d0-dble(i))*0.1d0
    endif
  end subroutine setup_initial_sigma



  !******************************************************************
  !******************************************************************


  subroutine search_mu(convergence)
    real(8)               :: naverage
    logical,intent(inout) :: convergence
    if(mpiID==0)then
       naverage=get_naverage(id=0)
       nindex1=nindex
       ndelta1=ndelta
       if((naverage >= nread+nerror))then
          nindex=-1
       elseif(naverage <= nread-nerror)then
          nindex=1
       else
          nindex=0
       endif
       if(nindex1+nindex==0)then !avoid loop forth and back
          ndelta=ndelta1/2.d0 !decreasing the step
          xmu=xmu+real(nindex,8)*ndelta
       else
          ndelta=ndelta1
          xmu=xmu+real(nindex,8)*ndelta
       endif
       write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",naverage," /",nread,&
            "| shift=",nindex*ndelta,"| xmu=",xmu
       write(*,"(A,f15.12)")"dn=",abs(naverage-nread)
       print*,""
       if(abs(naverage-nread)>nerror)convergence=.false.
       call splot(trim(adjustl(trim(name_dir)))//"/muVSiter.ipt",iloop,xmu,abs(naverage-nread),append=.true.)
    endif
    call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine search_mu



  !******************************************************************
  !******************************************************************


  function get_naverage(id) result(naverage2)
    integer :: is,id
    real(8) :: naverage(2),naverage2
    real(8) :: nii(2,Ns)
    if(mpiID==id)then
       naverage=0.d0
       do is=1,Ns
          nii(1,is) = 2.d0*sum(fermi(wr,beta)*dimag(fg(1,is,:)))/sum(dimag(fg(1,is,:)))
          nii(2,is) = 2.d0*sum(fermi(wr,beta)*dimag(fg(2,is,:)))/sum(dimag(fg(2,is,:)))
       enddo
       naverage(:) = sum(nii(:,:),dim=2)/real(Ns,8)
       naverage2=sum(naverage(:))/2.d0
    endif
  end function get_naverage

  !******************************************************************
  !******************************************************************



  subroutine get_gloc_mpi() 
    complex(8) :: Gloc(2*Ns,2*Ns),gf_tmp(2,Ns,L)
    integer    :: i,is,ispin
    call msg("Get local GF:",id=0)
    call start_timer
    gf_tmp=zero ; fg(ispin,:,:)=zero
    do i=1+mpiID,L,mpiSIZE
       Gloc=zero
       Gloc(1:Ns,1:Ns)          = -H0
       Gloc(Ns+1:2*Ns,Ns+1:2*Ns)= -H0
       do is=1,Ns
          Gloc(is,is) = cmplx(wr(i),eps,8) + xmu - sigma(1,is,i) - erandom(is)
          Gloc(Ns+is,Ns+is) = cmplx(wr(i),eps,8) + xmu - sigma(2,is,i) - erandom(is)
       enddo
       call mat_inversion_sym(Gloc,2*Ns)
       forall(is=1:Ns)
          gf_tmp(1,is,i) = Gloc(is,is)
          gf_tmp(2,is,i) = Gloc(Ns+is,Ns+is)
       end forall
       call eta(i,L,999)
    enddo
    call stop_timer
    call MPI_REDUCE(gf_tmp,fg(:,:,:),2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(fg(:,:,:),2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_gloc_mpi


  !******************************************************************
  !******************************************************************


  subroutine solve_impurity_mpi()
    integer    :: is
    complex(8) :: zsigma(2,Ns,L) !use size to be sure this has dimension == sigma
    call msg("Solve impurity:")
    if(Wdis/=0.d0)then
       zsigma=zero ; sigma=zero
       call start_timer
       do is=1+mpiID,Ns,mpiSIZE
          call eta(is,Ns,998)
          call solve_per_site(is)
       enddo
       call stop_timer
       call MPI_REDUCE(sigma,zsigma,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(zsigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       sigma=zsigma
    else
       call solve_per_site(is=1)
       forall(is=2:Ns)sigma(:,is,:)=sigma(:,1,:)
    endif
  end subroutine solve_impurity_mpi


  !******************************************************************
  !******************************************************************


  subroutine solve_per_site(is)
    integer                     :: is,ispin
    complex(8),allocatable,save :: Sold(:,:,:)
    siteId=is                  !this is for broyden
    if(.not.allocated(Sold))allocate(Sold(2,Ns,L))
    sold(ispin,is,:)  =  sigma(ispin,is,:)
    !This can't be part of the next loop: you need n_\sigma=0,1 before that
    forall(ispin=1:2)&
         n(ispin) = sum(fermi(wr,beta)*dimag(fg(ispin,is,:)))/sum(dimag(fg(ispin,is,:)))
    print*,n
    !broyden should be modified to workout the two spins at the same time: no loop!
    do spinId=1,2
       gamma(spinId,:) = one/(one/fg(spinId,is,:) + sigma(spinId,is,:))
       !gamma(spinId,:)= cmplx(wr,eps)+xmu-erandom(is)-sigma(ispin,is,:)-one/fg(ispin,is,:)
       xmu0(spinId)  = 0.d0     !this seems to be useless
       x(1)  = xmu
       call broydn(x,check)
       xmu0(spinId)  = x(1)
    enddo
    !
    do ispin=1,2
       sigma(ispin,is,:) =  solve_mpt_sopt(fg0(ispin,:),wr,n(3-ispin),n0(3-ispin),xmu0(ispin))
       sigma(ispin,is,:) =  weigth*sigma(ispin,is,:) + (1.d0-weigth)*sold(ispin,is,:)
    end do
  end subroutine solve_per_site


  !******************************************************************
  !******************************************************************


  subroutine print_out(converged)
    real(8)          :: nimp(2),nii(2,Ns)
    complex(8)       :: afg(2,L)
    logical          :: converged
    integer          :: i,is,ispin
    character(len=4) :: loop
    if(mpiID==0)then
       nimp=0.d0
       do ispin=1,2
          do is=1,Ns
             nii(ispin,is) = 2.d0*sum(fermi(wr,beta)*dimag(fg(ispin,is,:)))/sum(dimag(fg(ispin,is,:)))
          enddo
          nimp(ispin) = sum(nii(ispin,:))/dble(Ns)
       enddo
       print*,"nimp(1)  =",nimp(1)
       print*,"nimp(2)  =",nimp(2)
       call splot(trim(adjustl(trim(name_dir)))//"/navVSiloop.ipt",iloop,nimp(1),nimp(2),append=TT)


       !TO BE REMOVED
       write(loop,"(I4)")iloop
       do i=1,Ns
          call splot("Gloc1_site."//trim(adjustl(trim(loop)))//".ipt",wr,fg(1,i,:),append=TT)
          call splot("Sigma1_site."//trim(adjustl(trim(loop)))//".ipt",wr,sigma(1,i,:),append=TT)
          call splot("Gloc2_site."//trim(adjustl(trim(loop)))//".ipt",wr,fg(2,i,:),append=TT)
          call splot("Sigma2_site."//trim(adjustl(trim(loop)))//".ipt",wr,sigma(2,i,:),append=TT)
       enddo


       if(converged)then
          call splot(trim(adjustl(trim(name_dir)))//"/n1VSisite.ipt",nii(1,:))
          call splot(trim(adjustl(trim(name_dir)))//"/n2VSisite.ipt",nii(2,:))
          call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.ipt",erandom)
          call splot(trim(adjustl(trim(name_dir)))//"/LSigma1.ipt",sigma(1,:,:))
          call splot(trim(adjustl(trim(name_dir)))//"/LG1.ipt",fg(1,:,:))
          call splot(trim(adjustl(trim(name_dir)))//"/LSigma2.ipt",sigma(2,:,:))
          call splot(trim(adjustl(trim(name_dir)))//"/LG2.ipt",fg(2,:,:))
          afg(:,:)   =sum(fg(:,1:Ns,:),dim=2)/dble(Ns)
          call splot(trim(adjustl(trim(name_dir)))//"/DOS1.disorder.ipt",wr,-dimag(afg(1,:))/pi)
          call splot(trim(adjustl(trim(name_dir)))//"/DOS2.disorder.ipt",wr,-dimag(afg(2,:))/pi)
       endif
    endif
    return
  end subroutine print_out



end program hmmpt_spin






