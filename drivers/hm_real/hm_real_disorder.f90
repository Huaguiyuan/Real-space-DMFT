!########################################################
!     PURPOSE  : 
!     AUTHORS  : Adriano Amaricci
!########################################################
MODULE COMMON_BROYDN
  USE BROYDEN
  implicit none
  integer                :: siteId
  real(8)                :: xmu0,n,n0
  complex(8),allocatable :: fg(:,:),sigma(:,:)
  complex(8),allocatable :: fg0(:),gamma(:)
END MODULE COMMON_BROYDN



function funcv(x)
  USE RDMFT_VARS_GLOBAL
  USE COMMON_BROYDN
  USE TOOLS
  implicit none
  real(8),dimension(:),intent(in)  ::  x
  real(8),dimension(size(x))       ::  funcv
  xmu0=x(1)
  fg0 = one/(one/gamma +xmu-xmu0-U*(n-0.5d0))
  n0=sum(fermi(wr,beta)*aimag(fg0))/sum(aimag(fg0))
  funcv(1)=n-n0
  write(101+mpiID,"(3(f13.9))")n,n0,xmu0
end function funcv



program hmmpt
  USE RDMFT_VARS_GLOBAL
  USE COMMON_BROYDN
  implicit none
  integer    :: i,is
  real(8)    :: x(1),r
  logical    :: check,converged  
  complex(8),allocatable,dimension(:,:) :: sigma_tmp

  !GLOBAL INITIALIZATION:
  !=====================================================================
  include "init_global_disorder.f90"


  !ALLOCATE WORKING ARRAYS:
  !=====================================================================
  allocate(wr(L))
  allocate(fg(Ns,L),sigma(Ns,L),sigma_tmp(Ns,L))
  allocate(fg0(L),gamma(L))
  wr=linspace(-wmax,wmax,L,mesh=fmesh)



  !START DMFT LOOP SEQUENCE:
  !=====================================================================
  call setup_initial_sigma()
  iloop=0 ; converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !SOLVE G_II (GLOCAL) \FORALL FREQUENCY W\IN[-WMAX,WMAX]
     call get_gloc_mpi()

     !SOLVE IMPURITY MODEL, FOR ALL LATTICE SITES:
     call solve_impurity_mpi()

     converged=check_convergence(sigma(:,:),eps_error,Nsuccess,nloop,id=0)
     if(nread/=0.d0)call search_mu(converged)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_out(converged)
     call end_loop()
  enddo
  if(mpiID==0)call system("mv -vf *.err "//trim(adjustl(trim(name_dir)))//"/")
  call close_mpi()

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine setup_initial_sigma()
    logical :: check
    inquire(file="LSigma.ipt",exist=check)
    if(check)then
       if(mpiID==0)then
          call msg("Reading Sigma in input")
          call sread("LSigma.ipt",sigma(1:Ns,1:L))
       endif
       call MPI_BCAST(sigma,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    else
       sigma=u*(n-0.5d0)
    endif
  end subroutine setup_initial_sigma



  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine get_gloc_mpi() 
    complex(8) :: zeta,Gloc(Ns,Ns),gf_tmp(Ns,1:L)
    integer    :: i
    call msg("Get local GF:",id=0)
    call start_timer
    gf_tmp=zero ; fg=zero
    do i=1+mpiID,L,mpiSIZE
       zeta  = cmplx(wr(i),eps,8) + xmu
       Gloc  = -H0
       do is=1,Ns
          Gloc(is,is)=Gloc(is,is) + zeta - sigma(is,i) + erandom(is)
       enddo
       call mat_inversion_sym(Gloc,Ns)
       do is=1,Ns
          gf_tmp(is,i) = Gloc(is,is)
       enddo
       call eta(i,L,999)
    enddo
    call stop_timer
    call MPI_REDUCE(gf_tmp,fg,Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(fg,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_gloc_mpi



  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine solve_impurity_mpi()
    integer    :: is
    logical    :: disorder
    disorder=.false. ; if(Wdis/=0)disorder=.true.
    call msg("Solve impurity:")
    disorder=.true.
    if(disorder)then
       call start_timer
       sigma_tmp=zero
       do is=1+mpiID,Ns,mpiSIZE
          call solve_per_site(is)
          call eta(is,Ns,998)
       enddo
       call stop_timer
       call MPI_REDUCE(sigma_tmp,sigma,Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(sigma,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    else
       call solve_per_site(is=1)
       forall(is=1:Ns)sigma(is,:)=sigma_tmp(1,:)
    endif
  end subroutine solve_impurity_mpi



  !******************************************************************
  !******************************************************************






  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine solve_per_site(is)
    complex(8),dimension(:,:),allocatable,save :: sold
    integer :: is
    character(len=4) :: loop
    siteId=is
    if(.not.allocated(sold))allocate(sold(Ns,1:L))
    sold(is,:)=sigma(is,:)
    n=sum(fermi(wr,beta)*aimag(fg(is,:)))/sum(aimag(fg(is,:)))
    gamma = one/(one/fg(is,:) + sigma(is,:))
    x(1)=xmu
    call broydn(x,check)
    xmu0=x(1)
    !Evaluate self-energy and put it into Sigma_tmp to be mpi_reduced later on
    sigma_tmp(is,:) = solve_mpt_sopt(fg0(:),wr(:),n,n0,xmu0)
    sigma_tmp(is,:) =weigth*sigma_tmp(is,:) + (1.d0-weigth)*sold(is,:)
    !sold(is,:)  =sigma(is,:)
  end subroutine solve_per_site

  subroutine solve_per_site_debug(is)
    complex(8),dimension(:,:),allocatable,save :: sold
    integer :: is
    character(len=4) :: loop
    siteId=is
    if(.not.allocated(sold))then
       allocate(sold(Ns,1:L))
       sold=zero!sigma!zero
    endif
    if(is==1)then
       write(loop,"(I4)")iloop
       call splot("fg_prova_"//trim(adjustl(trim(loop)))//".ipt",wr,fg(is,:))
       call splot("old_fg_prova_"//trim(adjustl(trim(loop)))//".ipt",wr,sold(is,:))
    endif
    if(iloop>1)then
       fg(is,:)=weigth*fg(is,:) + (1.d0-weigth)*sold(is,:)
    endif
    sold(is,:)  =fg(is,:)!sigma(is,:)
    n=sum(fermi(wr,beta)*aimag(fg(is,:)))/sum(aimag(fg(is,:)))
    if(is==1)call splot("sigma_"//trim(adjustl(trim(loop)))//".ipt",wr,sigma(is,:))
    gamma = one/(one/fg(is,:) + sigma(is,:))
    if(is==1)call splot("gamma_"//trim(adjustl(trim(loop)))//".ipt",wr,gamma(:))
    x(1)=xmu
    call broydn(x,check)
    xmu0=x(1)!sospetto...
    if(is==1)call splot("fg0_"//trim(adjustl(trim(loop)))//".ipt",wr,fg0(:))
    sigma_tmp(is,:) = solve_mpt_sopt(fg0(:),wr(:),n,n0,xmu0)
    !sigma(is,:) =weigth*sigma(is,:) + (1.d0-weigth)*sold(is,:)
    !sold(is,:)  =gamma(:)!sigma(is,:)
  end subroutine solve_per_site_debug


  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : Print out results
  !+-------------------------------------------------------------------+
  subroutine print_out(converged)
    real(8)   :: nimp
    complex(8):: afg(1:L)
    real(8)   :: nii(Ns)
    logical   :: converged
    integer   :: i
    character(len=4) :: loop
    if(mpiID==0)then
       nimp=0.d0
       do is=1,Ns
          nii(is) = 2.d0*sum(fermi(wr,beta)*dimag(fg(is,1:L)))/sum(dimag(fg(is,1:L)))
       enddo
       nimp=sum(nii(:))/dble(Ns)
       print*,"nimp  =",nimp
       call splot(trim(adjustl(trim(name_dir)))//"/navVSiloop.ipt",iloop,nimp,append=TT)

       !TO BE REMOVED
       write(loop,"(I4)")iloop
       do i=1,Ns
          call splot("Gloc_site."//trim(adjustl(trim(loop)))//".ipt",wr,fg(i,:),append=TT)
          call splot("Sigma_site."//trim(adjustl(trim(loop)))//".ipt",wr,sigma(i,:),append=TT)
       enddo

       if(converged)then
          call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.ipt",nii)
          call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.ipt",erandom)
          call splot(trim(adjustl(trim(name_dir)))//"/LSigma.ipt",sigma(1:Ns,:))
          call splot(trim(adjustl(trim(name_dir)))//"/LG.ipt",fg(1:Ns,:))
          afg(:)   =sum(fg(1:Ns,1:L),dim=1)/dble(Ns)
          call splot(trim(adjustl(trim(name_dir)))//"/DOS.disorder.ipt",wr,-dimag(afg(:))/pi,append=TT)
       endif
    endif
    return
  end subroutine print_out



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


  function get_naverage(id) result(naverage)
    integer :: is,id
    real(8) :: naverage
    real(8) :: nii(Ns)
    if(mpiID==id)then
       naverage=0.d0
       do is=1,Ns
          nii(is) = 2.d0*sum(fermi(wr,beta)*dimag(fg(is,:)))/sum(dimag(fg(is,:)))
       enddo
       naverage = sum(nii(:))/real(Ns,8)
    endif
  end function get_naverage



  !******************************************************************
  !******************************************************************



end program hmmpt

