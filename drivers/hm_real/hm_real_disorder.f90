!########################################################
!PURPOSE  :solve the disordered (D) Hubbard model (HM) 
!using Modified Perturbation THeory (MPT), w/ DMFT
!COMMENTS :
!disorder realization depends on the parameter int idum:
!so that different realizations (statistics) are performed 
!calling this program many times while providing a *different* seed.
!the result of each calculation is stored in a different dir
!AUTHORS  : A.Amaricci
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



program hmmpt_real_disorder
  USE RDMFT_VARS_GLOBAL
  USE COMMON_BROYDN
  implicit none
  integer                               :: i,is
  real(8)                               :: x(1),r
  logical                               :: check,converged  
  complex(8),allocatable,dimension(:,:) :: sigma_tmp
  real(8),allocatable,dimension(:)      :: nii_tmp,dii_tmp


  !GLOBAL INITIALIZATION:
  !=====================================================================
  include "init_global_disorder.f90"


  !ALLOCATE WORKING ARRAYS:
  !=====================================================================
  allocate(fg(Ns,L),sigma(Ns,L))
  allocate(fg0(L),gamma(L))
  !
  allocate(sigma_tmp(Ns,L))
  allocate(nii_tmp(Ns),dii_tmp(Ns))
  !
  allocate(wr(L))
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
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)


contains

  !******************************************************************
  !******************************************************************

  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine setup_initial_sigma()
    logical :: check
    if(mpiID==0)then
       inquire(file="LSigma_realw.data",exist=check1)
       if(.not.check1)inquire(file="LSigma_realw.data.gz",exist=check)
       if(check)then
          call msg(bg_yellow("Reading Self-energy from file:"),lines=2)
          call sread("LSigma_realw.data",sigma(1:Ns,1:L),wr(1:L))
       endif
    else
       call msg(bg_yellow("Using Hartree-Fock self-energy"),lines=2)
       sigma=zero!\==u*(n-0.5d0)
    endif
    call MPI_BCAST(sigma,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
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
       Gloc  = zero-H0
       do is=1,Ns
          Gloc(is,is)=Gloc(is,is) + zeta - sigma(is,i) - erandom(is)
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
       nii_tmp  =zero
       do is=1+mpiID,Ns,mpiSIZE
          call solve_per_site(is)
          call eta(is,Ns,998)
       enddo
       call stop_timer
       call MPI_REDUCE(sigma_tmp,sigma,Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(sigma,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       !
       call MPI_REDUCE(nii_tmp,nii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(nii,Ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
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
    siteId=is
    if(.not.allocated(sold))allocate(sold(Ns,1:L))
    sold(is,:)=sigma(is,:)
    n=sum(fermi(wr,beta)*aimag(fg(is,:)))/sum(aimag(fg(is,:)))
    !
    nii_tmp(is)=2.d0*n
    !
    gamma = one/(one/fg(is,:) + sigma(is,:))
    x(1)=xmu
    call broydn(x,check)
    xmu0=x(1)
    !Evaluate self-energy and put it into Sigma_tmp to be mpi_reduced later on
    !
    sigma_tmp(is,:) = solve_mpt_sopt(fg0(:),wr(:),n,n0,xmu0)
    sigma_tmp(is,:) =weigth*sigma_tmp(is,:) + (1.d0-weigth)*sold(is,:)
    !
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
    integer                               :: i,j,is,M,row,col
    real(8)                               :: nimp
    real(8),dimension(Ns)                 :: cdwii,rii,sii,zii
    real(8)                               :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Ns)               :: data_covariance
    real(8),dimension(2)                  :: data_mean,data_sdev
    real(8),dimension(2,2)                :: covariance_nd
    logical                               :: converged
    real(8),allocatable,dimension(:)      :: wm
    complex(8),allocatable,dimension(:,:) :: fgm,sigm,fgm_tmp,sigm_tmp
    complex(8)                            :: afg(1:L),asig(1:L)

    if(mpiID==0)then
       nimp=sum(nii(:))/dble(Ns)
       print*,"nimp  =",nimp
       call splot(trim(adjustl(trim(name_dir)))//"/navVSiloop.data",iloop,nimp,append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/LSigma_realw.data",sigma(1:Ns,1:L),wr(1:L))
       call splot(trim(adjustl(trim(name_dir)))//"/LG_realw.data",fg(1:Ns,1:L),wr(1:L))
    endif

    if(converged)then
       !Plot averaged local functions
       if(mpiID==0)then
          afg(:)  = sum(fg(1:Ns,1:L),dim=1)/dble(Ns)
          asig(:) = sum(sigma(1:Ns,1:L),dim=1)/dble(Ns)
          call splot(trim(adjustl(trim(name_dir)))//"/aDOS.data",wr,-dimag(afg(:))/pi)
          call splot(trim(adjustl(trim(name_dir)))//"/aSigma_realw.data",wr,asig)
       endif


       !Get Matsubara GF at every sites, we used to get Z
       M = max(8192,8*L) ; allocate(wm(M))
       wm(:) = pi/beta*real(2*arange(1,M)-1,8)
       allocate(fgm(Ns,L),sigm(Ns,L))
       allocate(fgm_tmp(Ns,M),sigm_tmp(Ns,M))
       fgm_tmp =zero
       sigm_tmp=zero
       do is=1+mpiID,Ns,mpiSIZE
          call get_matsubara_gf_from_dos(wr,fg(is,1:L),fgm_tmp(is,1:M),beta)
          call get_matsubara_gf_from_dos(wr,sigma(is,1:L),sigm_tmp(is,1:M),beta)
       enddo
       deallocate(fg,sigma)
       call MPI_REDUCE(fgm_tmp(:,1:L),fgm,Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_REDUCE(sigm_tmp(:,1:L),sigm,Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(fgm,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       call MPI_BCAST(sigm,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       !
       deallocate(fgm_tmp,sigm_tmp)
       !
       if(mpiID==0)then
          call splot(trim(adjustl(trim(name_dir)))//"/LG_iw.data",fgm,wm)
          call splot(trim(adjustl(trim(name_dir)))//"/LSigma_iw.data",sigm,wm)
       endif


       !Plot observables: n,n_cdw,rho,sigma,zeta
       if(mpiID==0)then
          do is=1,Ns
             row=irow(is) 
             col=icol(is)
             cdwii(is) = (-1.d0)**(row+col)*nii(is)
             sii(is)   = dimag(sigm(is,1))-&
                  wm(1)*(dimag(sigm(is,2))-dimag(sigm(is,1)))/(wm(2)-wm(1))
             rii(is)   = dimag(fgm(is,1))-&
                  wm(1)*(dimag(fgm(is,2))-dimag(fgm(is,1)))/(wm(2)-wm(1))
             zii(is)   = 1.d0/( 1.d0 + abs( dimag(sigm(is,1))/wm(1) ))
          enddo
          rii=abs(rii)
          sii=abs(sii)
          zii=abs(zii)
          call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.data",nii)
          call splot(trim(adjustl(trim(name_dir)))//"/cdwVSisite.data",cdwii)
          call splot(trim(adjustl(trim(name_dir)))//"/rhoVSisite.data",rii)
          call splot(trim(adjustl(trim(name_dir)))//"/sigmaVSisite.data",sii)
          call splot(trim(adjustl(trim(name_dir)))//"/zetaVSisite.data",zii)
          call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.data",erandom)


          call get_moments(nii,mean,sdev,var,skew,kurt)
          data_mean(1)=mean ; data_sdev(1)=sdev
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.n.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(zii,mean,sdev,var,skew,kurt)
          data_mean(2)=mean ; data_sdev(2)=sdev
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.z.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(sii,mean,sdev,var,skew,kurt)
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.sigma.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(rii,mean,sdev,var,skew,kurt)
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.rho.data",mean,sdev,var,skew,kurt)

          data_covariance(1,:)=nii
          data_covariance(2,:)=zii
          covariance_nd = get_covariance(data_covariance,data_mean)
          open(10,file=trim(adjustl(trim(name_dir)))//"/covariance_n.z.data")
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)

          forall(i=1:2,j=1:2)covariance_nd(i,j) = covariance_nd(i,j)/(data_sdev(i)*data_sdev(j))
          open(10,file=trim(adjustl(trim(name_dir)))//"/correlation_n.z.data")
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)

       endif
    endif

    return
  end subroutine print_out



  !******************************************************************
  !******************************************************************


  subroutine search_mu(convergence)
    real(8)               :: naverage
    logical,intent(inout) :: convergence
    real(8)               :: ndelta1
    integer               :: nindex1    
    if(mpiID==0)then
       naverage=sum(nii(:))/real(Ns,8)
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
       else
          ndelta=ndelta1
       endif
       xmu=xmu+real(nindex,8)*ndelta
       write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",naverage," /",nread,&
            "| shift=",nindex*ndelta,"| xmu=",xmu
       write(*,"(A,f15.12)")"dn=",abs(naverage-nread)
       print*,""
       if(abs(naverage-nread)>nerror)convergence=.false.
       call splot(trim(adjustl(trim(name_dir)))//"/muVSiter.data",iloop,xmu,abs(naverage-nread),append=.true.)
    endif
    call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine search_mu

  !******************************************************************
  !******************************************************************



end program hmmpt_real_disorder


