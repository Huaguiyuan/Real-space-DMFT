!########################################################
!PURPOSE  :solve the attractive (A) disordered (D) Hubbard
! model (HM) using Modified Perturbation THeory (MPT), w/ DMFT
! disorder realization depends on the parameter int idum:
! so that different realizations (statistics) are performed 
! calling this program many times providing a *different* seed 
! +IDUM. The result of each calculation is stored different dir
! indexed by the seed itself.
!AUTHORS  : A.Amaricci, A.Priviter (CNR-IOM)
!########################################################
program ahm_matsubara_disorder
  USE RDMFT_VARS_GLOBAL
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma,sigma_tmp
  complex(8),allocatable,dimension(:,:)   :: fg0
  real(8),allocatable,dimension(:)        :: nii_tmp,dii_tmp
  logical                                 :: converged,converged1,converged2
  real(8)                                 :: r
  integer                                 :: i,is,Lerr

  !GLOBAL INITIALIZATION:
  !===============================================================
  include "init_global_disorder.f90"


  !ALLOCATE WORKING ARRAYS
  !===============================================================
  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)
  Lerr=min(1024,L)

  write(*,"(A,I9,A)")"Using ",L," frequencies"
  allocate(fg(2,Ns,L))
  allocate(fg0(2,L))
  allocate(sigma(2,Ns,L))
  !
  allocate(sigma_tmp(2,Ns,L))
  allocate(nii_tmp(Ns),dii_tmp(Ns))

  !START DMFT LOOP SEQUENCE:
  !==============================================================
  call setup_initial_sc_sigma()
  iloop=0 ; converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !SOLVE G_II (GLOCAL)
     call get_sc_gloc_mpi()      

     !SOLVE IMPURITY MODEL, \FORALL LATTICE SITES:
     call solve_sc_impurity_mpi()

     !converged=check_convergence_local(sigma(1,:,1:Lerr)+sigma(2,:,1:Lerr),eps_error,Nsuccess,nloop,id=0)
     converged = check_convergence_scalar(dii,eps_error,Nsuccess,nloop,id=0)
     ! !##ACTHUNG!!
     ! converged1=check_convergence_local(fg(1,:,1:Lerr),eps_error,Nsuccess,nloop,id=0,index=1,total=2)
     ! converged2=check_convergence_local(fg(2,:,1:Lerr),eps_error,Nsuccess,nloop,id=0,index=2,total=2)
     ! converged=converged1*converged2

     if(nread/=0.d0)call search_mu(converged)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_sc_out(converged)
     call end_loop()
  enddo
  if(mpiID==0)call system("mv -vf *.err "//trim(adjustl(trim(name_dir)))//"/")
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)



contains

  !******************************************************************
  !******************************************************************


  subroutine setup_initial_sc_sigma()    
    logical :: check1,check2,check
    if(mpiID==0)then
       inquire(file="LSigma_iw.data",exist=check1)
       if(.not.check1)inquire(file="LSigma_iw.data.gz",exist=check1)
       inquire(file="LSelf_iw.data",exist=check2)
       if(.not.check2)inquire(file="LSelf_iw.data.gz",exist=check2)
       check=check1.AND.check2
       if(check)then
          call msg("Reading Self-energy from file:",lines=2)
          call sread("LSigma_iw.data",sigma(1,1:Ns,1:L),wm)
          call sread("LSelf_iw.data",sigma(2,1:Ns,1:L),wm)
       else
          call msg("Using Hartree-Fock-Bogoliubov self-energy",lines=2)
          sigma(1,:,:)=zero ; sigma(2,:,:)=-deltasc
       endif
    endif
    call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  end subroutine setup_initial_sc_sigma


  !******************************************************************
  !******************************************************************


  subroutine get_sc_gloc_mpi()
    complex(8) :: Gloc(2*Ns,2*Ns),gf_tmp(2,Ns,L)
    integer    :: is
    call msg("Get local GF:",id=0)
    call start_timer
    fg=zero ; gf_tmp=zero
    do i=1+mpiID,L,mpiSIZE
       Gloc=zero
       Gloc(1:Ns,1:Ns)          = -H0
       Gloc(Ns+1:2*Ns,Ns+1:2*Ns)=  H0
       do is=1,Ns
          Gloc(is,is)      =  xi*wm(i)-sigma(1,is,i)       -erandom(is)+xmu
          Gloc(Ns+is,Ns+is)=  xi*wm(i)+conjg(sigma(1,is,i))+erandom(is)-xmu !-conjg(Gloc(is,is))
          Gloc(is,Ns+is)   = -sigma(2,is,i)
          Gloc(Ns+is,is)   = -sigma(2,is,i)
       enddo
       call matrix_inverse_sym(Gloc)
       forall(is=1:Ns)
          gf_tmp(1,is,i) = Gloc(is,is)
          !##ACTHUNG!!
          gf_tmp(2,is,i) = real(Gloc(is,Ns+is),8)
       end forall
       call eta(i,L,file="Glocal.eta")
    enddo
    call stop_timer
    call MPI_REDUCE(gf_tmp,fg,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(fg,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_sc_gloc_mpi


  !******************************************************************
  !******************************************************************



  subroutine solve_sc_impurity_mpi()
    integer :: is,i
    logical :: disorder
    disorder=.false. ; if(Wdis/=0.d0)disorder=.true.
    call msg("Solve impurity:")
    disorder=.true.
    if(disorder)then
       call start_timer
       sigma_tmp=zero
       nii_tmp  =zero
       dii_tmp  =zero
       do is=1+mpiID,Ns,mpiSIZE
          call solve_per_site(is)
          call eta(is,Ns)
       enddo
       call stop_timer
       call MPI_REDUCE(sigma_tmp,sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
       !
       call MPI_REDUCE(nii_tmp,nii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(nii,Ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
       !
       call MPI_REDUCE(dii_tmp,dii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(dii,Ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    else
       call solve_per_site(is=1)
       forall(is=1:Ns)sigma(:,is,:)=sigma_tmp(:,1,:)
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
    real(8)                                      :: n,n0,delta,delta0
    if(.not.allocated(sold))allocate(sold(2,Ns,L))
    !sold(:,is,:) = sigma(:,is,:)
    !
    call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
    call fftgf_iw2tau(fg(2,is,:),fgt(2,0:L),beta,notail=.true.)
    n = -fgt(1,L) ; delta = -u*fgt(2,L)
    !
    nii_tmp(is)=2.d0*n; dii_tmp(is)=delta
    !
    fg0=zero ; calG=zero
    det       = abs(fg(1,is,:))**2    + fg(2,is,:)**2
    fg0(1,:)  = conjg(fg(1,is,:))/det + sigma(1,is,:) + U*(n-0.5d0)
    fg0(2,:)  = fg(2,is,:)/det        + sigma(2,is,:) + delta
    det       = abs(fg0(1,:))**2      + fg0(2,:)**2
    calG(1,:) = conjg(fg0(1,:))/det
    calG(2,:) = fg0(2,:)/det
    !
    if(iloop>1)calG(:,:) =  weight*calG(:,:) + (1.d0-weight)*sold(:,is,:)
    sold(:,is,:)  =  calG(:,:)
    !
    call fftgf_iw2tau(calG(1,:),fg0t(1,:),beta)
    call fftgf_iw2tau(calG(2,:),fg0t(2,:),beta,notail=.true.)
    n0=-fg0t(1,L) ; delta0= -u*fg0t(2,L)

    write(750,"(2I4,4(f16.12))",advance="yes")mpiID,is,n,n0,delta,delta0
    !
    sigma_tmp(:,is,:) =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)
    !sigma_tmp(:,is,:) =  weight*sigma_tmp(:,is,:) + (1.d0-weight)*sold(:,is,:)
    !
    !!##ACTHUNG!!!
    sigma_tmp(2,is,:) = real(sigma_tmp(2,is,:),8)
  end subroutine solve_per_site



  !******************************************************************
  !******************************************************************



  subroutine print_sc_out(converged)
    integer                   :: i,j,is,row,col
    real(8)                   :: nimp,delta
    real(8),dimension(Ns)     :: cdwii,rii,sii,zii
    real(8)                   :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Ns)   :: data_covariance
    real(8),dimension(2,2)    :: covariance_nd
    real(8),dimension(2)      :: data_mean,data_sdev
    logical                   :: converged
    real(8),dimension(2,0:L)  :: fgt
    complex(8),dimension(2,L) :: afg,asigma
    character(len=4) :: loop


    if(mpiID==0)then
       write(loop,"(I4)")iloop

       nimp = sum(nii)/dble(Ns)
       delta= sum(dii)/dble(Ns)
       print*,"nimp  =",nimp
       print*,"delta =",delta
       call splot(trim(adjustl(trim(name_dir)))//"/nVSiloop.data",iloop,nimp,append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/deltaVSiloop.data",iloop,delta,append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.data",nii)
       call splot(trim(adjustl(trim(name_dir)))//"/deltaVSisite.data",dii)

       if(printf)then
          call splot(trim(adjustl(trim(name_dir)))//"/LSigma_iw.data",sigma(1,1:Ns,1:L),wm(1:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LSelf_iw.data",sigma(2,1:Ns,1:L),wm(1:L))
       endif

       if(converged)then
          if(.not.printf)then
             call splot(trim(adjustl(trim(name_dir)))//"/LSigma_iw.data",sigma(1,1:Ns,1:L),wm(1:L))
             call splot(trim(adjustl(trim(name_dir)))//"/LSelf_iw.data",sigma(2,1:Ns,1:L),wm(1:L))
          endif
          call splot(trim(adjustl(trim(name_dir)))//"/LG_iw.data",fg(1,1:Ns,1:L),wm(1:L))
          call splot(trim(adjustl(trim(name_dir)))//"/LF_iw.data",fg(2,1:Ns,1:L),wm(1:L))

          !Plot observables: n,delta,n_cdw,rho,sigma,zeta
          do is=1,Ns
             row=irow(is)
             col=icol(is)
             cdwii(is) = (-1.d0)**(row+col)*nii(is)
             sii(is)   = dimag(sigma(1,is,1))-&
                  wm(1)*(dimag(sigma(1,is,2))-dimag(sigma(1,is,1)))/(wm(2)-wm(1))
             rii(is)   = dimag(fg(1,is,1))-&
                  wm(1)*(dimag(fg(1,is,2))-dimag(fg(1,is,1)))/(wm(2)-wm(1))
             zii(is)   = 1.d0/( 1.d0 + abs( dimag(sigma(1,is,1))/wm(1) ))
          enddo
          rii=abs(rii)
          sii=abs(sii)
          zii=abs(zii)
          call splot(trim(adjustl(trim(name_dir)))//"/cdwVSisite.data",cdwii)
          call splot(trim(adjustl(trim(name_dir)))//"/rhoVSisite.data",rii)
          call splot(trim(adjustl(trim(name_dir)))//"/sigmaVSisite.data",sii)
          call splot(trim(adjustl(trim(name_dir)))//"/zetaVSisite.data",zii)
          call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.data",erandom)

          !Plot averaged local functions
          afg    = sum(fg,dim=2)/dble(Ns)
          asigma = sum(sigma,dim=2)/dble(Ns)
          call splot(trim(adjustl(trim(name_dir)))//"/aSigma_iw.data",wm,asigma(1,:))
          call splot(trim(adjustl(trim(name_dir)))//"/aSelf_iw.data",wm,asigma(2,:))
          call splot(trim(adjustl(trim(name_dir)))//"/aG_iw.data",wm,afg(1,:))
          call splot(trim(adjustl(trim(name_dir)))//"/aF_iw.data",wm,afg(2,:))


          call get_moments(nii,mean,sdev,var,skew,kurt)
          data_mean(1)=mean
          data_sdev(1)=sdev
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.n.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(dii,mean,sdev,var,skew,kurt)
          data_mean(2)=mean
          data_sdev(2)=sdev
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.delta.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(cdwii,mean,sdev,var,skew,kurt)
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.cdwn.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(zii,mean,sdev,var,skew,kurt)
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.zeta.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(sii,mean,sdev,var,skew,kurt)
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.sigma.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(rii,mean,sdev,var,skew,kurt)
          call splot(trim(adjustl(trim(name_dir)))//"/statistics.rho.data",mean,sdev,var,skew,kurt)

          data_covariance(1,:)=nii
          data_covariance(2,:)=dii
          covariance_nd = get_covariance(data_covariance,data_mean)
          open(10,file=trim(adjustl(trim(name_dir)))//"/covariance_n.delta.data")
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)

          forall(i=1:2,j=1:2)covariance_nd(i,j) = covariance_nd(i,j)/(data_sdev(i)*data_sdev(j))
          open(10,file=trim(adjustl(trim(name_dir)))//"/correlation_n.delta.data")
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)
       end if

    end if
    return
  end subroutine print_sc_out



  !******************************************************************
  !******************************************************************



  subroutine search_mu(convergence)
    integer, save         ::nindex
    integer               ::nindex1
    real(8)               :: naverage,ndelta1
    logical,intent(inout) :: convergence
    if(mpiID==0)then
       naverage=sum(nii)/dble(Ns)
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
          ndelta=real(ndelta1/2.d0,8) !decreasing the step
       else
          ndelta=ndelta1
       endif
       xmu=xmu+real(nindex,8)*ndelta
       write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",naverage,"/",nread,"| shift=",nindex*ndelta,"| mu=",xmu
       write(*,"(A,f15.12,A,f15.12)")"Density Error:",abs(naverage-nread),'/',nerror
       print*,""
       if(abs(naverage-nread)>nerror)convergence=.false.
       call splot(trim(adjustl(trim(name_dir)))//"/muVSiter.data",iloop,xmu,abs(naverage-nread),append=.true.)
    endif
    call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine search_mu


  !******************************************************************
  !******************************************************************

end program ahm_matsubara_disorder
