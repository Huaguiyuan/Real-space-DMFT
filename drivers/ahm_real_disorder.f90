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
program ahm_real_disorder
  USE RDMFT_VARS_GLOBAL
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma,sigma_tmp
  real(8),allocatable,dimension(:)        :: nii_tmp,dii_tmp
  logical                                 :: converged
  real(8)                                 :: r
  integer                                 :: i,is

  !GLOBAL INITIALIZATION:
  !===============================================================
  include "init_global_disorder.f90"


  !ALLOCATE WORKING ARRAYS
  !===============================================================
  allocate(wr(L))
  wr = linspace(-wmax,wmax,L,mesh=fmesh)

  allocate(fg(2,Ns,L))
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

     ! !##ACTHUNG!!
     converged = check_convergence_scalar(dii,eps_error,Nsuccess,nloop,&
          id=0,file=reg(name_dir)//"/error.err")

     if(nread/=0.d0)call search_mu(converged)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_sc_out(converged)
     call end_loop()
  enddo

  deallocate(fg,sigma,sigma_tmp,nii_tmp,dii_tmp)

  if(mpiID==0) then 
     open(10,file="used.inputRDMFT.in")
     write(10,nml=disorder)
     close(10)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)


contains

  !******************************************************************
  !******************************************************************


  subroutine setup_initial_sc_sigma()    
    logical :: check1,check2,check
    if(mpiID==0)then
       inquire(file="LSigma_realw.data",exist=check1)
       if(.not.check1)inquire(file="LSigma_realw.data.gz",exist=check1)
       inquire(file="LSelf_realw.data",exist=check2)
       if(.not.check2)inquire(file="LSelf_realw.data.gz",exist=check2)
       check=check1.AND.check2
       if(check)then
          call msg("Reading Self-energy from file:",lines=2)
          ! call sread("LSigma_realw.data",wr(1:L),sigma(1,1:Ns,1:L))
          ! call sread("LSelf_realw.data",wr(1:L),sigma(2,1:Ns,1:L))
          call sread("LSigma_realw.data",sigma(1,1:Ns,1:L),wr(1:L))
          call sread("LSelf_realw.data",sigma(2,1:Ns,1:L),wr(1:L))
       else
          call msg("Using Hartree-Fock-Bogoliubov self-energy",lines=2)
          sigma(1,:,:)=zero ; sigma(2,:,:)=-deltasc
       endif
    endif
    call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    !call MPI_BCAST(wr,L,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine setup_initial_sc_sigma


  !******************************************************************
  !******************************************************************


  subroutine get_sc_gloc_mpi()
    complex(8) :: Gloc(2*Ns,2*Ns),gf_tmp(2,Ns,L),zeta1,zeta2
    integer    :: i,is
    call msg("Get local GF:",id=0)
    call start_timer
    fg=zero ; gf_tmp=zero
    do i=1+mpiID,L,mpiSIZE
       Gloc=zero
       Gloc(1:Ns,1:Ns)          = -H0
       Gloc(Ns+1:2*Ns,Ns+1:2*Ns)=  H0
       do is=1,Ns
          zeta1=        cmplx(wr(i),eps,8)     + xmu - sigma(1,is,i)       - erandom(is)
          zeta2=-conjg( cmplx(wr(L+1-i),eps,8) + xmu - sigma(1,is,L+1-i) ) + erandom(is)
          Gloc(is,is)      = zeta1
          Gloc(Ns+is,Ns+is)= zeta2
          Gloc(is,Ns+is)   = -sigma(2,is,i)
          !S_12(w)=S^*(-w), by symmetry =S(w)=S_12(w)
          !we set this block to the correct function S^*(-w)
          !anyway with respect to the next call to *symmetric* 
          !matrix inversion this position is irrelevant 
          !as the routine only consider the upper blocks triangles (11,12,22)
          Gloc(Ns+is,is)   = -conjg(sigma(2,is,L+1-i))
       enddo
       !the call to *symmetry* routine enforces the symmetry condition
       !S^*(-w)=S(w)
       !this condition can be numerically broken using the generic 
       !inversion routine generating small (though sizeable) errors. 
       call matrix_inverse_sym(Gloc)
       forall(is=1:Ns)
          gf_tmp(1,is,i) = Gloc(is,is)
          gf_tmp(2,is,i) = Gloc(is,Ns+is)
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
    call msg("Solve impurity:")
    call start_timer
    sigma_tmp=zero
    nii_tmp  =zero
    dii_tmp  =zero
    do is=1+mpiID,Ns,mpiSIZE
       call solve_per_site(is)
       call eta(is,Ns,file="Impurity.eta")
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
  end subroutine solve_sc_impurity_mpi



  !******************************************************************
  !******************************************************************



  subroutine solve_per_site(is)
    integer                                      :: is
    complex(8)                                   :: det
    complex(8),dimension(:,:,:),allocatable,save :: sold
    complex(8),dimension(2,1:L)                  :: calG,fg0
    real(8)                                      :: n,n0,delta,delta0
    if(.not.allocated(sold))allocate(sold(2,Ns,1:L))
    !sold(:,is,:)=sigma(:,is,:)
    !
    n    = -sum(dimag(fg(1,is,:))*fermi(wr,beta))*fmesh/pi ! densita' per spin
    delta= -u*sum(dimag(fg(2,is,:))*fermi(wr,beta))*fmesh/pi
    !
    nii_tmp(is)=2.d0*n; dii_tmp(is)=delta
    !
    fg0=zero ; calG=zero 
    do i=1,L
       det     = fg(1,is,i)*conjg(fg(1,is,L+1-i)) + conjg(fg(2,is,L+1-i))*fg(2,is,i)
       fg0(1,i)= conjg(fg(1,is,L+1-i))/det  + sigma(1,is,i)      + u*(n-0.5d0)
       fg0(2,i)= fg(2,is,i)/det       + conjg(sigma(2,is,L+1-i)) + delta
    end do
    do i=1,L
       det     =  fg0(1,i)*conjg(fg0(1,L+1-i)) + conjg(fg0(2,L+1-i))*fg0(2,i)
       calG(1,i)=  conjg(fg0(1,L+1-i))/det
       calG(2,i)=  conjg(fg0(2,L+1-i))/det
    end do
    !
    if(iloop>1)calG(:,:) =  weight*calG(:,:) + (1.d0-weight)*sold(:,is,:)
    sold(:,is,:)  =  calG(:,:)
    !
    n0    = -sum(dimag(calG(1,:))*fermi(wr,beta))*fmesh/pi
    delta0= -u*sum(dimag(calG(2,:))*fermi(wr,beta))*fmesh/pi
    write(750,"(2I4,4(f16.12))",advance="yes")mpiID,is,n,n0,delta,delta0
    !
    sigma_tmp(:,is,:) =  solve_mpt_sc_sopt(calG,wr,n,n0,delta,delta0,L)
    !sigma_tmp(:,is,:) =  weight*sigma_tmp(:,is,:) + (1.d0-weight)*sold(:,is,:)
    !
  end subroutine solve_per_site



  !******************************************************************
  !******************************************************************



  subroutine print_sc_out(converged)
    integer                        :: i,j,is,row,col
    real(8)                        :: nimp,delta,ccdw
    real(8),dimension(Ns)          :: cdwii
    real(8),dimension(Nside,Nside) :: dij,nij,cij
    real(8),dimension(Nside)       :: grid_x,grid_y
    real(8)                        :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Ns)        :: data_covariance
    real(8),dimension(2,2)         :: covariance_nd
    real(8),dimension(2)           :: data_mean,data_sdev
    logical                        :: converged
    real(8),dimension(2,0:L)       :: fgt
    complex(8),dimension(2,L)      :: afg,asigma


    if(mpiID==0)then

       nimp = sum(nii)/dble(Ns)
       delta= sum(dii)/dble(Ns)
       ccdw = 0.d0
       do is=1,Ns
          row=irow(is)
          col=icol(is)
          ccdw = ccdw + (-1.d0)**(row+col)*(nii(is)-1.d0)
       enddo
       print*,"nimp  =",nimp
       print*,"delta =",delta
       print*,"ccdw  =",ccdw

       call splot(reg(name_dir)//"/nVSiloop.data",iloop,nimp,append=.true.)
       call splot(reg(name_dir)//"/deltaVSiloop.data",iloop,delta,append=.true.)
       call splot(reg(name_dir)//"/ccdwVSiloop.data",iloop,ccdw,append=.true.)
       call splot(reg(name_dir)//"/nVSisite.data",nii)
       call splot(reg(name_dir)//"/deltaVSisite.data",dii)


       call splot(reg(name_dir)//"/LSigma_realw.data",wr,sigma(1,1:Ns,1:L))
       call splot(reg(name_dir)//"/LSelf_realw.data",wr,sigma(2,1:Ns,1:L))
       call splot(reg(name_dir)//"/LG_realw.data",wr,fg(1,1:Ns,1:L))
       call splot(reg(name_dir)//"/LF_realw.data",wr,fg(2,1:Ns,1:L))

       !WHEN CONVERGED IS ACHIEVED PLOT ADDITIONAL INFORMATION:
       if(converged)then

          !Plot observables: n,delta,n_cdw,rho,sigma,zeta
          do is=1,Ns
             row=irow(is)
             col=icol(is)
             cdwii(is) = (-1.d0)**(row+col)*(nii(is)-1.d0)
          enddo
          do row=1,Nside
             grid_x(row)=row
             grid_y(row)=row
             do col=1,Nside
                i            = ij2site(row,col)
                nij(row,col) = nii(i)
                dij(row,col) = dii(i)
             enddo
          enddo


          call splot(reg(name_dir)//"/cdwVSisite.data",cdwii)
          call splot(reg(name_dir)//"/erandomVSisite.data",erandom)
          call splot3d(reg(name_dir)//"/3d_nVSij.ipt",grid_x,grid_y,nij)
          call splot3d(reg(name_dir)//"/3d_deltaVSij.ipt",grid_x,grid_y,dij)


          !Plot averaged local functions
          afg    = sum(fg,dim=2)/real(Ns,8)
          asigma = sum(sigma,dim=2)/real(Ns,8)

          call splot(reg(name_dir)//"/DOS.disorder.data",wr,-dimag(afg(1,:))/pi)
          call splot(reg(name_dir)//"/aSigma_realw.data",wr,asigma(1,:))
          call splot(reg(name_dir)//"/aSelf_realw.data",wr,asigma(2,:))
          call splot(reg(name_dir)//"/aG_realw.data",wr,afg(1,:))
          call splot(reg(name_dir)//"/aF_realw.data",wr,afg(2,:))


          call get_moments(nii,mean,sdev,var,skew,kurt)
          data_mean(1)=mean
          data_sdev(1)=sdev
          call splot(reg(name_dir)//"/statistics.n.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(dii,mean,sdev,var,skew,kurt)
          data_mean(2)=mean
          data_sdev(2)=sdev
          call splot(reg(name_dir)//"/statistics.delta.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(cdwii,mean,sdev,var,skew,kurt)
          call splot(reg(name_dir)//"/statistics.cdwn.data",mean,sdev,var,skew,kurt)
          !
          data_covariance(1,:)=nii
          data_covariance(2,:)=dii
          covariance_nd = get_covariance(data_covariance,data_mean)
          open(10,file=reg(name_dir)//"/covariance_n.delta.data")
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)

          forall(i=1:2,j=1:2)covariance_nd(i,j) = covariance_nd(i,j)/(data_sdev(i)*data_sdev(j))
          open(10,file=reg(name_dir)//"/correlation_n.delta.data")
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)
       end if

    end if
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
       if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
          ndelta=real(ndelta1/2.d0,8) !decreasing the step
       else
          ndelta=ndelta1
       endif
       xmu=xmu+real(nindex,8)*ndelta
       write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",naverage,"/",nread,"| shift=",nindex*ndelta,"| mu=",xmu
       write(*,"(A,f15.12,A,f15.12)")"Density Error:",abs(naverage-nread),'/',nerror
       print*,""
       if(abs(naverage-nread)>nerror)convergence=.false.
       call splot(reg(name_dir)//"/muVSiter.data",iloop,xmu,abs(naverage-nread),append=.true.)
    endif
    call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine search_mu


end program ahm_real_disorder
