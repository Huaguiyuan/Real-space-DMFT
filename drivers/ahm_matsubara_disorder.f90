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
  USE RDMFT
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma
  real(8),allocatable,dimension(:) :: erandom
  logical                                 :: converged
  real(8)                                 :: r
  integer                                 :: i,is

  !START MPI:
  !=====================================================================
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)

  !READ INPUT FILES:
  !=====================================================================
  call read_input("inputIPT.in")
  call rdmft_read_input("inputRDMFT.in")


  !BUILD THE LATTICE HAMILTONIAN:
  !=====================================================================
  call get_tb_hamiltonian


  !ALLOCATE WORKING ARRAYS:
  !=====================================================================
  wmax  =wmax+Wdis
  allocate(erandom(Ns))
  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)
  !
  allocate(fg(2,Ns,L))
  allocate(sigma(2,Ns,L))


  !BUILD RANDOM ENERGIES:
  !=====================================================================
  do i=1,100                     !get rid of few spurious random number in NR
     r=nrand(idum)
  enddo
  do is=1,Ns
     erandom(is)=(2.d0*nrand(idum)-1.d0)*Wdis/2.d0
  enddo


  !START DMFT LOOP SEQUENCE:
  !==============================================================
  call setup_sc_initial_sigma(sigma)
  iloop=0 ; converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !SOLVE G_II (GLOCAL)
     call get_sc_gloc_matsu_mpi(erandom,sigma,fg)      

     !SOLVE IMPURITY MODEL, \FORALL LATTICE SITES:
     call solve_sc_impurity_matsu_mpi(fg,sigma)

     converged = check_convergence_scalar(dii,eps_error,Nsuccess,nloop,&
          id=0,file="error.err")

     if(nread/=0.d0)call search_mu(converged)
     if(iloop>=nloop)converged=.true.
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_sc_out(converged)
     call end_loop()
  enddo

  deallocate(fg,sigma)

  if(mpiID==0) then 
     open(10,file="used.inputRDMFT.in")
     write(10,nml=disorder)
     close(10)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)


contains


  subroutine print_sc_out(converged)
    integer                         :: i,j,is,row,col
    real(8)                         :: nimp,delta,ccdw
    real(8),dimension(Ns)           :: cdwii,rii,sii,zii
    real(8),dimension(Nside,Nside) :: dij,nij,cij
    real(8),dimension(Nside)        :: grid_x,grid_y
    real(8)                         :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Ns)         :: data_covariance
    real(8),dimension(2,2)          :: covariance_nd
    real(8),dimension(2)            :: data_mean,data_sdev
    logical                         :: converged
    real(8),dimension(2,0:L)        :: fgt
    complex(8),dimension(2,L)       :: afg,asigma
    character(len=4)                :: loop


    if(mpiID==0)then
       write(loop,"(I4)")iloop

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

       call splot("nVSiloop.data",iloop,nimp,append=.true.)
       call splot("deltaVSiloop.data",iloop,delta,append=.true.)
       call splot("ccdwVSiloop.data",iloop,ccdw,append=.true.)
       call store_data("nVSisite.data",nii)
       call store_data("deltaVSisite.data",dii)

       !
       call splot("LSigma_iw.data",wm(1:L),sigma(1,1:Ns,1:L))
       call splot("LSelf_iw.data",wm(1:L),sigma(2,1:Ns,1:L))
       !

       !WHEN CONVERGED IS ACHIEVED PLOT ADDITIONAL INFORMATION:
       if(converged)then
          call splot("LG_iw.data",wm(1:L),fg(1,1:Ns,1:L))
          call splot("LF_iw.data",wm(1:L),fg(2,1:Ns,1:L))

          !Plot observables: n,delta,n_cdw,rho,sigma,zeta
          do is=1,Ns
             row=irow(is)
             col=icol(is)
             cdwii(is) = (-1.d0)**(row+col)*(nii(is)-1.d0)
             sii(is)   = dimag(sigma(1,is,1))-&
                  wm(1)*(dimag(sigma(1,is,2))-dimag(sigma(1,is,1)))/(wm(2)-wm(1))
             rii(is)   = dimag(fg(1,is,1))-&
                  wm(1)*(dimag(fg(1,is,2))-dimag(fg(1,is,1)))/(wm(2)-wm(1))
             zii(is)   = 1.d0/( 1.d0 + abs( dimag(sigma(1,is,1))/wm(1) ))
          enddo
          rii=abs(rii)
          sii=abs(sii)
          zii=abs(zii)
          do row=1,Nside
             grid_x(row)=row
             grid_y(row)=row
             do col=1,Nside
                i            = ij2site(row,col)
                nij(row,col) = nii(i)
                dij(row,col) = dii(i)
             enddo
          enddo


          call store_data("cdwVSisite.data",cdwii)
          call store_data("rhoVSisite.data",rii)
          call store_data("sigmaVSisite.data",sii)
          call store_data("zetaVSisite.data",zii)
          call store_data("erandomVSisite.data",erandom)
          call splot3d("3d_nVSij.ipt",grid_x,grid_y,nij)
          call splot3d("3d_deltaVSij.ipt",grid_x,grid_y,dij)


          !Plot averaged local functions
          afg    = sum(fg,dim=2)/real(Ns,8)
          asigma = sum(sigma,dim=2)/real(Ns,8)

          call splot("aSigma_iw.data",wm,asigma(1,:))
          call splot("aSelf_iw.data",wm,asigma(2,:))
          call splot("aG_iw.data",wm,afg(1,:))
          call splot("aF_iw.data",wm,afg(2,:))


          call get_moments(nii,mean,sdev,var,skew,kurt)
          data_mean(1)=mean
          data_sdev(1)=sdev
          call splot("statistics.n.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(dii,mean,sdev,var,skew,kurt)
          data_mean(2)=mean
          data_sdev(2)=sdev
          call splot("statistics.delta.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(cdwii,mean,sdev,var,skew,kurt)
          call splot("statistics.cdwn.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(zii,mean,sdev,var,skew,kurt)
          call splot("statistics.zeta.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(sii,mean,sdev,var,skew,kurt)
          call splot("statistics.sigma.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(rii,mean,sdev,var,skew,kurt)
          call splot("statistics.rho.data",mean,sdev,var,skew,kurt)

          data_covariance(1,:)=nii
          data_covariance(2,:)=dii
          covariance_nd = get_covariance(data_covariance,data_mean)
          open(10,file="covariance_n.delta.data")
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)

          forall(i=1:2,j=1:2)covariance_nd(i,j) = covariance_nd(i,j)/(data_sdev(i)*data_sdev(j))
          open(10,file="correlation_n.delta.data")
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
          ndelta=ndelta1/2.d0
       else
          ndelta=ndelta1
       endif
       xmu=xmu+real(nindex,8)*ndelta
       write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",naverage,"/",nread,"| shift=",nindex*ndelta,"| mu=",xmu
       write(*,"(A,f15.12,A,f15.12)")"Density Error:",abs(naverage-nread),'/',nerror
       print*,""
       if(abs(naverage-nread)>nerror)convergence=.false.
       call splot("muVSiter.data",iloop,xmu,abs(naverage-nread),append=.true.)
    endif
    call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine search_mu


end program ahm_matsubara_disorder
