!########################################################
!PURPOSE  :solve the attractive (A) disordered (D) Hubbard
! model (HM) using Modified Perturbation THeory (MPT), w/ DMFT
! disorder realization depends on the parameter int idum:
! so that different realizations (statistics) are performed 
! calling this program many times providing a *different* seed 
! +IDUM. The result of each calculation is stored different dir
! indexed by the seed itself.
!########################################################
program ed_slab
  USE RDMFT
  USE RANDOM,    only:nrand,init_random_number
  USE ERROR
  USE TOOLS
  USE IOTOOLS
  USE ARRAYS
  USE STATISTICS
  USE SQUARE_LATTICE
  USE DMFT_ED
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: Smats,Sreal !self_energies
  complex(8),allocatable,dimension(:,:,:) :: Gmats,Greal !local green's functions
  complex(8),allocatable,dimension(:,:,:) :: Delta      
  real(8),allocatable,dimension(:)        :: erandom,Usite
  real(8),allocatable,dimension(:,:,:)    :: bath
  logical                                 :: converged
  real(8)                                 :: Uperiod,Uamplitude
  real(8)                                 :: r
  integer                                 :: i,is,iloop
  integer                                 :: Nb(2),Nx,Lk
  

  !+-------------------------------------------------------------------+!
  ! START MPI !
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !+-------------------------------------------------------------------+!
  
  !+-------------------------------------------------------------------+!
  ! READ INPUT FILES !
  call rdmft_read_input("inputRDMFT.in")
  store_size=1024
  Nlat = Nside
  allocate(nii(Nlat))
  allocate(dii(Nlat))
  allocate(pii(Nlat))
  !+-------------------------------------------------------------------+!

  !+- input variable: nr of k points for each layer -+!
  call parse_input_variable(Nx,"Nx","inputRDMFT.in",default=10)
  call parse_input_variable(Uperiod,"Uperiod","inputRDMFT.in",default=dble(Nlat))
  call parse_input_variable(Uamplitude,"Uamplitude","inputRDMFT.in",default=0.d0)

  !+- build lattice hamiltonian -+!
  !call get_tb_hamiltonian 
  call get_slab_hamiltonian 
  !+-----------------------------+!

  !+- allocate matsubara and real frequencies -+!
  allocate(wm(Lmats),wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  wm(:)  = pi/beta*real(2*arange(1,Lmats)-1,8)
  !+----------------------------------+!


  !random energies tbr
  allocate(erandom(Nlat),Usite(Nlat))
  erandom=0.d0
  Usite=Uloc(1)
  do i=1,Nlat
     Usite(i) = Usite(i) + Uamplitude*dsin(pi*dble(i-1)/dble(Nlat-1)*Uperiod)
     write(77,*) i,Usite(i)
  end do


  !+- allocate a bath for each impurity -+!
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2),Nlat))
  allocate(Smats(2,Nlat,Lmats))
  allocate(Sreal(2,Nlat,Lreal))
  allocate(Gmats(2,Nlat,Lmats))
  allocate(Greal(2,Nlat,Lreal))
  allocate(Delta(2,Nlat,Lmats))
  !+- initialize baths -+!
  call init_lattice_baths(bath)

  
  !+- DMFT LOOP -+!
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     call ed_solve_sc_impurity_mats(bath,erandom,Usite,Nx,Delta,Gmats,Greal,Smats,Sreal)
     converged = check_convergence(pii,dmft_error,Nsuccess,nloop,id=0,file="error.err")
     call print_sc_out(converged)
     call end_loop()
  enddo
  !+-------------------------------------+!

  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)


contains


  subroutine print_sc_out(converged)
    integer                        :: i,j,is,row,col
    real(8)                        :: nimp,phi,ccdw,docc
    real(8),dimension(Nlat)        :: cdwii,rii,sii,zii
    real(8),dimension(Nside,Nside) :: dij,nij,cij,pij
    real(8),dimension(Nside)       :: grid_x,grid_y
    real(8)                        :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Nlat)      :: data_covariance
    real(8),dimension(2,2)         :: covariance_nd
    real(8),dimension(2)           :: data_mean,data_sdev
    logical                        :: converged
    complex(8),dimension(2,Lmats)  :: aGmats,aSmats
    complex(8),dimension(2,Lreal)  :: aGreal,aSreal
    character(len=4)               :: loop


    if(mpiID==0)then
       write(loop,"(I4)")iloop
       nimp = sum(nii)/dble(Nlat)
       phi  = sum(pii)/dble(Nlat)
       docc = sum(dii)/dble(Nlat)
       ccdw = 0.d0
       do is=1,Nlat
          ! row=irow(is)
          ! col=icol(is)
          ccdw = ccdw + (-1.d0)**(is)*(nii(is)-1.d0)
       enddo
       print*,"nimp  =",nimp
       print*,"phi   =",phi
       print*,"docc  =",docc
       print*,"ccdw  =",ccdw

       call splot("nVSiloop.data",iloop,nimp,append=.true.)
       call splot("phiVSiloop.data",iloop,phi,append=.true.)
       call splot("doccVSiloop.data",iloop,docc,append=.true.)
       call splot("ccdwVSiloop.data",iloop,ccdw,append=.true.)
       call store_data("nVSisite.data",nii)
       call store_data("phiVSisite.data",pii)
       call store_data("doccVSisite.data",dii)

       call splot("Delta_iw.data",wm(1:Lmats),Delta(1,1:Nlat,1:Lmats))
       call splot("Gamma_iw.data",wm(1:Lmats),Delta(2,1:Nlat,1:Lmats))

       !WHEN CONVERGED IS ACHIEVED PLOT ADDITIONAL INFORMATION:
       if(converged)then
          call splot("LG_iw.data",wm(1:Lmats),Gmats(1,1:Nlat,1:Lmats))
          call splot("LF_iw.data",wm(1:Lmats),Gmats(2,1:Nlat,1:Lmats))
          call splot("LG_realw.data",wr(1:Lreal),Greal(1,1:Nlat,1:Lreal))
          call splot("LF_realw.data",wr(1:Lreal),Greal(2,1:Nlat,1:Lreal))


          !Plot observables: n,delta,n_cdw,rho,sigma,zeta
          do is=1,Nlat
             ! row=irow(is)
             ! col=icol(is)
             cdwii(is) = (-1.d0)**(is)*(nii(is)-1.d0)
             sii(is)   = dimag(Smats(1,is,1))-&
                  wm(1)*(dimag(Smats(1,is,2))-dimag(Smats(1,is,1)))/(wm(2)-wm(1))
             rii(is)   = dimag(Gmats(1,is,1))-&
                  wm(1)*(dimag(Gmats(1,is,2))-dimag(Gmats(1,is,1)))/(wm(2)-wm(1))
             zii(is)   = 1.d0/( 1.d0 + abs( dimag(Smats(1,is,1))/wm(1) ))
          enddo
          rii=abs(rii)
          sii=abs(sii)
          zii=abs(zii)
          
          call store_data("cdwVSisite.data",cdwii)
          call store_data("rhoVSisite.data",rii)
          call store_data("sigmaVSisite.data",sii)
          call store_data("zetaVSisite.data",zii)
          call store_data("erandomVSisite.data",erandom)
          
          !Plot averaged local functions
          aGmats    = sum(Gmats,dim=2)/real(Nlat,8)
          aSmats = sum(Smats,dim=2)/real(Nlat,8)

          aGreal    = sum(Greal,dim=2)/real(Nlat,8)
          aSreal = sum(Sreal,dim=2)/real(Nlat,8)


          call splot("aSigma_iw.data",wm,aSmats(1,:))
          call splot("aSelf_iw.data",wm,aSmats(2,:))
          call splot("aSigma_realw.data",wr,aSreal(1,:))
          call splot("aSelf_realw.data",wr,aSreal(2,:))

          call splot("aG_iw.data",wm,aGmats(1,:))
          call splot("aF_iw.data",wm,aGmats(2,:))
          call splot("aG_realw.data",wr,aGreal(1,:))
          call splot("aF_realw.data",wr,aGreal(2,:))


          call get_moments(nii,mean,sdev,var,skew,kurt)
          data_mean(1)=mean
          data_sdev(1)=sdev
          call splot("statistics.n.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(dii,mean,sdev,var,skew,kurt)
          data_mean(2)=mean
          data_sdev(2)=sdev
          call splot("statistics.docc.data",mean,sdev,var,skew,kurt)
          !
          call get_moments(pii,mean,sdev,var,skew,kurt)
          data_mean(2)=mean
          data_sdev(2)=sdev
          call splot("statistics.phi.data",mean,sdev,var,skew,kurt)
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
          data_covariance(2,:)=pii
          covariance_nd = get_covariance(data_covariance,data_mean)
          open(10,file="covariance_n.phi.data")
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)

          forall(i=1:2,j=1:2)covariance_nd(i,j) = covariance_nd(i,j)/(data_sdev(i)*data_sdev(j))
          open(10,file="correlation_n.phi.data")
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
       naverage=sum(nii)/dble(Nlat)
       nindex1=nindex
       ndelta1=rdmft_ndelta
       if((naverage >= rdmft_nread+rdmft_nerror))then
          nindex=-1
       elseif(naverage <= rdmft_nread-rdmft_nerror)then
          nindex=1
       else
          nindex=0
       endif
       if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
          rdmft_ndelta=ndelta1/2.d0
       else
          rdmft_ndelta=ndelta1
       endif
       xmu=xmu+real(nindex,8)*rdmft_ndelta
       write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",naverage,"/",rdmft_nread,"| shift=",nindex*rdmft_ndelta,"| mu=",xmu
       write(*,"(A,f15.12,A,f15.12)")"Density Error:",abs(naverage-nread),'/',rdmft_nerror
       print*,""
       if(abs(naverage-rdmft_nread)>rdmft_nerror)convergence=.false.
       call splot("muVSiter.data",xmu,abs(naverage-rdmft_nread),append=.true.)
    endif
    call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine search_mu


end program ed_slab
