!########################################################
!PURPOSE  :solve the disordered (D) Hubbard
! model (HM) using  DMFT-ED
! disorder realization depends on the parameter int idum:
! so that different realizations (statistics) are performed 
! calling this program many times providing a *different* seed 
! +IDUM. 
!########################################################
program ed_hm_disorder
  USE RDMFT
  USE ERROR
  USE TOOLS
  USE IOTOOLS
  USE ARRAYS
  USE STATISTICS
  USE DMFT_ED
  implicit none
  complex(8),allocatable,dimension(:,:) :: Smats,Sreal !self_energies
  complex(8),allocatable,dimension(:,:) :: Gmats,Greal !local green's functions
  complex(8),allocatable,dimension(:,:) :: Delta,errDelta
  real(8),allocatable,dimension(:)        :: erandom
  real(8),allocatable,dimension(:,:,:)    :: bath,bath_old,errBath
  logical                                 :: converged
  real(8)                                 :: wmixing
  integer                                 :: i,is,iloop
  integer                                 :: Nb(2)
  real(8),dimension(:),allocatable        :: Verror

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
  call parse_input_variable(wmixing,"WMIXING","inputRDMFT.in",default=0.5d0)
  call rdmft_read_input("inputRDMFT.in")
  store_size=1024
  Nlat = Nside**2
  !+-------------------------------------------------------------------+!


  !+- build lattice hamiltonian -+!
  call get_lattice_hamiltonian(Nside,Nside)
  !+-----------------------------+!


  !+- allocate matsubara and real frequencies -+!
  allocate(wm(Lmats),wr(Lreal))
  wini  =wini-Wdis
  wfin  =wfin+Wdis
  wr = linspace(wini,wfin,Lreal)
  wm(:)  = pi/beta*real(2*arange(1,Lmats)-1,8)
  !+----------------------------------+!

  !random energies 
  allocate(erandom(Nlat))
  call random_seed(put=[idum])
  call random_number(erandom)
  erandom=(2.d0*erandom-1.d0)*Wdis/2.d0
  !erandom(Nlat) = -sum(erandom(1:Nlat-1))
  call store_data("erandomVSisite.data",erandom,(/(dble(i),i=1,Nlat)/))
  call splot("erandom.data",erandom,(/(1,i=1,size(erandom))/))

  !+- allocate a bath for each impurity -+!
  Nb=get_bath_size()
  allocate(nii(Nlat))
  allocate(dii(Nlat))
  allocate(bath(Nlat,Nb(1),Nb(2)))
  allocate(bath_old(Nlat,Nb(1),Nb(2)))
  allocate(Smats(Nlat,Lmats))
  allocate(Sreal(Nlat,Lreal))
  allocate(Gmats(Nlat,Lmats))
  allocate(Greal(Nlat,Lreal))
  allocate(Delta(Nlat,Lmats))

  !+- initialize baths -+!
  call init_lattice_baths(bath)
  bath_old=bath

  !+- DMFT LOOP -+!
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(mpiID==0)call start_loop(iloop,nloop,"DMFT-loop",unit=LOGfile)
     call ed_solve_impurity(bath,erandom,Delta,Gmats,Greal,Smats,Sreal)
     if(rdmft_phsym)then
        do i=1,Nlat
           call ph_symmetrize_bath(bath_(i,:,:))
        enddo
     endif
     bath=wmixing*bath + (1.d0-wmixing)*bath_old ; bath_old = bath
     if(mpiID==0)converged = check_convergence(Delta(:,:),dmft_error,Nsuccess,nloop,index=2,total=2,id=0,file="DELTAerror.err",reset=.false.)
     if(mpiID==0)converged = check_convergence_local(dii,dmft_error,Nsuccess,nloop,index=1,total=2,id=0,file="error.err")
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_out(converged)
     if(mpiID==0)call end_loop()
  enddo
  !+-------------------------------------+!

  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)


contains


  subroutine print_out(converged)
    integer                        :: i,j,is,row,col
    real(8)                        :: nimp,phi,ccdw,docc
    real(8),dimension(Nlat)        :: cdwii,rii,sii,zii
    real(8),dimension(Nside,Nside) :: dij,nij,cij
    real(8),dimension(Nside)       :: grid_x,grid_y
    real(8)                        :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Nlat)      :: data_covariance
    real(8),dimension(2,2)         :: covariance_nd
    real(8),dimension(2)           :: data_mean,data_sdev
    logical                        :: converged
    complex(8),dimension(Lmats)    :: aGmats,aSmats
    complex(8),dimension(Lreal)    :: aGreal,aSreal
    character(len=50)              :: suffix,cloop,cfoo


    if(mpiID==0)then
       suffix=".data"
       write(cfoo,'(I4.4)')iloop
       cloop = "_loop"//trim(cfoo)//".data"

       !Get CDW "order parameter"
       do is=1,Nlat
          row=irow(is)
          col=icol(is)
          cdwii(is) = (-1.d0)**(row+col)*(nii(is)-1.d0)
       enddo

       nimp = sum(nii)/dble(Nlat)
       docc = sum(dii)/dble(Nlat)
       ccdw = sum(cdwii)/dble(Nlat)
       print*,"<nimp>  =",nimp
       print*,"<docc>  =",docc
       print*,"<ccdw>  =",ccdw


       call splot("nVSiloop.data",iloop,nimp,append=.true.)
       call splot("doccVSiloop.data",iloop,docc,append=.true.)
       call splot("ccdwVSiloop.data",iloop,ccdw,append=.true.)
       call store_data("nVSisite"//trim(suffix),nii,(/(dble(i),i=1,Nlat)/))
       call store_data("doccVSisite"//trim(suffix),dii,(/(dble(i),i=1,Nlat)/))
       call store_data("cdwVSisite"//trim(suffix),cdwii,(/(dble(i),i=1,Nlat)/))


       !<DEBUG: to be removed or moved under converged section below
       ! call store_data("nVSisite"//trim(cloop),nii)
       ! call store_data("phiVSisite"//trim(cloop),pii)
       ! call store_data("LG_iw"//trim(suffix),wm(1:Lmats),Gmats(1,1:Nlat,1:Lmats))
       ! call store_data("LG_realw"//trim(suffix),wr(1:Lreal),Greal(1,1:Nlat,1:Lreal))
       !>DEBUG


       !WHEN CONVERGED IS ACHIEVED PLOT ADDITIONAL INFORMATION:
       if(converged)then
          call store_data("LG_iw"//trim(suffix),Gmats(1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LG_realw"//trim(suffix),Greal(1:Nlat,1:Lreal),wr(1:Lreal))
          call store_data("LSigma_iw"//trim(suffix),Smats(1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LSigma_realw"//trim(suffix),Sreal(1:Nlat,1:Lreal),wr(1:Lreal))
          ! call store_data("LDelta_iw"//trim(suffix),wm(1:Lmats),Delta(1,1:Nlat,1:Lmats))

          !Plot observables: n,delta,n_cdw,rho,sigma,zeta
          do is=1,Nlat
             row=irow(is)
             col=icol(is)
             sii(is)   = dimag(Smats(is,1))-&
                  wm(1)*(dimag(Smats(is,2))-dimag(Smats(is,1)))/(wm(2)-wm(1))
             rii(is)   = dimag(Gmats(is,1))-&
                  wm(1)*(dimag(Gmats(is,2))-dimag(Gmats(is,1)))/(wm(2)-wm(1))
             zii(is)   = 1.d0/( 1.d0 + abs( dimag(Smats(is,1))/wm(1) ))
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

          call store_data("rhoVSisite"//trim(suffix),rii,(/(dble(i),i=1,Nlat)/))
          call store_data("sigmaVSisite"//trim(suffix),sii,(/(dble(i),i=1,Nlat)/))
          call store_data("zetaVSisite"//trim(suffix),zii,(/(dble(i),i=1,Nlat)/))
          call splot3d("3d_nVSij"//trim(suffix),grid_x,grid_y,nij)
          call splot3d("3d_doccVSij"//trim(suffix),grid_x,grid_y,dij)

          !Plot averaged local functions
          aGmats    = sum(Gmats,dim=1)/real(Nlat,8)
          aSmats = sum(Smats,dim=1)/real(Nlat,8)

          aGreal    = sum(Greal,dim=1)/real(Nlat,8)
          aSreal = sum(Sreal,dim=1)/real(Nlat,8)


          call splot("aSigma_iw"//trim(suffix),wm,aSmats(:))
          call splot("aSigma_realw"//trim(suffix),wr,aSreal(:))

          call splot("aG_iw"//trim(suffix),wm,aGmats(:))
          call splot("aG_realw"//trim(suffix),wr,aGreal(:))


          call get_moments(nii,mean,sdev,var,skew,kurt)
          data_mean(1)=mean
          data_sdev(1)=sdev
          call splot("statistics.n"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(dii,mean,sdev,var,skew,kurt)
          data_mean(2)=mean
          data_sdev(2)=sdev
          call splot("statistics.docc"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(cdwii,mean,sdev,var,skew,kurt)
          call splot("statistics.cdwn"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(zii,mean,sdev,var,skew,kurt)
          call splot("statistics.zeta"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(sii,mean,sdev,var,skew,kurt)
          call splot("statistics.sigma"//trim(suffix),mean,sdev,var,skew,kurt)
          !
          call get_moments(rii,mean,sdev,var,skew,kurt)
          call splot("statistics.rho"//trim(suffix),mean,sdev,var,skew,kurt)

          data_covariance(1,:)=nii
          data_covariance(2,:)=dii
          covariance_nd = get_covariance(data_covariance,data_mean)
          open(10,file="covariance_n.docc"//trim(suffix))
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)

          forall(i=1:2,j=1:2)covariance_nd(i,j) = covariance_nd(i,j)/(data_sdev(i)*data_sdev(j))
          open(10,file="correlation_n.docc"//trim(suffix))
          do i=1,2
             write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
          enddo
          close(10)
       end if

    end if
  end subroutine print_out



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


end program ed_hm_disorder
