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
  USE FUNCTIONS
  USE RANDOM,    only:nrand,init_random_number
  USE ERROR
  USE TOOLS
  USE IOTOOLS
  USE ARRAYS
  USE STATISTICS
  USE SQUARE_LATTICE
  USE DMFT_ED
  implicit none
  complex(8),allocatable,dimension(:,:) :: Smats,Sreal !self_energies
  complex(8),allocatable,dimension(:,:) :: Gmats,Greal !local green's functions
  complex(8),allocatable,dimension(:,:) :: Delta      
  real(8),allocatable,dimension(:)        :: erandom,Usite,elocal
  real(8),allocatable,dimension(:,:,:)    :: bath,bath_old
  logical                                 :: converged
  real(8)                                 :: Uperiod,Uamplitude,DeltaV
  real(8)                                 :: r,de
  real(8)                                 :: wmixing
  real(8),allocatable,dimension(:)        :: epsik,wt
  integer                                 :: i,is,iloop
  integer                                 :: Nb(2),Nx,Lk,ik


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
  !+- input variable: nr of k points for each layer -+!
  call parse_input_variable(Nx,"Nx","inputRDMFT.in",default=10)
  call parse_input_variable(Uperiod,"Uperiod","inputRDMFT.in",default=dble(Nlat))
  call parse_input_variable(Uamplitude,"Uamplitude","inputRDMFT.in",default=0.d0)
  call parse_input_variable(wmixing,"WMIXING","inputRDMFT.in",default=1.d0)
  call parse_input_variable(DeltaV,"BIAS","inputRDMFT.in",default=0.d0)
  !
  call rdmft_read_input("inputRDMFT.in")
  !if(rdmft_phsym.AND.rdmft_lrsym)stop "error! make up your mind: phsym=T and lrsym=T is not possible."
  store_size=1024
  Nlat = Nside
  allocate(nii(Nlat))
  allocate(dii(Nlat))
  !+-------------------------------------------------------------------+!


  !+- build lattice hamiltonian -+!
  call get_lattice_hamiltonian(Nside)
  !+-----------------------------+!

  !+- allocate matsubara and real frequencies -+!
  allocate(wm(Lmats),wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  wm(:)  = pi/beta*real(2*arange(1,Lmats)-1,8)
  !+----------------------------------+!


  Lk   = square_lattice_dimension(Nx)
  allocate(epsik(Lk),wt(Lk))
  wt   = square_lattice_structure(Lk,Nx)
  epsik= square_lattice_dispersion_array(Lk,ts)
  call get_free_dos(epsik,wt)


  ! allocate(epsik(Nx),wt(Nx))
  ! epsik=linspace(-2.d0*ts,2.d0*ts,Nx,mesh=de)
  ! do ik=1,Nx
  !    wt(ik)=dens_bethe(epsik(ik),2.d0*ts)*de
  ! end do


  !random energies tbr
  allocate(erandom(Nlat),Usite(Nlat),elocal(Nlat))
  erandom=0.d0
  elocal = 0.d0
  Usite=Uloc(1)
  do i=1,Nlat
     Usite(i) = Usite(i) + Uamplitude*dsin(pi*dble(i-1)/dble(Nlat-1)*Uperiod)
     elocal(i)= DeltaV*0.5d0 - DeltaV/dble(Nlat-1)*dble(i-1) 
     write(77,*) i,Usite(i),elocal(i)
  end do
  elocal = elocal + erandom
  Usite=0.d0

  !+- allocate a bath for each impurity -+!
  Nb=get_bath_size()
  allocate(bath(Nlat,Nb(1),Nb(2)))
  allocate(bath_old(Nlat,Nb(1),Nb(2)))
  allocate(Smats(Nlat,Lmats))
  allocate(Sreal(Nlat,Lreal))
  allocate(Gmats(Nlat,Lmats))
  allocate(Greal(Nlat,Lreal))
  allocate(Delta(Nlat,Lmats))
  !+- initialize baths -+!
  call init_lattice_baths(bath)

  !+- DMFT LOOP -+!
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(mpiID==0) call start_loop(iloop,nloop,"DMFT-loop")
     bath_old=bath
     call ed_solve_impurity(bath,elocal,epsik,wt,Delta,Gmats,Greal,Smats,Sreal)
     bath=wmixing*bath + (1.d0-wmixing)*bath_old
     if(rdmft_phsym)then
        do i=1,Nlat
           call ph_symmetrize_bath(bath(i,:,:))
        enddo
     endif
     if(rdmft_lrsym)call lr_symmetrize_bath(bath)
     if(mpiID==0) converged = check_convergence_local(dii,dmft_error,Nsuccess,nloop,id=0,file="error.err")
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_out(converged)
     if(mpiID==0) call end_loop()
  enddo
  !+-------------------------------------+!

  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)


contains


  subroutine print_out(converged)
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
    complex(8),dimension(Lmats)  :: aGmats,aSmats
    complex(8),dimension(Lreal)  :: aGreal,aSreal
    character(len=4)               :: loop


    if(mpiID==0)then
       write(loop,"(I4)")iloop
       nimp = sum(nii)/dble(Nlat)
       docc = sum(dii)/dble(Nlat)
       print*,"<nimp>  =",nimp
       print*,"<docc>  =",docc

       call splot("nVSiloop.data",iloop,nimp,append=.true.)
       call splot("doccVSiloop.data",iloop,docc,append=.true.)
       call store_data("nVSisite.data",nii,(/(dble(i),i=1,Nlat)/))
       call store_data("doccVSisite.data",dii,(/(dble(i),i=1,Nlat)/))

       !WHEN CONVERGED IS ACHIEVED PLOT ADDITIONAL INFORMATION:
       if(converged)then
          !
          call store_data("LG_iw.data",Gmats(1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LG_realw.data",Greal(1:Nlat,1:Lreal),wr(1:Lreal))
          call store_data("LSigma_iw.data",Smats(1:Nlat,1:Lmats),wm(1:Lmats))
          call store_data("LSigma_realw.data",Sreal(1:Nlat,1:Lreal),wr(1:Lreal))
          !
          do is=1,Nlat
             sii(is)   = dimag(Smats(is,1))-&
                  wm(1)*(dimag(Smats(is,2))-dimag(Smats(is,1)))/(wm(2)-wm(1))
             rii(is)   = dimag(Gmats(is,1))-&
                  wm(1)*(dimag(Gmats(is,2))-dimag(Gmats(is,1)))/(wm(2)-wm(1))
             zii(is)   = 1.d0/( 1.d0 + abs( dimag(Smats(is,1))/wm(1) ))
          enddo
          rii=abs(rii)
          sii=abs(sii)
          zii=abs(zii)

          call store_data("rhoVSisite.data",rii,(/(dble(i),i=1,Nlat)/))
          call store_data("sigmaVSisite.data",sii,(/(dble(i),i=1,Nlat)/))
          call store_data("zetaVSisite.data",zii,(/(dble(i),i=1,Nlat)/))
          call store_data("erandomVSisite.data",erandom,(/(dble(i),i=1,Nlat)/))
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







end program ed_slab
