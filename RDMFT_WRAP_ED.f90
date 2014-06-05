!###############################################################
! PROGRAM  : RDMFT_WRAP_ED
! PURPOSE  : Contains the main function performin RDMFT calculation
! using the ED solver
! AUTHORS  : Adriano Amaricci, Giacomo Mazza
!###############################################################
module RDMFT_WRAP_ED
  USE RDMFT_INPUT_VARS
  USE RDMFT_VARS_GLOBAL
  USE RDMFT_AUX_FUNX
  USE TIMER
  !Impurity solver interface
  USE DMFT_ED
  implicit none
  private

  interface ed_solve_impurity
     module procedure ed_solve_impurity_lattice,ed_solve_impurity_slab
  end interface ed_solve_impurity

  interface ed_solve_sc_impurity
     module procedure ed_solve_sc_impurity_lattice,ed_solve_sc_impurity_slab
  end interface ed_solve_sc_impurity

  public :: ed_solve_impurity
  public :: ed_solve_sc_impurity
  public :: init_lattice_baths
  public :: lr_symmetrize_bath

contains


  !+------------------------------------------------------------+!
  ! PURPOSE  : parallel solution of impurity problem        
  !            for a purely real space problem
  !            Get Glocal
  !            Fit new bath
  !+------------------------------------------------------------+!
  subroutine ed_solve_impurity_lattice(bath_,eloc,Delta,Gmats,Greal,Smats,Sreal,Usite)
    real(8)                  :: bath_(:,:,:),bath_tmp(size(bath_,2),size(bath_,3))
    real(8)                  :: eloc(Nlat)
    real(8),optional         :: Usite(Nlat)
    complex(8),intent(inout) :: Delta(Nlat,Lmats)
    complex(8),intent(inout) :: Gmats(Nlat,Lmats)
    complex(8),intent(inout) :: Greal(Nlat,Lreal)
    complex(8),intent(inout) :: Smats(Nlat,Lmats)
    complex(8),intent(inout) :: Sreal(Nlat,Lreal)
    complex(8)               :: Smats_tmp(Nlat,Lmats)
    complex(8)               :: Sreal_tmp(Nlat,Lreal)
    complex(8)               :: Delta_tmp(Nlat,Norb,Norb,Lmats)
    complex(8)               :: calG(Lmats),cdet
    integer                  :: ilat,i
    real(8)                  :: nii_tmp(Nlat),dii_tmp(Nlat)
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    if(size(bath_,1).ne.Nlat) stop "init_lattice_bath: wrong bath size dimension 3 (Nlat)"
    do ilat=1+mpiID,Nlat,mpiSIZE
       check_dim = check_bath_dimension(bath_(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    if(mpiID==0)call start_timer
    Smats = zero ; Smats_tmp = zero
    Sreal = zero ; Sreal_tmp = zero
    Delta = zero ; Delta_tmp = zero
    nii   = 0.d0 ; nii_tmp   = 0.d0
    dii   = 0.d0 ; dii_tmp   = 0.d0
    !+- SOLVE SITE DEPENDENT IMPUTITY PROBLEM -+!
    if(mpiID/=0)LOGfile = 800+mpiID
    do ilat=1+mpiID,Nlat,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       Hloc(1,1,1,1)=eloc(ilat)
       call ed_solver(bath_(ilat,:,:))
       Smats_tmp(ilat,:) = impSmats(1,1,1,1,:)
       Sreal_tmp(ilat,:) = impSreal(1,1,1,1,:)
       nii_tmp(ilat)   = ed_dens(1)
       dii_tmp(ilat)   = ed_docc(1)
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(Smats_tmp,Smats,Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Sreal,Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    !+- GET GLOC -+!
    Gmats=zero
    Greal=zero
    call get_gloc_mats(eloc,Smats,Gmats)
    call get_gloc_real(eloc,Sreal,Greal)
    !+- GET SITE DEPENDENT HYBRIDIZATION FUNCTION AND FIT BATHS -+!
    do ilat=1+mpiID,Nlat,mpiSIZE
       Hloc(1,1,1,1)=eloc(ilat)
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       do i=1,Lmats
          if(cg_scheme=='weiss')then
             Delta_tmp(ilat,1,1,i)  =   one/(one/Gmats(ilat,i) + Smats(ilat,i))
          else
             Delta_tmp(ilat,1,1,i) = xi*wm(i) + xmu - Hloc(1,1,1,1) - Smats(ilat,i) - one/Gmats(ilat,i)
          endif
       end do
       call chi2_fitgf(Delta_tmp(ilat,:,:,:),bath_(ilat,:,:),ispin=1)
    end do
    call MPI_ALLREDUCE(Delta_tmp(1:Nlat,1,1,1:Lmats),Delta,Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_solve_impurity_lattice



  !+------------------------------------------------------------+!
  ! SC-version       
  !+------------------------------------------------------------+!
  subroutine ed_solve_sc_impurity_lattice(bath_,eloc,Delta,Gmats,Greal,Smats,Sreal,Usite)
    real(8)                  :: bath_(:,:,:),bath_tmp(size(bath_,2),size(bath_,3))
    real(8)                  :: eloc(Nlat)
    real(8),optional         :: Usite(Nlat)
    complex(8),intent(inout) :: Delta(2,Nlat,Lmats)
    complex(8),intent(inout) :: Gmats(2,Nlat,Lmats)
    complex(8),intent(inout) :: Greal(2,Nlat,Lreal)
    complex(8),intent(inout) :: Smats(2,Nlat,Lmats)
    complex(8),intent(inout) :: Sreal(2,Nlat,Lreal)
    complex(8)               :: Smats_tmp(2,Nlat,Lmats)
    complex(8)               :: Sreal_tmp(2,Nlat,Lreal)
    complex(8)               :: Delta_tmp(2,Nlat,Norb,Norb,Lmats)
    complex(8)               :: calG(2,Lmats),cdet
    integer                  :: ilat,i
    real(8)                  :: nii_tmp(Nlat),dii_tmp(Nlat),pii_tmp(Nlat)
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    if(size(bath_,1)/=Nlat) stop "init_lattice_bath: wrong bath size dimension 3 (Nlat)"
    do ilat=1+mpiID,Nlat,mpiSIZE
       check_dim = check_bath_dimension(bath_(ilat,:,:))
       if(.not.check_dim) stop "ed_solve_sc_impurity_lattice: wrong bath size dimension 1 or 2 "
    end do
    if(mpiID==0)call start_timer
    Smats = zero ; Smats_tmp = zero
    Sreal = zero ; Sreal_tmp = zero
    Delta = zero ; Delta_tmp = zero
    nii   = 0.d0 ; nii_tmp   = 0.d0
    dii   = 0.d0 ; dii_tmp   = 0.d0
    pii   = 0.d0 ; pii_tmp   = 0.d0
    !+- SOLVE SITE DEPENDENT IMPUTITY PROBLEM -+!
    if(mpiID/=0)LOGfile = 800+mpiID
    do ilat=1+mpiID,Nlat,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(present(Usite)) Uloc(1)=Usite(ilat)
       Hloc(1,1,1,1)=eloc(ilat)
       call ed_solver(bath_(ilat,:,:))
       Smats_tmp(1,ilat,:) = impSmats(1,1,1,1,:)
       Smats_tmp(2,ilat,:) = impSAmats(1,1,1,1,:)
       Sreal_tmp(1,ilat,:) = impSreal(1,1,1,1,:)
       Sreal_tmp(2,ilat,:) = impSAreal(1,1,1,1,:)
       nii_tmp(ilat)   = ed_dens(1)
       dii_tmp(ilat)   = ed_docc(1)
       pii_tmp(ilat)   = ed_phisc(1)
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(Smats_tmp,Smats,2*Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Sreal,2*Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(pii_tmp,pii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    !+- GET GLOC -+!
    Gmats=zero
    Greal=zero
    call get_sc_gloc_mats(eloc,Smats,Gmats)
    call get_sc_gloc_real(eloc,Sreal,Greal)
    !+- GET SITE DEPENDENT HYBRIDIZATION FUNCTION AND FIT BATHS -+!
    do ilat=1+mpiID,Nlat,mpiSIZE
       Hloc(1,1,1,1)=eloc(ilat)
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       do i=1,Lmats
          if(cg_scheme=='weiss')then
             cdet      = abs(Gmats(1,ilat,i))**2 + (Gmats(2,ilat,i))**2
             calG(1,i) = conjg(Gmats(1,ilat,i))/cdet + Smats(1,ilat,i)
             calG(2,i) =  Gmats(2,ilat,i)/cdet + Smats(2,ilat,i) 
             cdet                     =  abs(calG(1,i))**2 + (calG(2,i))**2
             Delta_tmp(1,ilat,1,1,i)  =  conjg(calG(1,i))/cdet
             Delta_tmp(2,ilat,1,1,i)  =  calG(2,i)/cdet
          else
             cdet       = abs(Gmats(1,ilat,i))**2 + (Gmats(2,ilat,i))**2
             Delta_tmp(1,ilat,1,1,i) = xi*wm(i) + xmu - Hloc(1,1,1,1) - Smats(1,ilat,i) - conjg(Gmats(1,ilat,i))/cdet 
             Delta_tmp(2,ilat,1,1,i) = -(Gmats(2,ilat,i)/cdet + Smats(2,ilat,i))
          endif
       end do
       call chi2_fitgf(Delta_tmp(:,ilat,:,:,:),bath_(ilat,:,:),ispin=1)
    end do
    call MPI_ALLREDUCE(Delta_tmp(:,:,1,1,:),Delta,2*Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_solve_sc_impurity_lattice








  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: solve impurity problem & get Glocal for a mixed real-k space problem 
  !          fit new bath                                                         
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_impurity_slab(bath_,eloc,epsik,wt,Delta,Gmats,Greal,Smats,Sreal,Usite)
    real(8)                          :: bath_(:,:,:),bath_tmp(size(bath_,2),size(bath_,3))
    real(8)                          :: eloc(Nlat)
    real(8),optional                 :: Usite(Nlat)
    real(8),dimension(:)             :: epsik
    real(8),dimension(size(epsik))   :: wt
    complex(8),intent(inout)         :: Delta(Nlat,Lmats)
    complex(8),intent(inout)         :: Gmats(Nlat,Lmats)
    complex(8),intent(inout)         :: Greal(Nlat,Lreal)
    complex(8),intent(inout)         :: Smats(Nlat,Lmats)
    complex(8),intent(inout)         :: Sreal(Nlat,Lreal)
    complex(8)                       :: Smats_tmp(Nlat,Lmats)
    complex(8)                       :: Sreal_tmp(Nlat,Lreal)
    complex(8)                       :: Gmats_k(Nlat,Lmats)
    complex(8)                       :: Greal_k(Nlat,Lreal)
    complex(8)                       :: Delta_tmp(Nlat,Norb,Norb,Lmats)
    complex(8)                       :: calG(Lmats),cdet
    integer                          :: ilat,i,ik
    integer                          :: Lk
    real(8)                          :: nii_tmp(Nlat),dii_tmp(Nlat)
    logical                          :: check_dim
    character(len=5)                 :: tmp_suffix
    if(size(bath_,1).ne.Nlat) stop "init_lattice_bath: wrong bath size dimension 3 (Nlat)"
    Lk=size(epsik)
    do ilat=1+mpiID,Nlat,mpiSIZE
       check_dim = check_bath_dimension(bath_(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    if(mpiID==0)call start_timer
    Smats = zero ; Smats_tmp = zero
    Sreal = zero ; Sreal_tmp = zero
    Delta = zero ; Delta_tmp = zero
    nii   = 0.d0 ; nii_tmp   = 0.d0
    dii   = 0.d0 ; dii_tmp   = 0.d0
    !+- SOLVE SITE DEPENDENT IMPURITY PROBLEM -+!
    LOGfile = 800+mpiID
    if(mpiID==0) LOGfile=6
    do ilat=1+mpiID,Nlat,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(present(Usite)) Uloc(1)=Usite(ilat)
       Hloc(1,1,1,1)=eloc(ilat)
       call ed_solver(bath_(ilat,:,:))
       Smats_tmp(ilat,:) = impSmats(1,1,1,1,:)
       Sreal_tmp(ilat,:) = impSreal(1,1,1,1,:)
       nii_tmp(ilat)   = ed_dens(1)
       dii_tmp(ilat)   = ed_docc(1)
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(Smats_tmp,Smats,Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Sreal,Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    !+- GET GLOC -+!
    Gmats=zero
    Greal=zero
    if(mpiID==0)write(LOGfile,*)"Get local GF:"
    if(mpiID==0)call start_timer
    do ik=1,Lk
       call get_gloc_mats(eloc,Smats,Gmats,epsik(ik),wt(ik))
       call get_gloc_real(eloc,Sreal,Greal,epsik(ik),wt(ik))
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)!,file="Glocal.eta")
    end do
    if(mpiID==0)call stop_timer
    !+- GET HYBRIDIZATION FUNCTION AND FIT BATHS -+!
    do ilat=1+mpiID,Nlat,mpiSIZE
       Hloc(1,1,1,1)=eloc(ilat)
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       do i=1,Lmats
          if(cg_scheme=='weiss')then
             Delta_tmp(ilat,1,1,i)  =   one/(one/Gmats(ilat,i) + Smats(ilat,i))
          else
             Delta_tmp(ilat,1,1,i) = xi*wm(i) + xmu - Hloc(1,1,1,1) - Smats(ilat,i) - one/Gmats(ilat,i)
          endif
       end do
       call chi2_fitgf(Delta_tmp(ilat,:,:,:),bath_(ilat,:,:),ispin=1)
    end do
    call MPI_ALLREDUCE(Delta_tmp(:,1,1,:),Delta,Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_solve_impurity_slab



  !+------------------------------------------------------------+!
  ! SC-version       
  !+------------------------------------------------------------+!
  subroutine ed_solve_sc_impurity_slab(bath_,eloc,epsik,wt,Delta,Gmats,Greal,Smats,Sreal,Usite)
    real(8)                          :: bath_(:,:,:),bath_tmp(size(bath_,2),size(bath_,3))
    real(8)                          :: eloc(Nlat)
    real(8),optional                 :: Usite(Nlat)
    real(8),dimension(:)             :: epsik
    real(8),dimension(size(epsik))   :: wt
    complex(8),intent(inout)         :: Delta(2,Nlat,Lmats)
    complex(8),intent(inout)         :: Gmats(2,Nlat,Lmats)
    complex(8),intent(inout)         :: Greal(2,Nlat,Lreal)
    complex(8),intent(inout)         :: Smats(2,Nlat,Lmats)
    complex(8),intent(inout)         :: Sreal(2,Nlat,Lreal)
    complex(8)                       :: Smats_tmp(2,Nlat,Lmats)
    complex(8)                       :: Sreal_tmp(2,Nlat,Lreal)
    complex(8)                       :: Gmats_k(2,Nlat,Lmats)
    complex(8)                       :: Greal_k(2,Nlat,Lreal)
    complex(8)                       :: Delta_tmp(2,Nlat,Norb,Norb,Lmats)
    complex(8)                       :: calG(2,Lmats),cdet
    integer                          :: ilat,i,ik
    integer                          :: Lk
    real(8)                          :: nii_tmp(Nlat),dii_tmp(Nlat),pii_tmp(Nlat)
    logical                          :: check_dim
    character(len=5)                 :: tmp_suffix
    if(size(bath_,1).ne.Nlat) stop "init_lattice_bath: wrong bath size dimension 3 (Nlat)"
    Lk=size(epsik)
    do ilat=1+mpiID,Nlat,mpiSIZE
       check_dim = check_bath_dimension(bath_(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    if(mpiID==0)call start_timer
    Smats = zero ; Smats_tmp = zero
    Sreal = zero ; Sreal_tmp = zero
    Delta = zero ; Delta_tmp = zero
    nii   = 0.d0 ; nii_tmp   = 0.d0
    dii   = 0.d0 ; dii_tmp   = 0.d0
    pii   = 0.d0 ; pii_tmp   = 0.d0
    !+- SOLVE SITE DEPENDENT IMPURITY PROBLEM -+!
    LOGfile = 800+mpiID
    if(mpiID==0) LOGfile=6
    do ilat=1+mpiID,Nlat,mpiSIZE
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       if(present(Usite)) Uloc(1)=Usite(ilat)
       Hloc(1,1,1,1)=eloc(ilat)
       call ed_solver(bath_(ilat,:,:))
       Smats_tmp(1,ilat,:) = impSmats(1,1,1,1,:)
       Smats_tmp(2,ilat,:) = impSAmats(1,1,1,1,:)
       Sreal_tmp(1,ilat,:) = impSreal(1,1,1,1,:)
       Sreal_tmp(2,ilat,:) = impSAreal(1,1,1,1,:)
       nii_tmp(ilat)   = ed_dens(1)
       dii_tmp(ilat)   = ed_docc(1)
       pii_tmp(ilat)   = ed_phisc(1)
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(Smats_tmp,Smats,2*Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(Sreal_tmp,Sreal,2*Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(pii_tmp,pii,Nlat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    !+- GET GLOC -+!
    Gmats=zero
    Greal=zero
    if(mpiID==0)call start_timer
    if(mpiID==0)write(LOGfile,*)"Get local GF (id=0):"
    do ik=1,Lk
       call get_sc_gloc_mats(eloc,Smats,Gmats,epsik(ik),wt(ik))
       call get_sc_gloc_real(eloc,Sreal,Greal,epsik(ik),wt(ik))
       if(mpiID==0)call eta(ik,Lk,unit=LOGfile)!,file="Glocal_mats.eta")
    end do
    if(mpiID==0)call stop_timer
    !+- GET HYBRIDIZATION FUNCTION AND FIT BATHS -+!
    do ilat=1+mpiID,Nlat,mpiSIZE
       Hloc(1,1,1,1)=eloc(ilat)
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       do i=1,Lmats
          if(cg_scheme=='weiss')then
             cdet       = abs(Gmats(1,ilat,i))**2 + (Gmats(2,ilat,i))**2
             calG(1,i) = conjg(Gmats(1,ilat,i))/cdet + Smats(1,ilat,i)
             calG(2,i) =  Gmats(2,ilat,i)/cdet + Smats(2,ilat,i) 
             cdet             =  abs(calG(1,i))**2 + (calG(2,i))**2
             Delta_tmp(1,ilat,1,1,i)  =  conjg(calG(1,i))/cdet
             Delta_tmp(2,ilat,1,1,i)  =  calG(2,i)/cdet
          else
             cdet       = abs(Gmats(1,ilat,i))**2 + (Gmats(2,ilat,i))**2
             Delta_tmp(1,ilat,1,1,i) = xi*wm(i) + xmu - Hloc(1,1,1,1) - Smats(1,ilat,i) - &
                  conjg(Gmats(1,ilat,i))/cdet 
             Delta_tmp(2,ilat,1,1,i) = -(Gmats(2,ilat,i)/cdet + Smats(2,ilat,i))
          endif
       end do
       call chi2_fitgf(Delta_tmp(:,ilat,:,:,:),bath_(ilat,:,:),ispin=1)
    end do
    call MPI_ALLREDUCE(Delta_tmp(:,:,1,1,:),Delta,2*Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine ed_solve_sc_impurity_slab









  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize lattice baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine init_lattice_baths(bath)
    real(8),dimension(:,:,:) :: bath
    integer :: ilat
    logical :: check_dim
    character(len=5) :: tmp_suffix
    if(size(bath,1).ne.Nlat) stop "init_lattice_bath: wrong bath size dimension 3 (Nlat)"
    do ilat=1,Nlat
       check_dim = check_bath_dimension(bath(ilat,:,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
       write(tmp_suffix,'(I4.4)') ilat
       ed_file_suffix="_site"//trim(tmp_suffix)
       call init_ed_solver(bath(ilat,:,:))
    end do
    call MPI_Barrier(MPI_COMM_WORLD,mpiERR)
  end subroutine init_lattice_baths




  !+-------------------------------------------------------------------+
  !PURPOSE  : given a bath array enforces the left-right symmetry 
  ! in finite and biased system by setting changing sign to the energies.
  !+-------------------------------------------------------------------+
  subroutine lr_symmetrize_bath(bath_)
    real(8),dimension(:,:,:)                      :: bath_
    type(effective_bath),dimension(size(bath_,1)) :: dmft_bath_
    integer                                       :: i,ilat
    if(size(bath_,1).ne.Nlat) stop "init_lattice_bath: wrong bath size dimension 3 (Nlat)"
    do ilat=1,Nlat
       call allocate_bath(dmft_bath_(ilat))
       call set_bath(bath_(ilat,:,:),dmft_bath_(ilat))
    enddo
    do ilat=1,Nlat
       do i=1,Nbath
          dmft_bath_(Nlat+1-ilat)%e(:,:,i)=-dmft_bath_(ilat)%e(:,:,i)
          dmft_bath_(Nlat+1-ilat)%v(:,:,i)= dmft_bath_(ilat)%v(:,:,i)
          if(ed_supercond)dmft_bath_(Nlat+1-ilat)%d(:,:,i)= dmft_bath_(ilat)%d(:,:,i)
       enddo
    enddo
    do ilat=1,Nlat
       call copy_bath(dmft_bath_(ilat),bath_(ilat,:,:))
       call deallocate_bath(dmft_bath_(ilat))
    enddo
  end subroutine lr_symmetrize_bath


end module RDMFT_WRAP_ED
