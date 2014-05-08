!###############################################################
! PROGRAM  : RDMFT_FUNX
! PURPOSE  : Contains the main function performin RDMFT calculation
! for a given solver.
!###############################################################
module RDMFT_AUX_FUNX
  USE RDMFT_INPUT_VARS
  USE RDMFT_VARS_GLOBAL
  USE TIMER
  USE IOTOOLS,   only:reg,sread
  USE MATRIX,    only:matrix_inverse,matrix_inverse_sym
  implicit none
  private

  public :: get_gloc_mats
  public :: get_gloc_real
  public :: get_sc_gloc_mats
  public :: get_sc_gloc_real

  !
  !NORMAL
  !
  interface get_gloc_mats
     module procedure get_gloc_mats_lattice,get_gloc_mats_slab
  end interface get_gloc_mats

  interface get_gloc_real
     module procedure get_gloc_real_lattice,get_gloc_real_slab
  end interface get_gloc_real

  !
  !SUPERCOND
  !
  interface get_sc_gloc_mats
     module procedure get_sc_gloc_mats_lattice,get_sc_gloc_mats_slab
  end interface get_sc_gloc_mats

  interface get_sc_gloc_real
     module procedure get_sc_gloc_real_lattice,get_sc_gloc_real_slab
  end interface get_sc_gloc_real


  public :: get_tb_hamiltonian  !now only 2D square
  public :: get_slab_hamiltonian !to be overloaded with get_tb_hamiltonian !!!
  public :: setup_initial_sigma, setup_sc_initial_sigma
  public :: get_indip_list


  public :: symmetrize
  public :: reshuffled

  interface symmetrize
     module procedure c_symmetrize,r_symmetrize
  end interface symmetrize

  interface reshuffled
     module procedure dv_reshuffled,zv_reshuffled,&
          dm_reshuffled,zm_reshuffled
  end interface reshuffled


contains



  !+----------------------------------------------------------------+
  !PURPOSE  : Evaluate the local GFs in Real-space formalism, 
  ! using matrix inversion 
  !+----------------------------------------------------------------+
  !
  !NORMAL
  !
  subroutine get_gloc_mats_lattice(elocal,sigma,fg) 
    real(8)    :: elocal(Nlat)
    complex(8) :: fg(Nlat,Lmats),sigma(Nlat,Lmats)
    complex(8) :: zeta,Gloc(Nlat,Nlat),gf_tmp(Nlat,1:Lmats)
    integer    :: i,is
    if(mpiID==0)write(LOGfile,*)"Get local GF Mats (id=0):"
    if(mpiID==0)call start_timer
    gf_tmp=zero
    fg=zero
    do i=1+mpiID,Lmats,mpiSIZE
       zeta  = xi*wm(i) + xmu
       Gloc  = -H0
       do is=1,Nlat
          Gloc(is,is) = Gloc(is,is) + zeta - elocal(is) - sigma(is,i) 
       enddo
       call matrix_inverse_sym(Gloc)
       do is=1,Nlat
          gf_tmp(is,i) = Gloc(is,is)
       enddo
       if(mpiID==0)call eta(i,Lmats,file="Glocal_mats.eta")
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(gf_tmp,fg,Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_gloc_mats_lattice

  subroutine get_gloc_real_lattice(elocal,sigma,fg) 
    real(8)    :: elocal(Nlat)
    complex(8) :: fg(Nlat,Lreal),sigma(Nlat,Lreal)
    complex(8) :: zeta,Gloc(Nlat,Nlat),gf_tmp(Nlat,1:Lreal)
    integer    :: i,is
    if(mpiID==0)write(LOGfile,*)"Get local GF Real (id=0):"
    if(mpiID==0)call start_timer
    gf_tmp=zero 
    fg=zero
    do i=1+mpiID,Lreal,mpiSIZE
       zeta  = cmplx(wr(i),eps,8) + xmu
       Gloc  = zero-H0
       do is=1,Nlat
          Gloc(is,is)=Gloc(is,is) + zeta - sigma(is,i) - elocal(is)
       enddo
       call matrix_inverse_sym(Gloc)
       do is=1,Nlat
          gf_tmp(is,i) = Gloc(is,is)
       enddo
       if(mpiID==0)call eta(i,Lreal,file="Glocal_real.eta")
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(gf_tmp,fg,Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_gloc_real_lattice


  !
  !SUPERCOND
  !
  subroutine get_sc_gloc_mats_lattice(elocal,sigma,fg)
    real(8)    :: elocal(Nlat)
    complex(8) :: fg(2,Nlat,Lmats),sigma(2,Nlat,Lmats)
    complex(8) :: Gloc(2*Nlat,2*Nlat),gf_tmp(2,Nlat,Lmats)
    integer    :: i,is
    if(mpiID==0)write(LOGfile,*)"Get local GF Mats (id=0):"
    if(mpiID==0)call start_timer
    fg=zero
    gf_tmp=zero
    do i=1+mpiID,Lmats,mpiSIZE
       Gloc=zero
       Gloc(1:Nlat,1:Nlat)          = -H0
       Gloc(Nlat+1:2*Nlat,Nlat+1:2*Nlat)=  H0
       do is=1,Nlat
          Gloc(is,is)           =  xi*wm(i)-sigma(1,is,i)        - elocal(is) + xmu
          Gloc(Nlat+is,Nlat+is) =  xi*wm(i)+conjg(sigma(1,is,i)) + elocal(is) - xmu !==-conjg(Gloc(is,is))
          Gloc(is,Nlat+is)      = -sigma(2,is,i)
          Gloc(Nlat+is,is)      = -sigma(2,is,i)                                    !==sigma(2,is,L+1-i) a simmetry in Matsubara!
       enddo
       call matrix_inverse_sym(Gloc)
       forall(is=1:Nlat)
          gf_tmp(1,is,i) = Gloc(is,is)
          !##ACTHUNG!!
          gf_tmp(2,is,i) = dreal(Gloc(is,Nlat+is))
       end forall
       if(mpiID==0)call eta(i,Lmats,file="Glocal_mats.eta")
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(gf_tmp,fg,2*Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_sc_gloc_mats_lattice

  subroutine get_sc_gloc_real_lattice(elocal,sigma,fg)
    real(8)    :: elocal(Nlat)
    complex(8) :: fg(2,Nlat,Lreal),sigma(2,Nlat,Lreal)
    complex(8) :: Gloc(2*Nlat,2*Nlat),gf_tmp(2,Nlat,Lreal),zeta1,zeta2
    integer    :: i,is
    if(mpiID==0)write(LOGfile,*)"Get local GF Real (id=0):"
    if(mpiID==0)call start_timer
    fg=zero ; gf_tmp=zero
    do i=1+mpiID,Lreal,mpiSIZE
       Gloc=zero
       Gloc(1:Nlat,1:Nlat)          = -H0
       Gloc(Nlat+1:2*Nlat,Nlat+1:2*Nlat)=  H0
       do is=1,Nlat
          zeta1=        cmplx(wr(i),eps,8)     + xmu - sigma(1,is,i)       - elocal(is)
          zeta2=-conjg( cmplx(wr(Lreal+1-i),eps,8) + xmu - sigma(1,is,Lreal+1-i) ) + elocal(is)
          Gloc(is,is)      = zeta1
          Gloc(Nlat+is,Nlat+is)= zeta2
          Gloc(is,Nlat+is)   = -sigma(2,is,i)
          !S_12(w)=S^*(-w), by symmetry =S(w)=S_12(w)
          !we set this block to the correct function S^*(-w)
          !anyway with respect to the next call to *symmetric* 
          !matrix inversion this position is irrelevant 
          !as the routine only consider the upper blocks triangles (11,12,22)
          Gloc(Nlat+is,is)   = -conjg(sigma(2,is,Lreal+1-i))
       enddo
       !the call to *symmetry* routine enforces the symmetry condition
       !S^*(-w)=S(w)
       !this condition can be numerically broken using the generic 
       !inversion routine generating small (though sizeable) errors. 
       call matrix_inverse_sym(Gloc)
       forall(is=1:Nlat)
          gf_tmp(1,is,i) = Gloc(is,is)
          gf_tmp(2,is,i) = Gloc(is,Nlat+is)
       end forall
       if(mpiID==0)call eta(i,Lreal,file="Glocal_real.eta")
    enddo
    if(mpiID==0)call stop_timer
    call MPI_ALLREDUCE(gf_tmp,fg,2*Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_sc_gloc_real_lattice






  !+----------------------------------------------------------------+
  !PURPOSE  : Evaluate the local GFs in Real-space formalism, 
  ! using matrix inversion 
  !+----------------------------------------------------------------+
  !
  !NORMAL
  !
  subroutine get_gloc_mats_slab(elocal,sigma,fg,ek,wtk)
    real(8)                  :: elocal(Nlat)
    complex(8),intent(inout) :: fg(Nlat,Lmats),sigma(Nlat,Lmats)
    real(8)                  :: ek,wtk
    complex(8)               :: zeta,Gloc(Nlat,Nlat),gf_tmp(Nlat,Lmats),fg_k(Nlat,Lmats)
    integer                  :: i,is
    fg_k=zero
    gf_tmp=zero
    do i=1+mpiID,Lmats,mpiSIZE
       Gloc  = zero
       zeta  = xi*wm(i) + xmu
       Gloc  = -H0
       do is=1,Nlat
          Gloc(is,is)           =  zeta  - elocal(is) - ek - sigma(is,i)
       enddo
       call matrix_inverse_sym(Gloc)
       forall(is=1:Nlat)
          gf_tmp(is,i) = Gloc(is,is)
       end forall
    enddo
    !+- reduce gf_tmp on fg_k -+!
    call MPI_ALLREDUCE(gf_tmp,fg_k,Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
    !+- k-sums -+!
    fg = fg + fg_k*wtk
  end subroutine get_gloc_mats_slab

  subroutine get_gloc_real_slab(elocal,sigma,fg,ek,wtk)
    real(8)                  :: elocal(Nlat)
    complex(8),intent(inout) :: fg(Nlat,Lreal),sigma(Nlat,Lreal)
    real(8)                  :: ek,wtk
    complex(8)               :: zeta,Gloc(Nlat,Nlat),gf_tmp(Nlat,Lreal),fg_k(Nlat,Lreal)
    integer                  :: i,is
    fg_k=zero
    gf_tmp=zero
    do i=1+mpiID,Lreal,mpiSIZE
       Gloc  = zero
       zeta  = cmplx(wr(i),eps,8) + xmu
       Gloc  = zero-H0
       do is=1,Nlat
          Gloc(is,is) =  zeta  - elocal(is) - ek - sigma(is,i)
       enddo
       call matrix_inverse_sym(Gloc)
       forall(is=1:Nlat)
          gf_tmp(is,i) = Gloc(is,is)
       end forall
    enddo
    !+- reduce gf_tmp on fg_k -+!
    call MPI_ALLREDUCE(gf_tmp,fg_k,Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
    !+- k-sums -+!
    fg = fg + fg_k*wtk
  end subroutine get_gloc_real_slab



  !
  !SUPERCOND
  !
  subroutine get_sc_gloc_mats_slab(elocal,sigma,fg,ek,wtk)
    real(8)    :: elocal(Nlat)
    complex(8),intent(inout) :: fg(2,Nlat,Lmats),sigma(2,Nlat,Lmats)
    real(8)    :: ek,wtk
    complex(8) :: Gloc(2*Nlat,2*Nlat),gf_tmp(2,Nlat,Lmats),fg_k(2,Nlat,Lmats)
    integer    :: i,is
    fg_k=zero
    gf_tmp=zero
    do i=1+mpiID,Lmats,mpiSIZE
       Gloc=zero
       Gloc(1:Nlat,1:Nlat)          = -H0
       Gloc(Nlat+1:2*Nlat,Nlat+1:2*Nlat)=  H0
       do is=1,Nlat
          Gloc(is,is)           =  xi*wm(i) - ek - sigma(1,is,i)        - elocal(is) + xmu
          Gloc(Nlat+is,Nlat+is) =  xi*wm(i) + ek + conjg(sigma(1,is,i)) + elocal(is) - xmu !==-conjg(Gloc(is,is))
          Gloc(is,Nlat+is)      = -sigma(2,is,i)
          Gloc(Nlat+is,is)      = -sigma(2,is,i)                                          !==sigma(2,is,L+1-i) a simmetry in Matsubara!
       enddo
       call matrix_inverse_sym(Gloc)
       forall(is=1:Nlat)
          gf_tmp(1,is,i) = Gloc(is,is)
          !##ACTHUNG!!
          gf_tmp(2,is,i) = dreal(Gloc(is,Nlat+is))
       end forall
    enddo
    !+- reduce gf_tmp on fg_k -+!
    call MPI_ALLREDUCE(gf_tmp,fg_k,2*Nlat*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
    !+- k-sums -+!
    fg = fg + fg_k*wtk
  end subroutine get_sc_gloc_mats_slab

  subroutine get_sc_gloc_real_slab(elocal,sigma,fg,ek,wtk)
    real(8)    :: elocal(Nlat)
    complex(8) :: fg(2,Nlat,Lreal),sigma(2,Nlat,Lreal)
    real(8)    :: ek,wtk
    complex(8) :: Gloc(2*Nlat,2*Nlat),gf_tmp(2,Nlat,Lreal),zeta1,zeta2,fg_k(2,Nlat,Lreal)
    integer    :: i,is
    fg_k=zero; 
    gf_tmp=zero
    do i=1+mpiID,Lreal,mpiSIZE
       Gloc=zero
       Gloc(1:Nlat,1:Nlat)          = -H0
       Gloc(Nlat+1:2*Nlat,Nlat+1:2*Nlat)=  H0
       do is=1,Nlat
          zeta1=        cmplx(wr(i),eps,8) - ek     + xmu - sigma(1,is,i)       - elocal(is)
          zeta2=-conjg( cmplx(wr(Lreal+1-i),eps,8) - ek + xmu - sigma(1,is,Lreal+1-i) ) + elocal(is)
          Gloc(is,is)      = zeta1
          Gloc(Nlat+is,Nlat+is)= zeta2
          Gloc(is,Nlat+is)   = -sigma(2,is,i)
          !S_12(w)=S^*(-w), by symmetry =S(w)=S_12(w)
          !we set this block to the correct function S^*(-w)
          Gloc(Nlat+is,is)   = -conjg(sigma(2,is,Lreal+1-i))
       enddo
       !the call to *symmetry* routine enforces the symmetry condition
       !S^*(-w)=S(w)
       call matrix_inverse_sym(Gloc)
       forall(is=1:Nlat)
          gf_tmp(1,is,i) = Gloc(is,is)
          gf_tmp(2,is,i) = Gloc(is,Nlat+is)
       end forall
    enddo
    call MPI_ALLREDUCE(gf_tmp,fg_k,2*Nlat*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
    fg = fg + fg_k*wtk
  end subroutine get_sc_gloc_real_slab











  !+----------------------------------------------------------------+
  !PURPOSE  : Build tight-binding Hamiltonian
  !+----------------------------------------------------------------+
  subroutine get_tb_hamiltonian(centered)
    integer          :: i,jj,j,k,row,col,link(4)
    logical,optional :: centered
    logical          :: symm
    symm=.false.;if(present(centered))symm=centered
    allocate(H0(Nlat,Nlat))
    allocate(icol(Nlat),irow(Nlat))
    allocate(ij2site(Nside,Nside))
    H0=0.d0
    do row=0,Nside-1
       do col=0,Nside-1
          i=col+row*Nside+1
          if(.not.symm)then
             irow(i)=row+1
             icol(i)=col+1
             ij2site(row+1,col+1)=i
          else
             irow(i)=-Nside/2+row                     ! cambio la tabella i -> isite,jsite
             icol(i)=-Nside/2+col                     ! per farla simmetrica.. aiutera' 
             ij2site(row-Nside/2,col-Nside/2)=i       ! a implementare le simmetrie
          endif
          !
          if(pbcflag)then ! PBC are implemented using the state labels and so they are mpt affected by symm
             !HOPPING w/ PERIODIC BOUNDARY CONDITIONS
             link(1)= row*Nside+1              + mod(col+1,Nside)  ;
             link(3)= row*Nside+1              + (col-1)           ; if((col-1)<0)link(3)=(Nside+(col-1))+row*Nside+1
             link(2)= mod(row+1,Nside)*Nside+1 + col               ; 
             link(4)= (row-1)*Nside+1          + col               ; if((row-1)<0)link(4)=col+(Nside+(row-1))*Nside+1
          else   
             !without PBC
             link(1)= row*Nside+1              + col+1   ; if((col+1)==Nside)link(1)=0
             link(3)= row*Nside+1              +(col-1)  ; if((col-1)<0)     link(3)=0
             link(2)= (row+1)*Nside+1 + col              ; if((row+1)==Nside)link(2)=0
             link(4)= (row-1)*Nside+1          + col     ; if((row-1)<0)     link(4)=0
          endif
          do jj=1,4
             if(link(jj)>0)H0(i,link(jj))=-ts !! ts must be negative.
          enddo
       enddo
    enddo
  end subroutine get_tb_hamiltonian


  !+- bulid linear chain hamiltonian -+!
  subroutine get_slab_hamiltonian
    integer          :: i
    allocate(H0(Nlat,Nlat))
    H0 = 0.d0
    do i=1,Nlat-1
       H0(i,i+1)=-ts
       H0(i+1,i)=-ts
    end do
  end subroutine get_slab_hamiltonian






  !+----------------------------------------------------------------+
  !PURPOSE  : Setup/Read the initial sigma function to start DMFT loop
  !+----------------------------------------------------------------+
  subroutine setup_initial_sigma(sigma)    
    complex(8),dimension(:,:) :: sigma
    real(8)                   :: foo_ome(size(sigma,2))
    logical                   :: check1,check2,check
    integer                   :: Lf,Ni
    Ni=size(sigma,1);if(Ni/=Nlat)stop"setup_initial_sigma: error in Sigma size 1"
    Lf=size(sigma,2);if(Lf/=Lmats.OR.Lf/=Lreal)stop"setup_initial_sigma: error in Sigma size 2"
    if(mpiID==0)then
       inquire(file=reg(fileSig),exist=check1)
       if(.not.check1)inquire(file=reg(fileSig)//".gz",exist=check1)
       check=check1
       if(check)then
          if(mpiID==0)write(LOGfile,*)"Reading Self-energy from file:"
          call sread(reg(fileSig),foo_ome,sigma(:,:))
       else
          if(mpiID==0)write(LOGfile,*)"Using Hartree-Fock-Bogoliubov self-energy"
          sigma(:,:)=zero
       endif
    endif
    call MPI_BCAST(sigma,Nlat*Lf,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  end subroutine setup_initial_sigma

  subroutine setup_sc_initial_sigma(sigma)    
    complex(8),dimension(:,:,:)  :: sigma
    real(8)                      :: foo_ome(size(sigma,3))
    logical                      :: check1,check2,check
    integer                     :: Lf,Ni
    if(size(sigma,1)/=2)stop"setup_sc_initial_sigma: error in Sigma size 1"
    Ni=size(sigma,2);if(Ni/=Nlat)stop"setup_sc_initial_sigma: error in Sigma size 2"
    Lf=size(sigma,3);if(Lf/=Lmats.OR.Lf/=Lreal)stop"setup_sc_initial_sigma: error in Sigma size 3"
    if(mpiID==0)then
       inquire(file=reg(fileSig),exist=check1)
       if(.not.check1)inquire(file=reg(fileSig)//".gz",exist=check1)
       inquire(file=reg(fileSelf),exist=check2)
       if(.not.check2)inquire(file=reg(fileSelf)//".gz",exist=check2)
       check=check1.AND.check2
       if(check)then
          if(mpiID==0)write(LOGfile,*)"Reading Self-energy from file:"
          call sread(reg(fileSig), foo_ome,sigma(1,:,:))
          call sread(reg(fileSelf),foo_ome,sigma(2,:,:))
       else
          if(mpiID==0)write(LOGfile,*)"Using Hartree-Fock-Bogoliubov self-energy"
          sigma(1,:,:)=zero ; sigma(2,:,:)=-deltasc
       endif
    endif
    call MPI_BCAST(sigma,2*Nlat*Lf,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  end subroutine setup_sc_initial_sigma







  !+----------------------------------------------------------------+
  !PURPOSE : build the list of the indipendent sites (1/8 of the square)
  !+----------------------------------------------------------------+
  subroutine get_indip_list()
    integer :: i,row,col,istate,jstate
    i=0
    do col=0,Nside/2
       do row=0,col
          i= i+1
          indipsites(i)=ij2site(row,col)
       enddo
    enddo
  end subroutine get_indip_list





  !+----------------------------------------------------------------+
  !PURPOSE : implement the trap simmetries on a real vector variable 
  ! with Nlat components 
  !+----------------------------------------------------------------+
  subroutine r_symmetrize(vec)
    integer                             :: row,col
    real(8), dimension(:),intent(INOUT) :: vec
    !assi cartesiani e diagonale  degeneracy=4
    do col=1,Nside/2
       vec(ij2site( col,   0))   =vec(ij2site(0,  col))
       vec(ij2site( 0,  -col))   =vec(ij2site(0,  col))
       vec(ij2site(-col,   0))   =vec(ij2site(0,  col))
       vec(ij2site(col, -col))   =vec(ij2site(col,col))
       vec(ij2site(-col, col))   =vec(ij2site(col,col))
       vec(ij2site(-col,-col))   =vec(ij2site(col,col))
    enddo
    !nel semipiano e fuori dalle linee sopramenzionate degeneracy =8 
    do col=2,Nside/2    
       do row=1,col-1
          vec(ij2site(-row, col))  =vec(ij2site(row,col)) ! riflessioni rispetto agli assi
          vec(ij2site( row,-col))  =vec(ij2site(row,col))
          vec(ij2site(-row,-col))  =vec(ij2site(row,col))
          vec(ij2site( col, row))  =vec(ij2site(row,col)) ! riflessione con la bisettrice 
          vec(ij2site(-col, row))  =vec(ij2site(row,col))
          vec(ij2site( col,-row))  =vec(ij2site(row,col))
          vec(ij2site(-col,-row))  =vec(ij2site(row,col))
       enddo
    enddo
  end subroutine r_symmetrize
  !+----------------------------------------------------------------+
  subroutine c_symmetrize(vec)
    integer                                :: row,col
    complex(8), dimension(:),intent(INOUT) :: vec
    !assi cartesiani e diagonale  degeneracy=4
    do col=1,Nside/2
       vec(ij2site( col,   0))   =vec(ij2site(0,  col))
       vec(ij2site( 0,  -col))   =vec(ij2site(0,  col))
       vec(ij2site(-col,   0))   =vec(ij2site(0,  col))
       vec(ij2site(col, -col))   =vec(ij2site(col,col))
       vec(ij2site(-col, col))   =vec(ij2site(col,col))
       vec(ij2site(-col,-col))   =vec(ij2site(col,col))
    enddo
    !nel semipiano e fuori dalle linee sopramenzionate degeneracy =8 
    do col=2,Nside/2    
       do row=1,col-1
          vec(ij2site(-row, col))  =vec(ij2site(row,col)) ! riflessioni rispetto agli assi
          vec(ij2site( row,-col))  =vec(ij2site(row,col))
          vec(ij2site(-row,-col))  =vec(ij2site(row,col))
          vec(ij2site( col, row))  =vec(ij2site(row,col)) ! riflessione con la bisettrice 
          vec(ij2site(-col, row))  =vec(ij2site(row,col))
          vec(ij2site( col,-row))  =vec(ij2site(row,col))
          vec(ij2site(-col,-row))  =vec(ij2site(row,col))
       enddo
    enddo
  end subroutine c_symmetrize






  !+----------------------------------------------------------------+
  !PURPOSE : Reshuffle the Nlat Lattice sites into Nindependent one
  !given the map indipsites
  !+----------------------------------------------------------------+
  function dv_reshuffled(m_in) result(m_out)
    integer                               :: i
    real(8), dimension(Nlat)           :: m_in
    real(8), dimension(Nindip)       :: m_out
    do i=1,Nindip
       m_out(i)=m_in(indipsites(i))
    enddo
  end function dv_reshuffled
  !
  function zv_reshuffled(m_in) result(m_out)
    integer                               :: i
    complex(8), dimension(Nlat)           :: m_in
    complex(8), dimension(Nindip)       :: m_out
    do i=1,Nindip
       m_out(i)=m_in(indipsites(i))
    enddo
  end function zv_reshuffled
  !
  function dm_reshuffled(m_in) result(m_out)
    integer                                 :: i
    real(8), dimension(:,:)                 :: m_in
    real(8), dimension(Nindip,size(m_in,2)) :: m_out
    if(size(m_in,1)<Nindip)stop"dm_reshuffled: error m_in size 1"
    do i=1,Nindip
       m_out(i,:)=m_in(indipsites(i),:)
    enddo
  end function dm_reshuffled
  !
  function zm_reshuffled(m_in) result(m_out)
    integer                                    :: i
    complex(8), dimension(:,:)                 :: m_in
    complex(8), dimension(Nindip,size(m_in,2)) :: m_out
    if(size(m_in,1)<Nindip)stop"dm_reshuffled: error m_in size 1"
    do i=1,Nindip
       m_out(i,:)=m_in(indipsites(i),:)
    enddo
  end function zm_reshuffled
  !

end module RDMFT_AUX_FUNX
