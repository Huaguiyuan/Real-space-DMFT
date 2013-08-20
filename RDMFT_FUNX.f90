!###############################################################
! PROGRAM  : RDMFT_FUNX
! PURPOSE  : Contains the main function performin RDMFT calculation
! for a given solver.
! main routines:
! - setup_<phase>_initial_sigma: setup the initial sigma to 
! start the R-DFMT calculation
! - get_<phase>_gloc_<formlism>_mpi: construct the local GF starting from 
! sigma function and tight-binding Hamiltonian via parallel matrix inversions
! - solve_<phase>_impurity_<formalism>_mpi: parallel solution of the 
! lattice sites impurity problems. This contains the call to specific
! solver routine. up to 20/08/2013 it contains only IPT interface.
!###############################################################
module RDMFT_FUNX
  USE RDMFT_VARS_GLOBAL
  !Impurity solver interface
  USE SOLVER_INTERFACE
  private

  interface symmetrize
     module procedure c_symmetrize,r_symmetrize
  end interface symmetrize

  interface reshuffled
     module procedure dv_reshuffled,zv_reshuffled,&
          dm_reshuffled,zm_reshuffled
  end interface reshuffled

  public :: get_tb_hamiltonian
  public :: setup_sc_initial_sigma
  public :: get_sc_gloc_matsu_mpi
  public :: get_sc_gloc_real_mpi
  public :: solve_sc_impurity_matsu_mpi
  public :: solve_sc_impurity_real_mpi


  public :: get_indip_list
  public :: symmetrize
  public :: reshuffled

contains

  !+----------------------------------------------------------------+
  !PURPOSE  : Build tight-binding Hamiltonian
  !+----------------------------------------------------------------+
  subroutine get_tb_hamiltonian(centered)
    integer          :: i,jj,j,k,row,col,link(4)
    logical,optional :: centered
    logical          :: symm
    symm=.false.;if(present(centered))symm=centered
    allocate(H0(Ns,Ns))
    allocate(icol(Ns),irow(Ns))
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


  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+----------------------------------------------------------------+
  !PURPOSE  : Setup/Read the initial sigma function to start DMFT loop
  !+----------------------------------------------------------------+
  subroutine setup_sc_initial_sigma(sigma)    
    complex(8),dimension(2,Ns,L) :: sigma
    real(8)                      :: foo_ome(L)
    logical                      :: check1,check2,check
    if(mpiID==0)then
       inquire(file=reg(fileSig),exist=check1)
       if(.not.check1)inquire(file=reg(fileSig)//".gz",exist=check1)
       inquire(file=reg(fileSelf),exist=check2)
       if(.not.check2)inquire(file=reg(fileSelf)//".gz",exist=check2)
       check=check1.AND.check2
       if(check)then
          call msg("Reading Self-energy from file:")
          call sread(reg(fileSig), foo_ome,sigma(1,1:Ns,1:L))
          call sread(reg(fileSelf),foo_ome,sigma(2,1:Ns,1:L))
       else
          call msg("Using Hartree-Fock-Bogoliubov self-energy")
          sigma(1,:,:)=zero ; sigma(2,:,:)=-deltasc
       endif
    endif
    call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  end subroutine setup_sc_initial_sigma


  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+----------------------------------------------------------------+
  !PURPOSE  : Evaluate the local GFs in Real-space formalism, 
  ! using matrix inversion
  !+----------------------------------------------------------------+
  subroutine get_sc_gloc_matsu_mpi(elocal,sigma,fg)
    real(8)    :: elocal(Ns)
    complex(8) :: fg(2,Ns,L),sigma(2,Ns,L)
    complex(8) :: Gloc(2*Ns,2*Ns),gf_tmp(2,Ns,L)
    integer    :: i,is
    call msg("Get local GF:",id=0)
    call start_timer
    fg=zero
    gf_tmp=zero
    do i=1+mpiID,L,mpiSIZE
       Gloc=zero
       Gloc(1:Ns,1:Ns)          = -H0
       Gloc(Ns+1:2*Ns,Ns+1:2*Ns)=  H0
       do is=1,Ns
          Gloc(is,is)      =  xi*wm(i)-sigma(1,is,i)        - elocal(is) + xmu
          Gloc(Ns+is,Ns+is)=  xi*wm(i)+conjg(sigma(1,is,i)) + elocal(is) - xmu !==-conjg(Gloc(is,is))
          Gloc(is,Ns+is)   = -sigma(2,is,i)
          Gloc(Ns+is,is)   = -sigma(2,is,i)!==sigma(2,is,L+1-i) a simmetry in Matsubara!
       enddo
       call matrix_inverse_sym(Gloc)
       forall(is=1:Ns)
          gf_tmp(1,is,i) = Gloc(is,is)
          !##ACTHUNG!!
          gf_tmp(2,is,i) = dreal(Gloc(is,Ns+is))
       end forall
       call eta(i,L,file="Glocal.eta")
    enddo
    call stop_timer
    call MPI_ALLREDUCE(gf_tmp,fg,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_sc_gloc_matsu_mpi

  subroutine get_sc_gloc_real_mpi(elocal,sigma,fg)
    real(8)    :: elocal(Ns)
    complex(8) :: fg(2,Ns,L),sigma(2,Ns,L)
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
          zeta1=        cmplx(wr(i),eps,8)     + xmu - sigma(1,is,i)       - elocal(is)
          zeta2=-conjg( cmplx(wr(L+1-i),eps,8) + xmu - sigma(1,is,L+1-i) ) + elocal(is)
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
    call MPI_ALLREDUCE(gf_tmp,fg,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    ! call MPI_REDUCE(gf_tmp,fg,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    ! call MPI_BCAST(fg,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_sc_gloc_real_mpi



  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+----------------------------------------------------------------+
  !PURPOSE  : parallel solution of impurity problem SC-version
  !+----------------------------------------------------------------+
  subroutine solve_sc_impurity_matsu_mpi(fg,sigma)
    integer    :: is,i
    complex(8) :: fg(2,Ns,L),sigma(2,Ns,L)
    real(8)    :: nii_tmp(Ns),dii_tmp(Ns)
    complex(8) :: sigma_tmp(2,Ns,L)
    call msg("Solve impurity:")
    call start_timer
    sigma_tmp=zero
    nii_tmp  =0.d0
    dii_tmp  =0.d0
    do is=1+mpiID,Ns,mpiSIZE
       call ipt_matsu_solve_per_site(is,fg,sigma,sigma_tmp(:,is,:),nii_tmp(is),dii_tmp(is))
       call eta(is,Ns,file="Impurity.eta")
    enddo
    call stop_timer
    call MPI_ALLREDUCE(sigma_tmp,sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine solve_sc_impurity_matsu_mpi

  subroutine solve_sc_impurity_real_mpi(fg,sigma)
    integer    :: is,i
    complex(8) :: fg(2,Ns,L),sigma(2,Ns,L)
    real(8)    :: nii_tmp(Ns),dii_tmp(Ns)
    complex(8) :: sigma_tmp(2,Ns,L)
    call msg("Solve impurity:")
    call start_timer
    sigma_tmp=zero
    nii_tmp  =0.d0
    dii_tmp  =0.d0
    do is=1+mpiID,Ns,mpiSIZE
       call ipt_real_solve_per_site(is,fg,sigma,sigma_tmp(:,is,:),nii_tmp(is),dii_tmp(is))
       call eta(is,Ns,file="Impurity.eta")
    enddo
    call stop_timer
    call MPI_ALLREDUCE(sigma_tmp,sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(nii_tmp,nii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    call MPI_ALLREDUCE(dii_tmp,dii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  end subroutine solve_sc_impurity_real_mpi



  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+----------------------------------------------------------------+
  !PURPOSE  : IPT solution of the is^th-impurity problem.
  !+----------------------------------------------------------------+
  subroutine ipt_matsu_solve_per_site(is,fg,sigma,sig_tmp,ntmp,dtmp)
    integer,intent(in)                           :: is
    complex(8),dimension(2,Ns,L),intent(in)      :: fg,sigma !Local GF,Self-energy
    complex(8),dimension(2,L),intent(out)        :: sig_tmp
    real(8),intent(out)                          :: ntmp,dtmp
    complex(8)                                   :: det(L)
    complex(8),dimension(:,:,:),allocatable,save :: Wold
    complex(8),dimension(2,L)                    :: calG,fg0
    real(8),dimension(2,0:L)                     :: fgt,fg0t
    real(8)                                      :: n,n0,delta,delta0
    if(.not.allocated(Wold))allocate(Wold(2,Ns,L))
    if(mix_type==1)Wold(:,is,:) = sigma(:,is,:)
    !
    call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
    call fftgf_iw2tau(fg(2,is,:),fgt(2,0:L),beta,notail=.true.)
    n = -fgt(1,L) ; delta = -u*fgt(2,L)
    !
    ntmp=2.d0*n; dtmp=delta
    !
    fg0=zero ; calG=zero
    det       = abs(fg(1,is,:))**2    + fg(2,is,:)**2
    fg0(1,:)  = conjg(fg(1,is,:))/det + sigma(1,is,:) + U*(n-0.5d0)
    fg0(2,:)  = fg(2,is,:)/det        + sigma(2,is,:) + delta
    det       = abs(fg0(1,:))**2      + fg0(2,:)**2
    calG(1,:) = conjg(fg0(1,:))/det
    calG(2,:) = fg0(2,:)/det
    !
    if(mix_type==0)then
       if(iloop>1)calG(:,:) =  weight*calG(:,:) + (1.d0-weight)*Wold(:,is,:)
       Wold(:,is,:)  =  calG(:,:)
    endif
    !
    call fftgf_iw2tau(calG(1,:),fg0t(1,:),beta)
    call fftgf_iw2tau(calG(2,:),fg0t(2,:),beta,notail=.true.)
    n0=-fg0t(1,L) ; delta0= -u*fg0t(2,L)
    write(750,"(2I4,4(f16.12))",advance="yes")mpiID,is,n,n0,delta,delta0
    !
    sig_tmp(:,:) =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)
    if(mix_type==1)sig_tmp(:,:) =  weight*sig_tmp(:,:) + (1.d0-weight)*Wold(:,is,:)
    !!##ACTHUNG!!!
    sig_tmp(2,:) = dreal(sig_tmp(2,:))+zero
  end subroutine ipt_matsu_solve_per_site

  subroutine ipt_real_solve_per_site(is,fg,sigma,sig_tmp,ntmp,dtmp)
    integer,intent(in)                           :: is
    complex(8),dimension(2,Ns,L),intent(in)      :: fg,sigma !Local GF,Self-energy
    complex(8),dimension(2,L),intent(out)        :: sig_tmp
    real(8),intent(out)                          :: ntmp,dtmp
    complex(8)                                   :: det
    complex(8),dimension(:,:,:),allocatable,save :: Wold
    complex(8),dimension(2,1:L)                  :: calG,fg0
    real(8)                                      :: n,n0,delta,delta0
    if(.not.allocated(Wold))allocate(Wold(2,Ns,1:L))
    if(mix_type==1)Wold(:,is,:) = sigma(:,is,:)
    !
    n    = -sum(dimag(fg(1,is,:))*fermi(wr,beta))*fmesh/pi ! densita' per spin
    delta= -u*sum(dimag(fg(2,is,:))*fermi(wr,beta))*fmesh/pi
    !
    ntmp=2.d0*n; dtmp=delta
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
    if(mix_type==0)then
       if(iloop>1)calG(:,:) =  weight*calG(:,:) + (1.d0-weight)*Wold(:,is,:)
       Wold(:,is,:)  =  calG(:,:)
    endif
    !
    n0    = -sum(dimag(calG(1,:))*fermi(wr,beta))*fmesh/pi
    delta0= -u*sum(dimag(calG(2,:))*fermi(wr,beta))*fmesh/pi
    write(750,"(2I4,4(f16.12))",advance="yes")mpiID,is,n,n0,delta,delta0
    !
    sig_tmp(:,:)  =  solve_mpt_sc_sopt(calG,wr,n,n0,delta,delta0,L)
    ! sigma_tmp(:,is,:) =  solve_mpt_sc_sopt(calG,wr,n,n0,delta,delta0,L)
    if(mix_type==1)sig_tmp(:,:) =  weight*sig_tmp(:,:) + (1.d0-weight)*Wold(:,is,:)
    !
  end subroutine ipt_real_solve_per_site






  !******************************************************************
  !******************************************************************
  !
  ! TRAP RELATED FUNCTIONS:
  !
  !******************************************************************
  !******************************************************************


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



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+----------------------------------------------------------------+
  !PURPOSE : implement the trap simmetries on a real vector variable 
  ! with Ns components 
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



  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+----------------------------------------------------------------+
  !PURPOSE : Reshuffle the Ns Lattice sites into Nindependent one
  !given the map indipsites
  !+----------------------------------------------------------------+
  function dv_reshuffled(m_in) result(m_out)
    integer                               :: i
    real(8), dimension(Ns)           :: m_in
    real(8), dimension(Nindip)       :: m_out
    do i=1,Nindip
       m_out(i)=m_in(indipsites(i))
    enddo
  end function dv_reshuffled
  !
  function zv_reshuffled(m_in) result(m_out)
    integer                               :: i
    complex(8), dimension(Ns)           :: m_in
    complex(8), dimension(Nindip)       :: m_out
    do i=1,Nindip
       m_out(i)=m_in(indipsites(i))
    enddo
  end function zv_reshuffled
  !
  function dm_reshuffled(m_in) result(m_out)
    integer                               :: i
    real(8), dimension(Ns,L)           :: m_in
    real(8), dimension(Nindip,L)       :: m_out
    do i=1,Nindip
       m_out(i,:)=m_in(indipsites(i),:)
    enddo
  end function dm_reshuffled
  !
  function zm_reshuffled(m_in) result(m_out)
    integer                               :: i
    complex(8), dimension(Ns,L)           :: m_in
    complex(8), dimension(Nindip,L)       :: m_out
    do i=1,Nindip
       m_out(i,:)=m_in(indipsites(i),:)
    enddo
  end function zm_reshuffled
  !



  !******************************************************************
  !******************************************************************
  !******************************************************************

end module RDMFT_FUNX
