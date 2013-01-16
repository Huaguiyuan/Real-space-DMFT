!########################################################
!     PURPOSE  :solve the attractive (A) Trapped (T) Hubbard
! model (HM) using Modified Perturbation THeory (MPT), w/ DMFT
!     AUTHORS  : A.Amaricci, A.Privitera (CNR-IOM)
!########################################################

program ahm_real_trap
  USE RDMFT_VARS_GLOBAL
  implicit none
  real(8)                                 :: r
  real(8)                                 :: n_tot,delta_tot
  integer                                 :: is,esp,lm
  logical                                 :: converged,convergedN,convergedD
  complex(8),allocatable,dimension(:,:)   :: fg0,fg_0
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma,sigma_tmp
  real(8),allocatable,dimension(:)        :: nii_tmp,dii_tmp,gap_ii_tmp

  real(8),allocatable,dimension(:) :: acheck

  !GLOBAL INITIALIZATION:
  !===============================================================
  include "init_global_trap.f90" 


  !ALLOCATE WORKING ARRAYS
  !===============================================================
  allocate(wr(L))
  wr = linspace(-wmax,wmax,L,mesh=fmesh)

  allocate(fg_0(Ns,L)) ! noninteracting Green Function ! only spin up !
  allocate(fg(2,Ns,L))
  allocate(fg0(2,L))
  allocate(sigma(2,Ns,L))
  !
  allocate(sigma_tmp(2,Ns,L))
  allocate(nii_tmp(Ns),dii_tmp(Ns),gap_ii_tmp(Ns))
  if(symmflag)then
     allocate(acheck(2*Nindip))
  else
     allocate(acheck(2*Ns))
  endif

  if(mpiID==0) print*,"Working on the real frequencies axis"

  ! calcolo la funzione di Green non interagente del sistema corrispondente 

  if(mpiID==0) print*,"Non-Interacting System"
  call get_gloc0_mpi()
  if(mpiID==0)call splot(trim(adjustl(trim(name_dir)))//"/LG0_realw.data",fg_0(1:Ns,1:L),wr(1:L))  


  !==============================================================
  !START DMFT LOOP SEQUENCE:
  !==============================================================
  iloop=0
  converged=.false.
  if(U==0.d0)converged=.true.   ! salto il loop DMFT se U=0
  call setup_initial_sc_sigma()
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !SOLVE G_II (GLOCAL) \FORALL FREQUENCIES: 
     call get_sc_gloc_mpi()        

     !SOLVE IMPURITY MODEL,\FORALL LATTICE SITES: DA CAMBIARE USANDO LE SIMMETRIE
     call solve_sc_impurity_mpi()  

     n_tot    = sum(nii)
     delta_tot= sum(dii)
     if(mpiID==0)write(*,"(A,F15.9)")"Number of Particles=",n_tot

     if (symmflag) then
        acheck(1:Nindip)=reshuffled(nii)
        acheck(Nindip+1:2*Nindip)=reshuffled(dii)
     else  
        acheck(1:Ns)=nii
        acheck(Ns+1:2*Ns)=dii
     endif
     converged=check_convergence_scalar(acheck,eps_error,Nsuccess,nloop,id=0,strict=.true.)

     if (densfixed) call search_mu_trap(converged) ! search_mu(converged)

     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     if(converged)call compute_gapmap_mpi  
     call print_sc_out(converged)
     call end_loop()
  enddo

  deallocate(fg,fg0,fg_0,sigma,sigma_tmp,nii_tmp,dii_tmp,gap_ii_tmp)

  if(mpiID==0) then 
     call system("mv -vf *.err "//trim(adjustl(trim(name_dir)))//"/")
     open(10,file="used.inputRDMFT.in")
     write(10,nml=disorder)
     close(10)
     call system("cp -vf used.*.in "//trim(adjustl(trim(name_dir)))//"/ ")
  endif


  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)




contains


  !******************************************************************
  !******************************************************************


  pure function trap_distance_square(is) result(dist)
    integer,intent(in) :: is
    integer            :: center
    real(8)            :: dist
    ! The center is in (0,0)
    dist=dble(irow(is))**2 + dble(icol(is))**2  ! questa ora non era cambiata
  end function trap_distance_square



  !******************************************************************
  !******************************************************************


  subroutine setup_initial_sc_sigma()
    logical :: check1,check2,check
    if(mpiID==0)then
       inquire(file="LSigma_realw.ipt",exist=check1)
       if(.not.check1)inquire(file="LSigma_realw.ipt.gz",exist=check1)
       inquire(file="LSelf_realw.ipt",exist=check2)
       if(.not.check2)inquire(file="LSelf_realw.ipt.gz",exist=check2)
       check=check1.AND.check2
       if(check)then
          call msg("Reading Self-energy from file:",lines=2)
          call sread("LSigma_realw.ipt",sigma(1,1:Ns,1:L),wr)
          call sread("LSelf_realw.ipt",sigma(2,1:Ns,1:L),wr)
       else
          call msg("Using Hartree-Fock-Bogoliubov self-energy:",lines=2)
          sigma(1,:,:)=zero ; sigma(2,:,:)=-deltasc
       endif
    endif
    call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  end subroutine setup_initial_sc_sigma



  !******************************************************************
  !******************************************************************

  subroutine get_sc_gloc_mpi()
    complex(8) :: Gloc(2*Ns,2*Ns),gf_tmp(2,Ns,L),zeta1,zeta2
    integer    :: i,is
    call msg("Get local GF: (ETA --> fort.999)",id=0)
    call start_timer
    fg=zero ; gf_tmp=zero
    do i=1+mpiID,L,mpiSIZE          !parallelizziamo sulle frequenze... TODO memory proximity
       Gloc=zero
       Gloc(1:Ns,1:Ns)          = -H0
       Gloc(Ns+1:2*Ns,Ns+1:2*Ns)=  H0
       do is=1,Ns
          zeta1=       cmplx(wr(i),eps,8)     -a0trap + xmu - sigma(1,is,i)      - etrap(is)
          zeta2=-conjg(cmplx(wr(L+1-i),eps,8) -a0trap + xmu - sigma(1,is,L+1-i)) + etrap(is)
          Gloc(is,is)      =  zeta1
          Gloc(Ns+is,Ns+is)=  zeta2
          Gloc(is,Ns+is)   = -sigma(2,is,i)
          Gloc(Ns+is,is)   = -sigma(2,is,L+1-i)   ! check
       enddo
       call matrix_inverse_sym(Gloc)!,2*Ns)
       forall(is=1:Ns)
          gf_tmp(1,is,i) = Gloc(is,is)
          gf_tmp(2,is,i) = Gloc(is,Ns+is)
       end forall
       call eta(i,L,999)
    enddo
    call stop_timer
    call MPI_REDUCE(gf_tmp,fg,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(fg,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_sc_gloc_mpi


  !******************************************************************
  !******************************************************************

  subroutine compute_gapmap_mpi   ! evaluate spectral gap for every lattice site
    integer    :: is
    gap_ii=0.d0; gap_ii_tmp=0.d0
    call msg("Computing Local Spectral Gap: (ETA --> fort.969)",id=0)
    call start_timer
    do is=1+mpiID,Ns,mpiSIZE         
       call get_gap(-dimag(fg(1,is,:))/pi,gap_ii_tmp(is),wr) 
       call eta(is,L,969)
    enddo
    call stop_timer
    call MPI_REDUCE(gap_ii_tmp,gap_ii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(gap_ii,Ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine compute_gapmap_mpi

  !******************************************************************
  !******************************************************************


  subroutine get_gloc0_mpi()   ! in the noninteracting case the spins are independent
    complex(8) :: Gloc0(Ns,Ns),gf0_tmp(Ns,L),zeta   ! tutte locali
    integer    :: i,is
    call msg("Get local noninteracting GF: (ETA --> fort.999)",id=0)
    call start_timer
    fg_0=zero ; gf0_tmp=zero      ! fg_0 e' l'unica in uscita
    do i=1+mpiID,L,mpiSIZE          
       Gloc0=zero
       Gloc0(1:Ns,1:Ns)          = -H0        
       do is=1,Ns
          zeta=       cmplx(wr(i),eps,8)     - a0trap + xmu    - etrap(is)
          Gloc0(is,is)      =  zeta
       enddo
       call matrix_inverse_sym(Gloc0)!,Ns)
       forall(is=1:Ns)
          gf0_tmp(is,i) = Gloc0(is,is)
       end forall
       call eta(i,L,999)
    enddo
    call stop_timer
    call MPI_REDUCE(gf0_tmp,fg_0,Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(fg_0,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine get_gloc0_mpi



  !******************************************************************
  !******************************************************************


  !Usiamo le simmetrie per risolvere solo siti non equivalenti
  subroutine solve_sc_impurity_mpi()
    integer    :: is,i
    call msg("Solve impurity:")
    sigma_tmp = zero      ! dummy variable for reduce purposes.. 
    nii_tmp   = 0.d0      ! density profile to be calculated on the fly 
    dii_tmp   = 0.d0      ! delta profile to be calculated on the fly 
    !
    call start_timer
    if(.not.symmflag)then       ! tutti i siti sono inequivalenti
       do is=1+mpiID,Ns,mpiSIZE
          call solve_per_site(is)
          call eta(is,Ns,998)
       enddo
    else                        ! uso solo quelli indipendenti e li simmetrizzo
       do is=1+mpiID,Nindip,mpiSIZE
          call solve_per_site(indipsites(is))
          call eta(is,Nindip,998)
       enddo
       !ho calcolato le self-energie per i siti non equivalenti
       do i=1,L
          call symmetrize(sigma_tmp(1,:,i))
          call symmetrize(sigma_tmp(2,:,i))
       enddo
       !ho implementato le simmetrie del cubo sulle self-energie 
       !adesso implementa le simmetrie del cubo sul profilo di densita' e delta
       call symmetrize(nii_tmp)
       call symmetrize(dii_tmp)
    endif
    call stop_timer
    !
    sigma=zero ; nii=zero ; dii=zero
    call MPI_REDUCE(sigma_tmp,sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    !
    call MPI_REDUCE(nii_tmp,nii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(nii,Ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    !
    call MPI_REDUCE(dii_tmp,dii,Ns,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(dii,Ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    !nii=nii_ ; dii=dii_
  end subroutine solve_sc_impurity_mpi



  !******************************************************************
  !******************************************************************



  subroutine solve_per_site(is)
    integer                                      :: is,i
    complex(8)                                   :: det
    complex(8),dimension(:,:,:),allocatable,save :: sold
    complex(8),dimension(2,L)                    :: calG
    real(8),dimension(2,0:L)                     :: fgt,fg0t
    real(8)                                      :: n,n0,delta,delta0

    if(.not.allocated(sold))allocate(sold(2,Ns,L))
    sold(:,is,:)  =  sigma(:,is,:)
    !
    n    = -sum(dimag(fg(1,is,:))*fermi(wr,beta))*fmesh/pi ! densita' per spin
    delta= -(u*sum(dimag(fg(2,is,:))*fermi(wr,beta))*fmesh/pi) 
    !
    nii_tmp(is)=2.d0*n; dii_tmp(is)=delta
    !
    fg0=zero ; calG=zero
    do i=1,L
       det      = fg(1,is,i)*conjg(fg(1,is,L+1-i)) + conjg(fg(2,is,L+1-i))*fg(2,is,i)
       calG(1,i)= conjg(fg(1,is,L+1-i))/det  + sigma(1,is,i)      + u*(n-0.5d0)
       calG(2,i)= fg(2,is,i)/det       + conjg(sigma(2,is,L+1-i)) + delta
    end do
    do i=1,L
       det     =  calG(1,i)*conjg(calG(1,L+1-i)) + conjg(calG(2,L+1-i))*calG(2,i)
       fg0(1,i)=  conjg(calG(1,L+1-i))/det
       fg0(2,i)=  conjg(calG(2,L+1-i))/det
    end do
    n0    = -sum(dimag(fg0(1,:))*fermi(wr,beta))*fmesh/pi!/sum(dimag(fg(1,is,:)))
    delta0= -u*sum(dimag(fg0(2,:))*fermi(wr,beta))*fmesh/pi
    !    write(750,"(2I4,4(f16.12))",advance="yes")mpiID,is,n,n0,delta,delta0
    !
    sigma_tmp(:,is,:) =  solve_mpt_sc_sopt(fg0,wr,n,n0,delta,delta0,L)
    sigma_tmp(:,is,:) =  weight*sigma_tmp(:,is,:) + (1.d0-weight)*sold(:,is,:)
    !
  end subroutine solve_per_site



  !******************************************************************
  !******************************************************************



  subroutine print_sc_out(converged)
    integer                  :: is,i,row,col,occupied,center,border,corner
    real(8)                  :: n_center,n_av,n_corner,n_border,n_min
    real(8)                  :: delta_center,delta_max,delta_av
    real(8)                  :: C         ! CDW order parameter 
    real(8)                  :: e_corner
    real(8),parameter        :: density_threshold=0.01d0
    complex(8)               :: afg(2,L),asig(2,L)
    real(8)                  :: cdw(Ns)
    logical                  :: converged
    character(len=4)         :: loop
    real(8),dimension(2,0:L) :: fgt
    real(8),allocatable      :: grid_x(:),grid_y(:)
    real(8),allocatable      :: dij(:,:),nij(:,:),eij(:,:)
    real(8),allocatable      :: gap_ij(:,:)
    if(mpiID==0)then

       !BUILD A GRID FOR  LATTICE PLOTS:                 
       if (.not.allocated(grid_x)) then 
          allocate(grid_x(-Nside/2:Nside/2),grid_y(-Nside/2:Nside/2))
          allocate(nij(-Nside/2:Nside/2,-Nside/2:Nside/2))
          allocate(dij(-Nside/2:Nside/2,-Nside/2:Nside/2))
          allocate(eij(-Nside/2:Nside/2,-Nside/2:Nside/2))
          allocate(gap_ij(-Nside/2:Nside/2,-Nside/2:Nside/2))
          do row=-Nside/2,Nside/2
             grid_x(row)=dble(row)
             grid_y(row)=dble(row)
          enddo
       endif

       !Evaluate the CDW order-parameter and threshold for occupation
       occupied=0 
       C=0.d0
       do is=1,Ns
          cdw(is) = (nii(is)-1.d0)*((-1.d0)**abs(irow(is)+icol(is))) 
          C=C+cdw(is)
          if (nii(is).ge.density_threshold) occupied=occupied+1
       enddo

       !Special points
       corner=ij2site(Nside/2,Nside/2)
       center=ij2site(0,0)
       border=ij2site(0,Nside/2)

       !Average and total observables:
       n_av     = n_tot/dble(occupied)
       delta_av = delta_tot/dble(occupied)
       n_corner = nii(corner)
       n_border = nii(border)
       n_min    = minval(nii)
       e_corner = a0trap+etrap(corner)
       call splot(trim(adjustl(trim(name_dir)))//"/ntotVSiloop.ipt",iloop,n_tot,append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/davVSiloop.ipt",iloop,delta_av,append=TT)
       call system("rm -fv *site.ipt *.ipt.gz")   ! se append=false non serve ...

       !Print some information to user:
       print*,"========================================"
       print*,"Average density =",n_av
       print*,"Delta_av =",delta_av
       print*,"Minimum density =",n_min
       print*,"Residual density at the border=",n_border 
       if (n_border.gt.density_threshold) print*,"WARNING: MAYBE TOUCHING THE TRAP BOUNDARIES"
       print*,"========================================"
       print*,"Residual density at the corner=",n_corner
       print*,"Trap energy at the corner =",e_corner
       if (n_corner.gt.n_min+density_threshold) print*,"ACHTUNG: NON-MONOTONIC DENSITY PROFILE"
       print*,"========================================"       



       !Evaluate 3D distribution of density and order-parameter:
       do row=-Nside/2,Nside/2
          do col=-Nside/2,Nside/2
             i            = ij2site(row,col)
             nij(row,col) = nii(i)
             dij(row,col) = dii(i)
          enddo
       enddo
       call splot("3d_nVSij.ipt",grid_x,grid_y,nij)
       call splot("3d_deltaVSij.ipt",grid_x,grid_y,dij)


       !STORE GREEN's FUNCTIONS AND SELF-ENERGY IN COMPACT FORM TO SAVE SPACE
       call splot(trim(adjustl(trim(name_dir)))//"/LSigma_realw.ipt",sigma(1,1:Ns,1:L),wr(1:L))
       call splot(trim(adjustl(trim(name_dir)))//"/LSelf_realw.ipt",sigma(2,1:Ns,1:L),wr(1:L))
       call splot(trim(adjustl(trim(name_dir)))//"/LG_realw.ipt",fg(1,1:Ns,1:L),wr(1:L))
       call splot(trim(adjustl(trim(name_dir)))//"/LF_realw.ipt",fg(2,1:Ns,1:L),wr(1:L))


       !plotting selected green-function and self-energies for quick
       !data processing and debugging
       call splot(trim(adjustl(trim(name_dir)))//"/Gloc_realw_center.ipt",wr,fg(1,center,1:L),append=FF)
       call splot(trim(adjustl(trim(name_dir)))//"/Floc_realw_center.ipt",wr,fg(2,center,1:L),append=FF)
       call splot(trim(adjustl(trim(name_dir)))//"/Sigma_realw_center.ipt",wr,sigma(1,center,1:L),append=TT) 
       call splot(trim(adjustl(trim(name_dir)))//"/Self_realw_center.ipt",wr,sigma(2,center,1:L),append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/Gloc_realw_border.ipt",wr,fg(1,border,1:L),append=FF)
       call splot(trim(adjustl(trim(name_dir)))//"/Floc_realw_border.ipt",wr,fg(2,border,1:L),append=FF)
       call splot(trim(adjustl(trim(name_dir)))//"/Sigma_realw_border.ipt",wr,sigma(1,border,1:L),append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/Self_realw_border.ipt",wr,sigma(2,border,1:L),append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/Gloc_realw_corner.ipt",wr,fg(1,corner,1:L),append=FF)
       call splot(trim(adjustl(trim(name_dir)))//"/Floc_realw_corner.ipt",wr,fg(2,corner,1:L),append=FF)
       call splot(trim(adjustl(trim(name_dir)))//"/Sigma_realw_corner.ipt",wr,sigma(1,corner,1:L),append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/Self_realw_corner.ipt",wr,sigma(2,corner,1:L),append=TT)




       !IF CONVERGENCE IS ACHIEVED PRINT SHORT SUMMARY OF INFO and local observables:
       if(converged) then
          call msg("*********** SUMMARY *************")
          n_center= nii(center); delta_center=dii(center); delta_max=maxval(dii)
          write(*,"(A12,F18.12)")"trap_occupancy",dble(occupied)/dble(Ns)
          write(*,"(A12)")"Density-----------------"
          write(*,"(A12,F18.12)")"center  =",n_center
          write(*,"(A12,F18.12)")"Total   =",n_tot
          write(*,"(A12,F18.12)")"average =",n_av
          write(*,"(A12)")"Delta-------------------"
          write(*,"(A12,F18.12)")"center  =",delta_center
          write(*,"(A12,F18.12)")"average =",delta_av
          write(*,"(A12,F18.12)")"max     =",delta_max
          print*,"CDW ORDER PARAMETER",C
          print*,"**********************     END  SUMMARY        *************************"





          !Get the trap shape and print converged 2d/3d-plots
          do row=-Nside/2,Nside/2
             do col=-Nside/2,Nside/2
                i            = ij2site(row,col)
                eij(row,col) = etrap(i)
                gap_ij(row,col) = gap_ii(i)
             enddo
          enddo

          call splot(trim(adjustl(trim(name_dir)))//"/2d_nVSij.data",grid_y,nij(0,:))
          call splot(trim(adjustl(trim(name_dir)))//"/2d_deltaVSij.data",grid_y,dij(0,:))
          call splot(trim(adjustl(trim(name_dir)))//"/2d_etrapVSij.data",grid_y,eij(0,:))
          call splot(trim(adjustl(trim(name_dir)))//"/2d_gapVSij.data",grid_y,gap_ij(0,:))

          call splot(trim(adjustl(trim(name_dir)))//"/3d_nVSij.data",grid_x,grid_y,nij)
          call splot(trim(adjustl(trim(name_dir)))//"/3d_deltaVSij.data",grid_x,grid_y,dij)
          call splot(trim(adjustl(trim(name_dir)))//"/3d_etrapVSij.data",grid_x,grid_y,eij)
          call splot(trim(adjustl(trim(name_dir)))//"/3d_gapVSij.data",grid_x,grid_y,gap_ij)

          call system("mv -vf used.*.in "//trim(adjustl(trim(name_dir)))//"/ ")
       endif

    end if
    return
  end subroutine print_sc_out


  !******************************************************************
  !******************************************************************


  subroutine search_mu(convergence)
    integer, save         ::nindex
    integer               ::nindex1
    real(8),save          ::muold,N_old
    real(8)               :: ndelta1
    logical,intent(inout) ::convergence
    if(mpiID==0)then
       nindex1=nindex  !! save old value of the increment sign
       ndelta1=ndelta
       !      muold=a0trap    !! implement to extrapulate a sound value of chitrap
       !      N_old=N_old     !! from the same procedure to fix the density
       if((n_tot >= n_wanted+n_tol))then
          nindex=-1
       elseif(n_tot <= n_wanted-n_tol)then
          nindex=1
       else
          nindex=0
       endif

       !       if (nindex1+nindex==0) then      
       if((nindex1+nindex==0).and.(nindex.ne.0))then       !avoid loop back and back [CHANGE, avoid reducing the step for consecutive convergencies
          ndelta=ndelta1/2.d0
       endif
       a0trap=a0trap-chitrap*real(nindex,8)*ndelta  ! a0trap=-mu_tot

       ! in this way chitrap is the inverse trap compressibility chitrap=dmu/d(N_tot)

       print*,"=======================  DENSITY-LOOP   ========================="
       write(*,"(A,f15.12,A,f15.12)")"A0TRAP=",a0trap," step =",ndelta
       write(*,"(A,f15.12,A,f15.12)")"density error=",abs(n_tot-n_wanted)," vs",n_tol
       if(abs(n_tot-n_wanted) > n_tol) then
          print*,"********* density loop not yet converged ***********" 
          if (iloop < nloop) then
             convergence=.false.
          else
             print*,"       FORCED DENSITY LOOP EXIT !!        "
             print*,"CONSIDER INCREASING NLOOP OR CHANGE CHITRAP"
             convergence=.true.
          endif
       else
          print*,"*********     density-loop CONVERGED      ***********"
       endif
       print*,"================================================================="

       call splot(trim(adjustl(trim(name_dir)))//"/a0trapVSiter.ipt",iloop,a0trap,abs(n_tot-n_wanted),append=.true.)
    endif
    call MPI_BCAST(a0trap,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine search_mu

  !******************************************************************
  !******************************************************************



  !****************************************************************************************!
  !  Improved search mu routine which adapt both density and DMFT threshold in order       !
  !  achieve BOTH convergencies simultaneously                                             !
  !  This version is slightly modified with respect to search_mu_2 used in single site code!
  !  version 3 will be portable and generic to be included in a separate modulus           !                                                 ! 
  !****************************************************************************************!
  subroutine search_mu_trap(convergence) 
    integer, save         :: nindex
    integer               :: nindex_old
    integer,save          :: nfactor,dmftfactor  ! factors to adapt threshold for the DMFT and density loop 
    real(8),save          :: n_threshold,dmft_threshold,step_init   
    integer               :: idritto             ! avoid too many steps in the same direction [NOT USED]
    real(8),save          :: chi,chi_phys        ! local and estimate of the physical compressibility
    real(8),save          :: xmu_old(2),n_old
    logical,intent(inout) :: convergence

    if(mpiID==0)then

       if (iloop==1) then 
          print*,"First Time Inside Search_Mu"
          dmftfactor=50  ! initial threshold for the DMFT loop     = dmftfactor*eps_error
          nfactor=100     ! initial threshold for the density loop  = nfactor*nerror
          idritto=0      ! to avoid too many steps in one direction [NOT USED SO FAR]
          !      inverto il ruolo delle variabili cosi' uso le stesse al di fuori della routine
          dmft_threshold=eps_error ; n_threshold=n_tol   ! saved
          step_init=ndelta                               ! saved
          eps_error=max(real(dmftfactor,8)*dmft_threshold,dmft_threshold) 
          n_tol=max(real(nfactor,8)*n_threshold,n_threshold) 
          print*,"Initial DMFT threshold",eps_error
          print*,"Initial density threshold",n_tol
          chi=chitrap     ! initial guess for chi
       endif


       print*,"                                                                        "
       print*,"========================     Density loop    ==========================="
       !    NB se xmu non varia perche il density loop e' convergiuto la chi e' zero... !!!!!   
       !    per cui devo metterci questo if
       if (xmu_old(1).ne.(-a0trap)) then
          xmu_old(2)=xmu_old(1)  ![CHANGE WITH VECTORS]
          xmu_old(1)=-a0trap  ! save current and previous value of the chemical potential to estimate chi 
          chi=(xmu_old(1)-xmu_old(2))/(n_tot-n_old)   !! stimo la compressibilita' cosi' funziona solo dopo la seconda convergenza
          print*,"Actual Trap Compressibility =",chi
          !         TO BE IMPLEMENTED
          !         NOT CHANGE MU WHEN chi < 0 (pensa alla molla forzata) 
       endif

       ! updato le soglie per n e DMFT solo se sono convergiuti tutti e due contemporaneamente
       if((convergence).and.(abs(n_tot-n_wanted) < n_tol)) then 
          print*,"PREVIOUS DMFT loop converged withing the density threshold"
          print*,"Decreasing Density and DMFT Thresholds"
          nfactor=nfactor/2 ; dmftfactor=dmftfactor/2
          eps_error=max(real(dmftfactor,8)*dmft_threshold,dmft_threshold)   ! ora cambio eps_error che viene letta da check convergence mentre
          n_tol=max(nfactor*n_threshold,n_threshold)                        ! ho salvato la soglia finale in dmft_threshold, stessa cosa su n_tol per semplicita'
          ndelta=ndelta/2                                                   ! TO BE IMPROVED questa sembra funzionare meglio  

          print*,"NEW density threshold",n_tol
          print*,"NEW   DMFT  threshold",eps_error
          print*,"NEW step =",ndelta
          if (chi > 0) then !.and.(chi<=10.d0)) then
             chi_phys=chi 
             print*,"UPDATED PHYSICAL CHI =",chi_phys                ! idealmente dovrei usarlo per fare una wise choice dell ndelta iniziale.., l'update e' un po' complicato..
             !             ndelta=min(ndelta/2,abs(n_tot-n_wanted)*chi_phys/5.d0)  ! uso chi per un guess sullo step, usando sempre lo step piu' conservativa 
             !          else
             !             print*,"Unphysical compressibility"
          endif
          !      real(ndelta/2.d0,8)
          !      uso solo i valori dei sistemi a convergenza per stimare chi 
          !      in modo da avere delle compressibilita' fisiche
          !      non funziona... :-( 
       endif


       nindex_old=nindex     ! save old value of the error sign for bracketing purposes
       n_old=n_tot               ! save old value of the density to estimate chi
       !       if (chi > 0) then ! not updating mu if chi <0   it doesn't work...
       if((n_tot >= n_wanted + n_tol))then    
          nindex=-1
          if (idritto.le.0) idritto=idritto-1        
       elseif(n_tot <= n_wanted - n_tol)then
          nindex=1
          if (idritto.le.0) idritto=idritto+1
       else
          print*,"Density loop converged within the actual threshold",n_tol
          nindex=0   ! sono dentro la threshold attuale e non mi muovo
          idritto=0
       endif
       ! else 
       !    print*,"Negative compressibility, mu not updated"
       !    nindex=0   ;   idritto=0
       ! endif

       !   decido come fare lo step 
       !   voglio evitare sia step troppo brevi 
       !   che troppo lunghi per una data threshold
       !   inoltre quando risetto la threshold devo risettare gli step [TODO] maybe not needed
       if((nindex_old+nindex==0).and.(nindex.ne.0))then       !  avoid loop forth and back by  decreasing the step. 
          ndelta=real(ndelta/2.d0,8)        ! The second condition avoid  further decreasing the step for consecutive convergence of the density loop
          idritto=0 
          print*,"decreasing step to avoid loop back and forth, actual step =",ndelta
          !PER ADESSO LO LASCIO COSI' SENZA IMPLEMENTARE increasing steps
          !     elseif ((abs(idritto)>=10).and.(converged)) then           ! conservative version [useless]
          !     elseif (abs(idritto)>=10) then           ! too many steps in the same direction -> increase steps if the DMFT loop is converged, decrease otherwise! 
          !        idritto=0 
          ! !       if (converged) then
          !           ndelta=real(ndelta*2.d0,8)
          !           print*,"Increasing step, actual step =",ndelta
          ! !       else
          ! !          ndelta=real(ndelta/2.d0,8)
          ! !          print*,"Decreasing step, actual step =",ndelta
          ! !       endif       

       endif


       !============================= update chemical potential========================
       ! xmu=xmu+chi*real(nindex,8)*ndelta                   
       a0trap=a0trap-real(nindex,8)*ndelta                   

       write(*,"(A,2f15.8,A,I5,A,I3,f15.8)")"mu,n=",xmu_old(1),n_tot,"/",n_wanted,"| ",nindex,ndelta
       print*,"========================================================================"
       if((abs(n_tot-n_wanted) > n_threshold).OR.(eps_error > dmft_threshold)) then 
          if (iloop<=nloop) convergence=.false.  ! avoid too many loops 
       elseif(convergence) then
          ! save variables in used.inputRDMFT.in for future runs 
          ! dovrei essere convergiuto ! NB devo convergere 3 volte problema            
          eps_error=dmft_threshold ; n_tol=n_threshold   
          ndelta=step_init ; chitrap=chi_phys            
       endif
       call splot(trim(adjustl(trim(name_dir)))//"/a0trapVSiter.ipt",iloop,xmu_old(1),abs(n_tot-n_wanted),chi,append=.true.)
    endif
    call MPI_BCAST(a0trap,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine search_mu_trap








  !******************************************************************
  !******************************************************************

  subroutine get_gap_old(g_imag,gap_)
    real(8)                  :: min,gap_,threshold
    real(8)                  :: g_imag(:)
    integer                  :: i
    threshold=0.1d0   ! the gap value is not computed below this threshold 
    ! for the abs(g_imag) value. 
    min=0.d0
    gap_=0.d0
    do i=L/2-1,1,-1                       ! ho un problema di soglia se la gap non e' perfettamente simmetrica
       if(g_imag(i).lt.min) then                      ! potrei avere un massimo spostato
          min=g_imag(i)
       elseif(abs(g_imag(i)).gt.threshold)then    ! metto una soglia sui valori della DOS ai picchi di coerenza 
          gap_=-wr(i+1)                            ! per evitarla
          return
          !goto 100
       endif
    enddo
    !100 continue
  end subroutine get_gap_old


  !**********************************************************************!
  ! finds the closest maximum to the axis origin,
  ! which given the shape of the DOS for a superconducting 
  ! system has the energy gap/2 as a x-coordinate 
  ! Date : 24-09-2012
  !**********************************************************************!
  subroutine get_gap(dos_,gap_,wr_)
    real(8)                  :: gap_plus,gap_minus,gap_
    real(8)                  :: dos_(:),wr_(:)
    integer, dimension(1)    :: i_plus,i_minus   
    ! trick needed since maxloc always results in a vector variable
    i_plus  =maxloc(dos_(L/2+1:L))
    i_minus =maxloc(dos_(1:L/2))
    gap_plus  = wr_(i_plus(1)+L/2)
    gap_minus =-wr_(i_minus(1))
    print*,"gap_+ =", gap_plus
    print*,"gap_ =" , gap_minus
    gap_=min(gap_plus,gap_minus)
  end subroutine get_gap

  !******************************************************************
  !******************************************************************


end program ahm_real_trap
