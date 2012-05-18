!########################################################
!     PURPOSE  :solve the attractive (A) Trapped (T) Hubbard
! model (HM) using Modified Perturbation THeory (MPT), w/ DMFT
!     AUTHORS  : A.Amaricci, A.Priviter (CNR-IOM)
!########################################################

program ahm_matsubara_trap
  USE RDMFT_VARS_GLOBAL
  implicit none
  real(8)                                 :: r
  real(8)                                 :: n_tot,delta_tot
  integer                                 :: is,esp,lm
  logical                                 :: converged
  complex(8),allocatable,dimension(:,:)   :: fg0
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma,sigma_tmp
  real(8),allocatable,dimension(:)        :: nii_tmp,dii_tmp


  !GLOBAL INITIALIZATION:
  !===============================================================
  include "init_global_trap.f90" 


  !ALLOCATE WORKING ARRAYS
  !===============================================================
  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)

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

     !SOLVE G_II (GLOCAL) \FORALL FREQUENCIES: 
     call get_sc_gloc_mpi()        

     !SOLVE IMPURITY MODEL,\FORALL LATTICE SITES: DA CAMBIARE USANDO LE SIMMETRIE
     call solve_sc_impurity_mpi()  

     n_tot    = sum(nii)
     delta_tot= sum(dii)
     if(mpiID==0)write(*,"(A,F15.9)")"Number of Particles=",n_tot

     if (symmflag) then
        converged=check_convergence(reshuffled(sigma(1,:,:)+sigma(2,:,:)),eps_error,Nsuccess,nloop,id=0,tight=.true.)
     else  
        converged=check_convergence(sigma(1,:,:)+sigma(2,:,:),eps_error,Nsuccess,nloop,id=0,tight=.true.) 
     endif
     if (densfixed) call search_mu(converged)

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
       inquire(file="LSigma.ipt",exist=check1)
       if(.not.check1)inquire(file="LSigma.ipt.gz",exist=check1)
       inquire(file="LSelf.ipt",exist=check2)
       if(.not.check2)inquire(file="LSelf.ipt.gz",exist=check2)
       check=check1.AND.check2
       if(check)then
          call msg("Reading Self-energy from file:",lines=2)
          call sread("LSigma.ipt",sigma(1,1:Ns,1:L))
          call sread("LSelf.ipt",sigma(2,1:Ns,1:L))
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
    complex(8) :: Gloc(2*Ns,2*Ns),gf_tmp(2,Ns,L)
    integer    :: i,is
    call msg("Get local GF: (ETA --> fort.999)",id=0)
    call start_timer
    fg=zero ; gf_tmp=zero
    do i=1+mpiID,L,mpiSIZE          ! parallelizziamo sulle frequenze... TODO memory proximity
       Gloc=zero
       Gloc(1:Ns,1:Ns)          = -H0
       Gloc(Ns+1:2*Ns,Ns+1:2*Ns)=  H0
       do is=1,Ns
          Gloc(is,is)      =  xi*wm(i)-a0trap-sigma(1,is,i)-etrap(is)
          Gloc(Ns+is,Ns+is)= -conjg(Gloc(is,is))
          Gloc(is,Ns+is)   = -sigma(2,is,i)
          Gloc(Ns+is,is)   = -sigma(2,is,i)!sigma(2,is,L+1-i)!this should be a simmetry in Matsubara!
       enddo
       call mat_inversion_sym(Gloc,2*Ns)
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
    complex(8)                                   :: det(L)
    complex(8),dimension(:,:,:),allocatable,save :: sold
    complex(8),dimension(2,L)                    :: calG
    real(8),dimension(2,0:L)                     :: fgt,fg0t
    real(8)                                      :: n,n0,delta,delta0

    if(.not.allocated(sold))allocate(sold(2,Ns,L))
    sold(:,is,:)  =  sigma(:,is,:)
    call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
    call fftgf_iw2tau(fg(2,is,:),fgt(2,0:L),beta,notail=.true.)
    n    = -fgt(1,L) ; delta= -u*fgt(2,L)
    !
    nii_tmp(is)=2.d0*n; dii_tmp(is)=delta
    !
    fg0=zero ; calG=zero
    det       = abs(fg(1,is,:))**2    + fg(2,is,:)**2
    fg0(1,:)  = conjg(fg(1,is,:))/det + sigma(1,is,:) + U*(n-0.5d0)
    fg0(2,:)  = fg(2,is,:)/det        + sigma(2,is,:) + delta
    det       =  abs(fg0(1,:))**2     + fg0(2,:)**2
    calG(1,:) =  conjg(fg0(1,:))/det
    calG(2,:) =  fg0(2,:)/det
    call fftgf_iw2tau(calG(1,:),fg0t(1,:),beta)
    call fftgf_iw2tau(calG(2,:),fg0t(2,:),beta,notail=.true.)
    n0=-fg0t(1,L) ; delta0= -u*fg0t(2,L)
    write(750,"(2I4,4(f16.12))",advance="yes")mpiID,is,n,n0,delta,delta0
    !
    sigma_tmp(:,is,:) =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)
    sigma_tmp(:,is,:) =  weight*sigma_tmp(:,is,:) + (1.d0-weight)*sold(:,is,:)
    !
  end subroutine solve_per_site



  !******************************************************************
  !******************************************************************



  subroutine print_sc_out(converged)
    integer                  :: is,i,row,col,occupied,center,border,corner
    real(8)                  :: n_center,n_av,n_corner,n_border,n_min
    real(8)                  :: delta_center,delta_max,delta_av
    real(8)                  :: e_corner
    real(8),parameter        :: density_threshold=0.01d0
    complex(8)               :: afg(2,L),asig(2,L)
    real(8)                  :: cdw(Ns)
    logical                  :: converged
    character(len=4)         :: loop
    real(8),dimension(2,0:L) :: fgt
    real(8),allocatable                     :: grid_x(:),grid_y(:)
    real(8),allocatable                     :: dij(:,:),nij(:,:),eij(:,:),cdwij(:,:)

    if(mpiID==0)then

       occupied=0 

       do is=1,Ns
          cdw(is) = (nii(is)-1.d0)*((-1.d0)**abs(irow(is)+icol(is))) 
          if (nii(is).ge.density_threshold) occupied=occupied+1
       enddo

       corner=ij2site(Nside/2,Nside/2)
       center=ij2site(0,0)
       border=ij2site(0,Nside/2)

       n_av  = n_tot/dble(occupied); delta_av=delta_tot/dble(occupied)
       n_corner=nii(corner); n_border=nii(border); n_min=minval(nii)
       e_corner=a0trap+etrap(corner)

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

       call splot(trim(adjustl(trim(name_dir)))//"/ntotVSiloop.ipt",iloop,n_tot,append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/davVSiloop.ipt",iloop,delta_av,append=TT)

       call system("rm -fv *site.ipt *.ipt.gz")   ! se append=false non serve ...
       !       call system("rm -fv *site.ipt *loc*.ipt *.ipt.gz")   ! se append=false non serve ... !

       !-------------------------------------------------------------------------------------------       
       !      to be removed.. after debugging store only compact forms of
       !      Self-energy and Green functions to save space 

       do i=1,Ns
          call splot("Gloc_iw_site.ipt",wm,fg(1,i,1:L),append=TT)
          call splot("Floc_iw_site.ipt",wm,fg(2,i,1:L),append=TT)
          call splot("Sigma_iw_site.ipt",wm,sigma(1,i,1:L),append=TT)
          call splot("Self_iw_site.ipt",wm,sigma(2,i,1:L),append=TT)
       enddo
       !******************************************************************************

       !   plotting selected green-function and self-energies for quick
       !   data processing and debugging

       call splot("Gloc_iw_center.ipt",wm,fg(1,center,1:L),append=FF)
       call splot("Floc_iw_center.ipt",wm,fg(2,center,1:L),append=FF)
       call splot("Sigma_iw_center.ipt",wm,sigma(1,center,1:L),append=TT) ! controllo come evolvono le selfenergie con le iterazioni 
       call splot("Self_iw_center.ipt",wm,sigma(2,center,1:L),append=TT)


       call splot("Gloc_iw_border.ipt",wm,fg(1,border,1:L),append=FF)
       call splot("Floc_iw_border.ipt",wm,fg(2,border,1:L),append=FF)
       call splot("Sigma_iw_border.ipt",wm,sigma(1,border,1:L),append=TT)
       call splot("Self_iw_border.ipt",wm,sigma(2,border,1:L),append=TT)


       call splot("Gloc_iw_corner.ipt",wm,fg(1,corner,1:L),append=FF)
       call splot("Floc_iw_corner.ipt",wm,fg(2,corner,1:L),append=FF)
       call splot("Sigma_iw_corner.ipt",wm,sigma(1,corner,1:L),append=TT)
       call splot("Self_iw_corner.ipt",wm,sigma(2,corner,1:L),append=TT)


       !          STORE GREEN's FUNCTIONS AND SELF-ENERGY IN COMPACT FORM TO SAVE SPACE
       !          I do it at every interaction (even though is much slower) 
       !          so that in case of chrashes we can automatically restart with 
       !          a meaningful self-energy (just unzip LSigma.ipt.gz)
       !
       call splot("LSigma.ipt",sigma(1,1:Ns,1:L))
       call splot("LSelf.ipt",sigma(2,1:Ns,1:L))
       call splot("LG.ipt",fg(1,1:Ns,1:L))
       call splot("LF.ipt",fg(2,1:Ns,1:L))

       !--------------------------------------------------------------------------------

       if(converged) then

          print*,"**********************       SUMMARY        *************************"

          n_center= nii(center); delta_center=dii(center); delta_max=maxval(dii)

          print*,"trap_occupancy",dble(occupied)/dble(Ns)
          print*,"Density-----------------"
          print*,"center =",n_center
          print*,"Total =",n_tot
          print*,"average =",n_av
          print*,"Delta-------------------"
          print*,"center =",delta_center
          print*,"average =",delta_av
          print*,"max =",delta_max

          print*,"**********************     END  SUMMARY        *************************"

          !      redundant..

          !          call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.ipt",nii)
          !          call splot(trim(adjustl(trim(name_dir)))//"/deltaVSisite.ipt",dii)
          !          call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.ipt",erandom)


          !       do i=1,Ns
          !          call splot(trim(adjustl(trim(name_dir)))//"/Gloc_iw_site.ipt",wm,fg(1,i,1:L),append=TT)
          !          call splot(trim(adjustl(trim(name_dir)))//"/Floc_iw_site.ipt",wm,fg(2,i,1:L),append=TT)
          !          call splot(trim(adjustl(trim(name_dir)))//"/Sigma_iw_site.ipt",wm,sigma(1,i,1:L),append=TT)
          !          call splot(trim(adjustl(trim(name_dir)))//"/Self_iw_site.ipt",wm,sigma(2,i,1:L),append=TT)
          !       enddo

          !          call splot(trim(adjustl(trim(name_dir)))//"/LSigma.ipt",sigma(1,1:Ns,1:L))
          !          call splot(trim(adjustl(trim(name_dir)))//"/LSelf.ipt",sigma(2,1:Ns,1:L))
          !          call splot(trim(adjustl(trim(name_dir)))//"/LG.ipt",fg(1,1:Ns,1:L))
          !          call splot(trim(adjustl(trim(name_dir)))//"/LF.ipt",fg(2,1:Ns,1:L))


          ! afg(:,:)   =sum(fg(1:2,1:Ns,1:L),dim=2)/dble(Ns)
          ! asig(:,:)  =sum(sigma(1:2,1:Ns,1:L),dim=2)/dble(Ns)
          ! call splot(trim(adjustl(trim(name_dir)))//"/Gave_iw.ipt",wm,afg(1,1:L))
          ! call splot(trim(adjustl(trim(name_dir)))//"/Sigmaave_iw.ipt",wm,asig(1,1:L))


          !BUILD A GRID FOR  LATTICE PLOTS:

          allocate(grid_x(-Nside/2:Nside/2),grid_y(-Nside/2:Nside/2))
          allocate(nij(-Nside/2:Nside/2,-Nside/2:Nside/2))
          allocate(dij(-Nside/2:Nside/2,-Nside/2:Nside/2))
          allocate(eij(-Nside/2:Nside/2,-Nside/2:Nside/2))
          allocate(cdwij(-Nside/2:Nside/2,-Nside/2:Nside/2))

          do row=-Nside/2,Nside/2
             grid_x(row)=dble(row)
             grid_y(row)=dble(row)
          enddo

          do row=-Nside/2,Nside/2
             do col=-Nside/2,Nside/2
                i            = ij2site(row,col)
                nij(row,col) = nii(i)
                dij(row,col) = dii(i)
                eij(row,col) = etrap(i)
                cdwij(row,col)=cdw(i)
             enddo
          enddo
          call splot("2d_nVSij.data",grid_y,nij(0,:))
          call splot("2d_deltaVSij.data",grid_y,dij(0,:))
          call splot("2d_etrapVSij.data",grid_y,eij(0,:))
          call splot("2d_cdwVSij.data",grid_y,cdwij(0,:))
          call splot("3d_nVSij.data",grid_x,grid_y,nij)
          call splot("3d_deltaVSij.data",grid_x,grid_y,dij)
          call splot("3d_etrapVSij.data",grid_x,grid_y,eij)
          call splot("3d_cdwVSij.data",grid_x,grid_y,cdwij)
          call system("mv -vf 3d_*.data plot_3d* "//trim(adjustl(trim(name_dir)))//"/ ")
          call system("mv -vf *.ipt used.*.in "//trim(adjustl(trim(name_dir)))//"/ ")
          call system("cp -vf *.ipt.gz "//trim(adjustl(trim(name_dir)))//"/ ")
          call system("gunzip -v LSigma.ipt.gz LSelf.ipt.gz")
          call system("rm -vf *ipt.gz")
       endif
    end if
    return
  end subroutine print_sc_out




  !******************************************************************
  !******************************************************************



  function reshuffled(m_in) result(m_out)
    integer                               :: i
    complex(8), dimension(Ns,L)           :: m_in
    complex(8), dimension(Nindip,L)       :: m_out
    do i=1,Nindip
       m_out(i,:)=m_in(indipsites(i),:)
    enddo
  end function reshuffled


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
       if(nindex1+nindex==0)then       !avoid loop back and back
          ndelta=ndelta1/2.d0
          !deltan=deltan/2.d0         !by decreasing the step         
       endif
       a0trap=a0trap-chitrap*real(nindex,8)*ndelta  ! a0trap=-mu_tot

       !         in this way chitrap is the inverse trap compressibility chitrap=dmu/d(N_tot)

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
end program ahm_matsubara_trap
