!########################################################
!PURPOSE  :solve the attractive (A) disordered (D) Hubbard
! model (HM) using Modified Perturbation THeory (MPT), w/ DMFT
! disorder realization depends on the parameter int idum:
! so that different realizations (statistics) are performed 
! calling this program many times providing a *different* seed 
! +IDUM. The result of each calculation is stored different dir
! indexed by the seed itself.
!AUTHORS  : A.Amaricci
!########################################################
program ahm_matsubara_disorder
  USE RDMFT_VARS_GLOBAL
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma,sigma_tmp
  complex(8),allocatable,dimension(:,:)   :: fg0
  logical                                 :: converged
  real(8)                                 :: r
  integer                                 :: i,is

  !GLOBAL INITIALIZATION:
  !===============================================================
  include "init_global_disorder.f90"


  !ALLOCATE WORKING ARRAYS
  !===============================================================
  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)

  write(*,"(A,I9,A)")"Using ",L," frequencies"
  allocate(fg(2,Ns,L))
  allocate(fg0(2,L))
  allocate(sigma(2,Ns,L))
  allocate(sigma_tmp(2,Ns,L))



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

     converged=check_convergence(sigma(1,:,:)+sigma(2,:,:),eps_error,Nsuccess,nloop,id=0)
     if(nread/=0.d0)call search_mu(converged)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
     call print_sc_out(converged)
     call end_loop()
  enddo
  if(mpiID==0)call system("mv -vf *.err "//trim(adjustl(trim(name_dir)))//"/")
  call close_mpi()



contains

  !******************************************************************
  !******************************************************************


  subroutine setup_initial_sc_sigma()
    logical :: check1,check2,check
    if(mpiID==0)then
       inquire(file="LSigma.ipt",exist=check1)
       inquire(file="LSelf.ipt",exist=check2)
       check=check1*check2
       if(check)then
          write(*,*)"Reading Sigma in input:"
          call sread("LSigma.ipt",sigma(1,1:Ns,1:L))
          call sread("LSelf.ipt",sigma(2,1:Ns,1:L))
       else
          sigma(1,:,:)=zero ; sigma(2,:,:)=-deltasc
       endif
    endif
    call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  end subroutine setup_initial_sc_sigma


  !******************************************************************
  !******************************************************************


  subroutine get_sc_gloc_mpi()
    complex(8) :: Gloc(2*Ns,2*Ns),gf_tmp(2,Ns,L)
    integer    :: is
    call msg("Get local GF:",id=0)
    call start_timer
    fg=zero ; gf_tmp=zero
    do i=1+mpiID,L,mpiSIZE
       Gloc=zero
       Gloc(1:Ns,1:Ns)          = -H0
       Gloc(Ns+1:2*Ns,Ns+1:2*Ns)=  H0
       do is=1,Ns
          Gloc(is,is)      =  xi*wm(i)+xmu-sigma(1,is,i)-erandom(is)
          Gloc(Ns+is,Ns+is)= -conjg(Gloc(is,is))
          Gloc(is,Ns+is)   = -sigma(2,is,i)
          Gloc(Ns+is,is)   = -sigma(2,is,i)!sigma(2,is,-i)!this should be a simmetry in Matsubara!
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



  subroutine solve_sc_impurity_mpi()
    integer :: is,Nchunk,R
    logical :: disorder
    disorder=.false. ; if(Wdis/=0)disorder=.true.
    call msg("Solve impurity:")
    disorder=.true.
    if(disorder)then
       call start_timer
       sigma_tmp=zero
       do is=1+mpiID,Ns,mpiSIZE
          call solve_per_site(is)
          call eta(is,Ns,998)
       enddo
       call stop_timer
       call MPI_REDUCE(sigma_tmp,sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
       call MPI_BCAST(sigma,2*Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
    else
       call solve_per_site(is=1)
       forall(is=1:Ns)sigma(:,is,:)=sigma_tmp(:,1,:)
    end if
  end subroutine solve_sc_impurity_mpi



  !******************************************************************
  !******************************************************************



  subroutine solve_per_site(is)
    integer                                      :: is
    complex(8)                                   :: det(L)
    complex(8),dimension(:,:,:),allocatable,save :: sold
    complex(8),dimension(2,L)                    :: calG
    real(8),dimension(2,0:L)                     :: fgt,fg0t
    real(8)                                      :: n,n0,delta,delta0
    if(.not.allocated(sold))allocate(sold(2,Ns,L))
    sold(:,is,:)  =  sigma(:,is,:)
    call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
    call fftgf_iw2tau(fg(2,is,:),fgt(2,0:L),beta,notail=.true.)
    n    = -real(fgt(1,L),8) ; delta= -u*fgt(2,0)
    fg0=zero ; calG=zero
    det       = abs(fg(1,is,:))**2    + fg(2,is,:)**2
    fg0(1,:)  = conjg(fg(1,is,:))/det + sigma(1,is,:) - U*(n-0.5d0)
    fg0(2,:)  = fg(2,is,:)/det        + sigma(2,is,:) + delta
    det       =  abs(fg0(1,:))**2 + fg0(2,:)**2
    calG(1,:) =  conjg(fg0(1,:))/det
    calG(2,:) =  fg0(2,:)/det
    call fftgf_iw2tau(calG(1,:),fg0t(1,:),beta)
    call fftgf_iw2tau(calG(2,:),fg0t(2,:),beta,notail=.true.)
    n0=-real(fg0t(1,L)) ; delta0= -u*fg0t(2,0)
    write(750,"(2I4,4(f16.12))",advance="yes")mpiID,is,n,n0,delta,delta0
    sigma_tmp(:,is,:) =  solve_mpt_sc_matsubara(calG,n,n0,delta,delta0)
    sigma_tmp(:,is,:) =  weigth*sigma_tmp(:,is,:) + (1.d0-weigth)*sold(:,is,:)
  end subroutine solve_per_site



  !******************************************************************
  !******************************************************************



  subroutine print_sc_out(converged)
    integer                  :: i,is
    real(8)                  :: nimp,delta
    real(8)                  :: nii(Ns),dii(Ns)
    logical                  :: converged
    character(len=4)         :: loop
    real(8),dimension(2,0:L) :: fgt

    if(mpiID==0)then
       nimp=0.d0 ; delta=0.d0
       do is=1,Ns
          call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
          call fftgf_iw2tau(fg(2,is,:),fgt(2,0:L),beta,notail=.true.)
          nii(is) = -2.d0*real(fgt(1,L))
          dii(is) = -u*fgt(2,0)
       enddo
       nimp = sum(nii)/dble(Ns)
       delta= sum(dii)/dble(Ns)
       print*,"nimp  =",nimp
       print*,"delta =",delta
       call splot(trim(adjustl(trim(name_dir)))//"/navVSiloop.ipt",iloop,nimp,append=TT)
       call splot(trim(adjustl(trim(name_dir)))//"/davVSiloop.ipt",iloop,delta,append=TT)

       call splot(trim(adjustl(trim(name_dir)))//"/LSigma_iw.ipt",sigma(1,1:Ns,1:L),wm(1:L))
       call splot(trim(adjustl(trim(name_dir)))//"/LSelf_iw.ipt",sigma(2,1:Ns,1:L),wm(1:L))
       call splot(trim(adjustl(trim(name_dir)))//"/LG_iw.ipt",fg(1,1:Ns,1:L),wm(1:L))
       call splot(trim(adjustl(trim(name_dir)))//"/LF_iw.ipt",fg(2,1:Ns,1:L),wm(1:L))

       if(converged)then
          call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.ipt",nii)
          call splot(trim(adjustl(trim(name_dir)))//"/deltaVSisite.ipt",dii)
          call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.ipt",erandom)
       endif
    end if
    return
  end subroutine print_sc_out



  !******************************************************************
  !******************************************************************



  subroutine search_mu(convergence)
    real(8)               :: naverage
    logical,intent(inout) :: convergence
    if(mpiID==0)then
       naverage=get_naverage(id=0)
       nindex1=nindex
       ndelta1=ndelta
       if((naverage >= nread+nerror))then
          nindex=-1
       elseif(naverage <= nread-nerror)then
          nindex=1
       else
          nindex=0
       endif
       if(nindex1+nindex==0)then !avoid loop forth and back
          ndelta=real(ndelta1/2.d0,8) !decreasing the step
          xmu=xmu+real(nindex,8)*ndelta
       else
          ndelta=ndelta1
          xmu=xmu+real(nindex,8)*ndelta
       endif
       write(*,"(f15.12,1x,f15.12,1x,f15.12)")xmu,real(nindex,8),ndelta
       write(*,"(A,f15.12,A,f15.12,A,f15.12)")" n=",naverage,"/",nread,"| ",ndelta
       print*,""
       print*,abs(naverage-nread),xmu
       print*,""
       if(abs(naverage-nread)>nerror)convergence=.false.
       call splot(trim(adjustl(trim(name_dir)))//"/muVSiter.ipt",iloop,xmu,abs(naverage-nread),append=.true.)
    endif
    call MPI_BCAST(xmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
  end subroutine search_mu


  !******************************************************************
  !******************************************************************


  function get_naverage(id) result(naverage)
    integer :: is,id
    real(8) :: naverage
    real(8) :: nii(Ns)
    real(8),dimension(2,0:L) :: fgt
    if(mpiID==id)then
       naverage=0.d0
       do is=1,Ns
          call fftgf_iw2tau(fg(1,is,:),fgt(1,0:L),beta)
          nii(is) = -2.d0*real(fgt(1,L))
       enddo
       naverage = sum(nii)/dble(Ns)
    endif
  end function get_naverage

  !******************************************************************
  !******************************************************************

end program ahm_matsubara_disorder
