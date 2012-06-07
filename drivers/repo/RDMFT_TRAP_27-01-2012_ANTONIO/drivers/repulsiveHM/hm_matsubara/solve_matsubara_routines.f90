!******************************************************************
!******************************************************************


subroutine setup_initial_sigma()
  logical :: check1,check2,check
  inquire(file="LSigma.ipt",exist=check)
  if(check)then
     if(mpiID==0)then
        write(*,*)"Reading Sigma in input:"
        call sread("LSigma.ipt",sigma(1:Ns,1:L))
     endif
     call MPI_BCAST(sigma,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  else
     n=0.5d0 ; sigma=zero!u*(n-0.5d0)
  endif
end subroutine setup_initial_sigma



!******************************************************************
!******************************************************************



subroutine get_gloc_mpi() 
  complex(8) :: zeta,Gloc(Ns,Ns),gf_tmp(Ns,L)
  integer    :: i,is
  call msg("Get local GF: (ETA --> fort.999)",id=0)
  call start_timer
  fg=zero ; gf_tmp=zero
  do i=1+mpiID,L,mpiSIZE
     Gloc  = zero
     Gloc  = -H0
     do is=1,Ns
        Gloc(is,is) = xi*wm(i) + xmu - sigma(is,i) - erandom(is)
     enddo
     call mat_inversion_sym(Gloc,Ns)
     do is=1,Ns
        gf_tmp(is,i) = Gloc(is,is)
     enddo
     call eta(i,L,999)
  enddo
  call stop_timer
  call MPI_REDUCE(gf_tmp,fg,Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
  call MPI_BCAST(fg,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
end subroutine get_gloc_mpi


!******************************************************************
!******************************************************************




subroutine solve_impurity_mpi()
  integer    :: is
  complex(8) :: zsigma(Ns,1:L)
  call msg("Solve impurity: (ETA --> fort.998)")
  if(Wdis/=0.d0)then
     call start_timer
     zsigma=zero ; sigma=zero !sigma must be set to zero or mpi_red will ends summing old sigmas
     do is=1+mpiID,Ns,mpiSIZE
        call solve_per_site(is)
        call eta(is,Ns,998)
     enddo
     call stop_timer
     call MPI_REDUCE(sigma,zsigma,Ns*L,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,MPIerr)
     call MPI_BCAST(zsigma,Ns*L,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
     sigma=zsigma
  else
     call solve_per_site(is=1)
     forall(is=2:Ns)sigma(is,:)=sigma(1,:)
  endif
end subroutine solve_impurity_mpi


!******************************************************************
!******************************************************************


subroutine solve_per_site(is)
  integer                     :: is
  complex(8),allocatable,save :: Sold(:,:)
  isite=is
  if(.not.allocated(Sold))allocate(Sold(Ns,1:L))
  sold(is,:)  =  sigma(is,:)

  call fftgf_iw2tau(fg(is,:),fgt(is,:),beta)
  n=-real(fgt(is,L))
  gamma = one/(one/fg(is,:) + sigma(is,:))
  xmu0=0.d0
  x(1)=xmu
  call broydn(x,check)
  xmu0=x(1)
  sigma(is,:) =  solve_mpt_matsubara(fg0,n,n0,xmu0)
  sigma(is,:) =  weigth*sigma(is,:) + (1.d0-weigth)*sold(is,:)
end subroutine solve_per_site


!******************************************************************
!******************************************************************


subroutine print_out(success)
  real(8)   :: nimp,nii(Ns)
  complex(8):: afg(1:L),asig(1:L)
  logical   :: success
  integer   :: i,is
  character(len=4) :: loop
  if(mpiID==0)then
     nimp=0.d0
     do is=1,Ns
        call fftgf_iw2tau(fg(is,1:L),fgt(is,0:L),beta)
        nii(is) =-2.d0*real(fgt(is,L))
     enddo
     nimp=sum(nii)/real(Ns,8)
     print*,"nimp  =",nimp
     call splot(trim(adjustl(trim(name_dir)))//"/navVSiloop.ipt",iloop,nimp,append=TT)


     write(loop,"(I4)")iloop
     do i=1,Ns
        call splot("Gloc_iw_site."//trim(adjustl(trim(loop)))//".ipt",wm,fg(i,1:L),append=TT)
        call splot("Sigma_iw_site."//trim(adjustl(trim(loop)))//".ipt",wm,sigma(i,1:L),append=TT)
     enddo

     if(success)then
        call splot(trim(adjustl(trim(name_dir)))//"/nVSisite.ipt",nii)
        call splot(trim(adjustl(trim(name_dir)))//"/erandomVSisite.ipt",erandom)
        call splot(trim(adjustl(trim(name_dir)))//"/LSigma.ipt",sigma)
        call splot(trim(adjustl(trim(name_dir)))//"/LG.ipt",fg)
        afg(:)  = sum(fg(1:Ns,1:L),dim=1)/dble(Ns)
        asig(:) = sum(sigma(1:Ns,1:L),dim=1)/dble(Ns)
        call splot(trim(adjustl(trim(name_dir)))//"/Gav_iw.ipt",wm,afg)
        call splot(trim(adjustl(trim(name_dir)))//"/Sigmaav_iw.ipt",wm,asig)
     endif
  endif
  return
end subroutine print_out

