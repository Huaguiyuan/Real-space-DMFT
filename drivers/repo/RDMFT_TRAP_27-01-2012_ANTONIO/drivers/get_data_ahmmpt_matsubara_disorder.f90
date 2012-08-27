program get_data
  USE RDMFT_VARS_GLOBAL
  USE TOOLS
  USE MATRIX
  USE FGSL
  implicit none
  logical                                 :: check
  integer                                 :: list_idum_length
  integer,allocatable                     :: list_idum(:)
  integer                                 :: iread,row,col,is
  character(len=16)                       :: DIR
  real(8),allocatable                     :: nii(:),dii(:),rii(:),sii(:),zii(:)
  real(8),allocatable                     :: grid_x(:),grid_y(:)
  !
  real(8),allocatable                     :: dij(:,:),adij(:,:),gdij(:,:)
  real(8),allocatable                     :: nij(:,:),anij(:,:),gnij(:,:)
  real(8),allocatable                     :: rij(:,:),arij(:,:),grij(:,:)
  real(8),allocatable                     :: sij(:,:),asij(:,:),gsij(:,:)
  real(8),allocatable                     :: zij(:,:),azij(:,:),gzij(:,:)
  !
  real(8),allocatable                     :: ddn(:,:),dddelta(:,:)
  real(8),allocatable                     :: cov_nd(:,:)
  !
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma,gf_tmp
  complex(8),allocatable,dimension(:,:)   :: fg0
  complex(8),allocatable,dimension(:,:)   :: avG,avS
  !
  complex(8),allocatable,dimension(:,:,:) :: afg,asigma
  !
  logical                                 :: converged
  real(8)                                 :: r,delta,delta0


  !READ INPUT FILES:
  !=====================================================================
  call read_input("inputIPT.in")
  call rdmft_read_input("inputRDMFT.in")
  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)
  write(*,"(A,I9,A)")"Using ",L," frequencies"

  !     ALLOCATE WORKING ARRAYS:
  !=====================================================================
  Ns    =Nside**2
  n     =0.5d0
  wmax  =wmax+Wdis
  allocate(erandom(Ns))
  allocate(H0(Ns,Ns))
  allocate(icol(Ns),irow(Ns))
  allocate(ij2site(Nside,Nside))
  !
  allocate(fg(2,Ns,L))
  allocate(fg0(2,L))
  allocate(sigma(2,Ns,L))
  allocate(gf_tmp(2,Ns,L))
  !
  allocate(afg(2,Ns,L),asigma(2,Ns,L))
  !
  allocate(nii(Ns),dii(Ns),sii(Ns),rii(Ns),zii(Ns))
  !
  allocate(dij(Nside,Nside),adij(Nside,Nside),gdij(Nside,Nside))
  allocate(nij(Nside,Nside),anij(Nside,Nside),gnij(Nside,Nside))
  allocate(rij(Nside,Nside),arij(Nside,Nside),grij(Nside,Nside))
  allocate(sij(Nside,Nside),asij(Nside,Nside),gsij(Nside,Nside))
  allocate(zij(Nside,Nside),azij(Nside,Nside),gzij(Nside,Nside))
  !
  allocate(grid_x(Nside),grid_y(Nside))
  !
  allocate(avG(2,L),avS(2,L))

  !     BUILD THE LATTICE HAMILTONIAN:
  !=====================================================================
  call get_tb_hamiltonian

  !INQUIRE && READ RANDOM NUMBER SEEDS:
  inquire(file="list_idum",exist=check)
  if(.not.check)call abort("I can not find list_idum! ciao")
  list_idum_length=file_length("list_idum")
  allocate(list_idum(list_idum_length))
  call sread("list_idum",list_idum)

  !BUILD A GRID FOR  LATTICE PLOTS:
  do row=1,Nside
     grid_x(row)=dble(row)
     grid_y(row)=dble(row)
  enddo

  !INIT AVERAGE ARRAYS:
  adij= zero ; gdij=zero
  anij= zero ; gnij=zero
  arij= zero ; grij=zero 
  asij= zero ; gsij=zero
  azij= zero ; gzij=zero
  afg = zero ; asigma=zero

  !CREATE DIRECTORY FOR THE AVERAGED RESULTS:
  DIR="RESULTS" ; call create_data_dir(DIR)

  !START LOOP OVER DISORDER REALIZATIONS: 
  do iread=1,list_idum_length
     idum=list_idum(iread)
     write(name_dir,"(I12)")idum
     name_dir="idum_"//trim(adjustl(trim(name_dir)))

     !READ IDUM-REALIZED Q.TIES: [N,DELTA,G,SIGMA]
     call read_all(trim(adjustl(trim(name_dir))))

     do row=1,Nside
        do col=1,Nside
           i      = ij2site(row,col)
           sii(i) = sij(row,col)
           rii(i) = rij(row,col)
           zii(i) = zij(row,col)
        enddo
     enddo
     call splot("sigVSisite.data",sii)
     call splot("rhoVSisite.data",rii)
     call splot("zVSisite.data",zii)


     avG(:,:) = sum(fg(1:2,1:Ns,1:L),dim=2)/dble(Ns)
     avS(:,:) = sum(sigma(1:2,1:Ns,1:L),dim=2)/dble(Ns)
     call splot("lat.ave.G_iw.data",wm,avG(1,1:L))
     call splot("lat.ave.F_iw.data",wm,avG(2,1:L))
     call splot("lat.ave.Sigma_iw.data",wm,avS(1,1:L))
     call splot("lat.ave.Self_iw.data",wm,avS(2,1:L))


     !BUILD ARITHEMTIC & GEOMETRIC AVERAGES OVER DISORDER:
     adij = adij + dij ; gdij=gdij+log(dij) 
     anij = anij + nij ; gnij=gnij+log(nij) 
     arij = arij + rij ; grij=grij+log(rij) 
     asij = asij + sij ; gsij=gsij+log(sij) 
     azij = azij + zij ; gzij=gzij+log(zij) 

     !BUILD ARITHEMTIC AVERAGES OVER DISORDER OF LOCAL GF:
     afg   = afg + fg
     asigma= asigma+sigma

     call system("mv -vf *.data "//trim(adjustl(trim(name_dir)))//"/ 2>/dev/null")
  enddo
  close(100)

  !PLOT \EE(F,G,Sigma,Self):
  afg=afg/real(list_idum_length,8)
  asigma=asigma/real(list_idum_length,8)
  do is=1,Ns
     call splot("G.arithmetic_iw.data",wm,afg(1,is,1:L),append=.true.)
     call splot("F.arithmetic_iw.data",wm,afg(2,is,1:L),append=.true.)
     call splot("Sigma.arithmetic_iw.data",wm,asigma(1,is,1:L),append=.true.)
     call splot("Self.arithmetic_iw.data",wm,asigma(2,is,1:L),append=.true.)
  enddo

  !PLOT \EE(n,d,rho,sigma)_arit,geom
  adij=adij/real(list_idum_length,8) ; gdij=exp(gdij/real(list_idum_length,8))
  anij=anij/real(list_idum_length,8) ; gnij=exp(gnij/real(list_idum_length,8))
  arij=arij/real(list_idum_length,8) ; grij=exp(grij/real(list_idum_length,8))
  asij=asij/real(list_idum_length,8) ; gsij=exp(gsij/real(list_idum_length,8))
  azij=azij/real(list_idum_length,8) ; gzij=exp(gzij/real(list_idum_length,8))
  call splot("3d_nVSij_arithmetic.data",grid_x,grid_y,anij)
  call splot("3d_deltaVSij_arithmetic.data",grid_x,grid_y,adij)
  call splot("3d_sigVSij_arithmetic.data",grid_x,grid_y,asij)
  call splot("3d_zetaVSij_arithmetic.data",grid_x,grid_y,azij)
  call splot("3d_rhoVSij_arithmetic.data",grid_x,grid_y,arij)
  call splot("3d_nVSij_geometric.data",grid_x,grid_y,gnij)
  call splot("3d_deltaVSij_geometric.data",grid_x,grid_y,gdij)
  call splot("3d_sigVSij_geometric.data",grid_x,grid_y,gsij)
  call splot("3d_zetaVSij_geometric.data",grid_x,grid_y,gzij)
  call splot("3d_rhoVSij_geometric.data",grid_x,grid_y,grij)

  !BUILD UP THE COVARIANCE MATRIX
  call msg("Get Covariance:")
  call get_covariance(anij,adij)

  !MOVE EVERYTHING INTO RESULTS/
  call system("mv -vf *.data "//trim(adjustl(trim(DIR)))//"/ 2>/dev/null")


contains


  !*************************************************************************************
  !*************************************************************************************
  !*************************************************************************************


  subroutine get_covariance(anij_,adij_)
    real(8),dimension(Nside,Nside) :: anij_,adij_  !input \bar{n(i,j)},\bar{\delta(i,j)}
    real(8),dimension(Ns)          :: nbar,dbar
    real(8),dimension(Ns)          :: den,ded
    integer                        :: ird,is,js,rowi,coli,rowj,colj
    character(len=20)              :: idum_dir
    real(8)                        :: cov_dn_dD(Ns,Ns),cov_dn_dn(Ns,Ns),cov_dD_dD(Ns,Ns)
    real(8)                        :: var_dn_dd(Ns),var_dn_dn(Ns),var_dd_dd(Ns)
    real(8)                        :: r_correlation(Ns,Ns),r_local(Ns)
    !
    real(8),dimension(3*Ns)        :: work
    real(8),dimension(Ns)          :: evalues
    integer                        :: info,lwork
    !
    !FIRST: transform the inputs from ij-grid to Ns-array
    do is=1,Ns
       rowi=irow(is)
       coli=icol(is)
       nbar(is)=anij_(rowi,coli)
       dbar(is)=adij_(rowi,coli)
    enddo
    !SECOND: loop over disorder realizations and build the covariance matrix
    !this is done in three steps: a) build the fluctuations df_i=(f-fbar)_i
    !b) build the products df_i*df_j
    !c) sum up the matrix C_ij^(n) = (df_i*df_j)^(n)
    cov_dn_dd=0.d0
    cov_dn_dn=0.d0
    cov_dd_dd=0.d0
    do ird=1,list_idum_length
       idum=list_idum(ird)
       write(idum_dir,"(I12)")idum
       idum_dir="idum_"//trim(adjustl(trim(idum_dir)))
       call sread(trim(adjustl(trim(idum_dir)))//"/nVSisite.ipt",nii)
       call sread(trim(adjustl(trim(idum_dir)))//"/deltaVSisite.ipt",dii)
       den=nii-nbar
       ded=dii-dbar
       do is=1,Ns
          do js=1,Ns
             cov_dn_dd(is,js)=cov_dn_dd(is,js)+den(is)*ded(js)
             cov_dn_dn(is,js)=cov_dn_dn(is,js)+den(is)*den(js)
             cov_dd_dd(is,js)=cov_dd_dd(is,js)+ded(is)*ded(js)
          enddo
       enddo
    enddo
    cov_dn_dd=cov_dn_dd/real(list_idum_length,8)
    cov_dn_dn=cov_dn_dn/real(list_idum_length,8)
    cov_dd_dd=cov_dd_dd/real(list_idum_length,8)

    forall(is=1:Ns)var_dn_dd(is)=cov_dn_dd(is,is)
    forall(is=1:Ns)var_dn_dn(is)=cov_dn_dn(is,is)
    forall(is=1:Ns)var_dd_dd(is)=cov_dd_dd(is,is)
    call splot("cov_dn_dd.diagonal.data",var_dn_dd)
    call splot("cov_dn_dn.diagonal.data",var_dn_dn)
    call splot("cov_dd_dd.diagonal.data",var_dd_dd)

    forall(is=1:Ns,js=1:Ns)&
         r_correlation(is,js) = cov_dn_dd(is,js)/sqrt(var_dn_dn(is)*var_dd_dd(js))
    forall(is=1:Ns)r_local(is)=r_correlation(is,is)
    call splot("r_correlation.local.data",r_local)

    call dsyev('V','U',Ns,cov_dn_dd,Ns,evalues,work,size(work),info)
    call splot("cov_dn_dd.eigenvalues.data",evalues)

    call dsyev('V','U',Ns,cov_dn_dn,Ns,evalues,work,size(work),info)
    call splot("cov_dn_dn.eigenvalues.data",evalues)

    call dsyev('V','U',Ns,cov_dd_dd,Ns,evalues,work,size(work),info)
    call splot("cov_dd_dd.eigenvalues.data",evalues)

  end subroutine get_covariance


  !*************************************************************************************
  !*************************************************************************************
  !*************************************************************************************


  subroutine read_all(idum_dir)
    character(len=*) :: idum_dir
    real(8)          :: w
    call sread(idum_dir//"/nVSisite.ipt",nii)
    call sread(idum_dir//"/deltaVSisite.ipt",dii)
    call sread(idum_dir//"/LSigma.ipt",sigma(1,1:Ns,1:L))
    call sread(idum_dir//"/LSelf.ipt",sigma(2,1:Ns,1:L))
    call sread(idum_dir//"/LG.ipt",fg(1,1:Ns,1:L))
    call sread(idum_dir//"/LF.ipt",fg(2,1:Ns,1:L))
    do row=1,Nside
       do col=1,Nside
          i            = ij2site(row,col)
          nij(row,col) = nii(i)
          dij(row,col) = dii(i)
          sij(row,col) = dimag(sigma(1,i,1))-&
               wm(1)*(dimag(sigma(1,i,2))-dimag(sigma(1,i,1)))/(wm(2)-wm(1))
          rij(row,col) = dimag(fg(1,i,1))-&
               wm(1)*(dimag(fg(1,i,2))-dimag(fg(1,i,1)))/(wm(2)-wm(1))
          zij          = 1.d0/( 1.d0 + abs( dimag(sigma(1,i,1))/wm(1) ))
       enddo
    enddo
    rij=abs(rij)
    sij=abs(sij)
    zij=abs(zij)
  end subroutine read_all


  !****************************************************************************
  !****************************************************************************
  !****************************************************************************


  ! subroutine get_histogram(y,nbin,filename)
  !   real(8),dimension(:)            :: y
  !   character(len=*)                :: filename
  !   character(len=10)               :: cnbin
  !   character(len=len(filename)+30) :: pname
  !   integer                         :: i,nbin,Ls
  !   real(fgsl_double)               :: a, b, x
  !   integer(fgsl_size_t)            :: n
  !   integer(fgsl_int)               :: status
  !   type(fgsl_histogram)            :: h
  !   type(fgsl_file)                 :: pfile
  !   real(8),allocatable,dimension(:):: ra,rb
  !   integer,allocatable,dimension(:):: hh
  !   !
  !   Ls     = size(y)   ; n = nbin ;  write(cnbin,"(I10)")nbin
  !   pname  = trim(adjustl(trim(filename)))//"_"//trim(adjustl(trim(cnbin)))
  !   a      = minval(y) ; b = maxval(y)
  !   pfile  = fgsl_open(trim(adjustl(trim(pname)))//".data",'a')
  !   h      = fgsl_histogram_alloc(n)
  !   status = fgsl_histogram_set_ranges_uniform(h, a, b)
  !   do i=1,Ls
  !      x=y(i)
  !      status = fgsl_histogram_increment (h, x)
  !   enddo
  !   status = fgsl_histogram_fprintf(pfile, h,'%12.7f','%7.0f')
  !   call fgsl_histogram_free(h)
  !   status = fgsl_close(pfile)
  !   !
  !   allocate(ra(nbin),rb(nbin),hh(nbin))
  !   open(20,file=trim(adjustl(trim(pname)))//".data")
  !   open(30,file=trim(adjustl(trim(pname)))//".histogram")
  !   do i=1,nbin
  !      read(20,*)ra(i),rb(i),hh(i)
  !      write(30,*)ra(i),hh(i)
  !      write(30,*)rb(i),hh(i)
  !   enddo
  !   close(20);close(30)
  ! end subroutine get_histogram


end program get_data
