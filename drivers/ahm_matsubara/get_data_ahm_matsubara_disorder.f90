program get_data
  USE RDMFT_VARS_GLOBAL
  USE TOOLS
  USE MATRIX
  USE FGSL
  implicit none
  logical                                 :: check
  integer                                 :: list_idum_length
  integer,allocatable,dimension(:)        :: list_idum
  integer                                 :: iread,row,col,is
  character(len=16)                       :: DIR
  real(8),allocatable,dimension(:)        :: nii,dii,rii,sii,zii,cdwnii
  !
  real(8),allocatable,dimension(:,:)      :: dij,adij,gdij
  real(8),allocatable,dimension(:,:)      :: nij,anij,gnij
  real(8),allocatable,dimension(:,:)      :: rij,arij,grij
  real(8),allocatable,dimension(:,:)      :: sij,asij,gsij
  real(8),allocatable,dimension(:,:)      :: zij,azij,gzij
  !
  real(8),allocatable,dimension(:,:)      :: ddn,dddelta
  real(8),allocatable,dimension(:,:)      :: cov_nd
  !
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma
  complex(8),allocatable,dimension(:,:)   :: fg0
  complex(8),allocatable,dimension(:,:)   :: avG,avS
  !
  complex(8),allocatable,dimension(:,:,:) :: afg,asigma
  !
  logical                                 :: wgf
  real(8)                                 :: r,delta,delta0


  !READ INPUT FILES:
  !=====================================================================
  wgf=.false.
  n_command_arg=command_argument_count()
  if(n_command_arg/=0)then
     do i=1,n_command_arg
        call get_command_argument(i,arg_buffer)
        call cmd_help(arg_buffer,"get_data_ahmmpt_matsubara_disorder wgf=[.false.]")
        call cmd_var(arg_buffer)
        if(nml_name=="NOGF")read(nml_value,*)wgf
     enddo
  endif
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

  if(wgf)then
     allocate(fg(2,Ns,L),fg0(2,L),sigma(2,Ns,L))
     allocate(afg(2,Ns,L),asigma(2,Ns,L))
     allocate(avG(2,L),avS(2,L))
     allocate(sii(Ns),rii(Ns),zii(Ns))
     allocate(rij(Nside,Nside),arij(Nside,Nside),grij(Nside,Nside))
     allocate(sij(Nside,Nside),asij(Nside,Nside),gsij(Nside,Nside))
     allocate(zij(Nside,Nside),azij(Nside,Nside),gzij(Nside,Nside))
  else
     allocate(nii(Ns),dii(Ns),cdwnii(Ns))
     allocate(dij(Nside,Nside),adij(Nside,Nside),gdij(Nside,Nside))
     allocate(nij(Nside,Nside),anij(Nside,Nside),gnij(Nside,Nside))
  endif


  !     BUILD THE LATTICE HAMILTONIAN:
  !=====================================================================
  call get_tb_hamiltonian

  !INQUIRE && READ RANDOM NUMBER SEEDS:
  inquire(file="list_idum",exist=check)
  if(.not.check)call abort("I can not find list_idum! ciao")
  list_idum_length=file_length("list_idum")
  allocate(list_idum(list_idum_length))
  call sread("list_idum",list_idum)

  !INIT AVERAGE ARRAYS:

  if(wgf)then
     arij= zero ; grij=zero 
     asij= zero ; gsij=zero
     azij= zero ; gzij=zero
     afg = zero ; asigma=zero
  else
     adij= zero ; gdij=zero
     anij= zero ; gnij=zero
  endif

  !CREATE DIRECTORY FOR THE AVERAGED RESULTS:
  DIR="RESULTS" ; call create_data_dir(DIR)

  !START LOOP OVER DISORDER REALIZATIONS: 
  do iread=1,list_idum_length
     idum=list_idum(iread)
     write(name_dir,"(I12)")idum
     name_dir="idum_"//trim(adjustl(trim(name_dir)))

     if(wgf)then
        call read_gf(trim(adjustl(trim(name_dir))))
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
        arij = arij + rij ; grij=grij+log(rij) 
        asij = asij + sij ; gsij=gsij+log(sij) 
        azij = azij + zij ; gzij=gzij+log(zij) 

        !BUILD ARITHEMTIC AVERAGES OVER DISORDER OF LOCAL GF:
        afg   = afg + fg
        asigma= asigma+sigma
     else
        call read_n_delta(trim(adjustl(trim(name_dir))))
        !BUILD ARITHEMTIC & GEOMETRIC AVERAGES OVER DISORDER:
        adij = adij + dij ; gdij=gdij+log(dij) 
        anij = anij + nij ; gnij=gnij+log(nij) 
        call splot("cdwnVSisite.data",cdwnii)
     endif

     call system("mv -vf *.data *.data.gz "//trim(adjustl(trim(name_dir)))//"/ 2>/dev/null")
  enddo
  close(100)



  if(wgf)then
     !PLOT \EE(F,G,Sigma,Self):
     afg=afg/real(list_idum_length,8)
     asigma=asigma/real(list_idum_length,8)
     do is=1,Ns
        call splot("G.arithmetic_iw.data",wm,afg(1,is,1:L),append=.true.)
        call splot("F.arithmetic_iw.data",wm,afg(2,is,1:L),append=.true.)
        call splot("Sigma.arithmetic_iw.data",wm,asigma(1,is,1:L),append=.true.)
        call splot("Self.arithmetic_iw.data",wm,asigma(2,is,1:L),append=.true.)
     enddo
     arij=arij/real(list_idum_length,8) ; grij=exp(grij/real(list_idum_length,8))
     asij=asij/real(list_idum_length,8) ; gsij=exp(gsij/real(list_idum_length,8))
     azij=azij/real(list_idum_length,8) ; gzij=exp(gzij/real(list_idum_length,8))
     call splot("sigVSij_arithmetic.data",asij)
     call splot("zetaVSij_arithmetic.data",azij)
     call splot("rhoVSij_arithmetic.data",arij)
     call splot("sigVSij_geometric.data",gsij)
     call splot("zetaVSij_geometric.data",gzij)
     call splot("rhoVSij_geometric.data",grij)
  else
     adij=adij/real(list_idum_length,8) ; gdij=exp(gdij/real(list_idum_length,8))
     anij=anij/real(list_idum_length,8) ; gnij=exp(gnij/real(list_idum_length,8))
     call splot("nVSij_arithmetic.data",anij)
     call splot("deltaVSij_arithmetic.data",adij)
     call splot("nVSij_geometric.data",gnij)
     call splot("deltaVSij_geometric.data",gdij)
     call get_covariance(anij,adij)
     call get_correlations()
  endif


  !MOVE EVERYTHING INTO RESULTS/
  call system("mv -vf *.data "//trim(adjustl(trim(DIR)))//"/ 2>/dev/null")

contains

  subroutine read_n_delta(idum_dir)
    character(len=*) :: idum_dir
    real(8)          :: w
    call sread(idum_dir//"/nVSisite.ipt",nii)
    call sread(idum_dir//"/deltaVSisite.ipt",dii)
    do row=1,Nside
       do col=1,Nside
          i            = ij2site(row,col)
          nij(row,col) = nii(i)
          dij(row,col) = dii(i)
          cdwnii(i)    = (-1.d0)**(row+col)*nij(row,col)
       enddo
    enddo
  end subroutine read_n_delta


  !*************************************************************************************
  !*************************************************************************************
  !*************************************************************************************


  subroutine read_gf(idum_dir)
    character(len=*) :: idum_dir
    real(8)          :: w
    call sread(idum_dir//"/LSigma.ipt",sigma(1,1:Ns,1:L))
    call sread(idum_dir//"/LSelf.ipt",sigma(2,1:Ns,1:L))
    call sread(idum_dir//"/LG.ipt",fg(1,1:Ns,1:L))
    call sread(idum_dir//"/LF.ipt",fg(2,1:Ns,1:L))
    do row=1,Nside
       do col=1,Nside
          i            = ij2site(row,col)
          sij(row,col) = dimag(sigma(1,i,1))-&
               wm(1)*(dimag(sigma(1,i,2))-dimag(sigma(1,i,1)))/(wm(2)-wm(1))
          rij(row,col) = dimag(fg(1,i,1))-&
               wm(1)*(dimag(fg(1,i,2))-dimag(fg(1,i,1)))/(wm(2)-wm(1))
          zij(row,col) = 1.d0/( 1.d0 + abs( dimag(sigma(1,i,1))/wm(1) ))
       enddo
    enddo
    rij=abs(rij)
    sij=abs(sij)
    zij=abs(zij)
  end subroutine read_gf


  !*************************************************************************************
  !*************************************************************************************
  !*************************************************************************************



  subroutine get_covariance(anij_,adij_)
    real(8),dimension(Nside,Nside) :: anij_,adij_
    character(len=20)              :: idum_dir
    real(8),dimension(Ns)          :: nbar,dbar
    real(8),dimension(Ns)          :: ni,di,den,ded
    integer                        :: i,j,is,js,ird,idum
    real(8),dimension(Ns,Ns)       :: cov_dn_dD,cov_dn_dn,cov_dD_dD
    real(8),dimension(Ns)          :: var_dn_dd,var_dn_dn,var_dd_dd,r_local
    real(8),dimension(Ns,Ns)       :: idum_fluct_dn_dD,idum_fluct_dn_dn,idum_fluct_dD_dD
    do is=1,Ns
       i=irow(is);j=icol(is)
       nbar(is)=anij_(i,j)
       dbar(is)=adij_(i,j)
    enddo
    cov_dn_dD=0.d0 ; cov_dn_dn=0.d0 ; cov_dD_dD=0.d0
    !loop over all realizations
    do ird=1,list_idum_length
       idum=list_idum(ird)
       write(idum_dir,"(I12)")idum
       idum_dir="idum_"//trim(adjustl(trim(idum_dir)))
       call sread(trim(adjustl(trim(idum_dir)))//"/nVSisite.ipt",ni)
       call sread(trim(adjustl(trim(idum_dir)))//"/deltaVSisite.ipt",Di)
       den  = ni-nbar ; deD=Di-Dbar
       !Build-up fluctuations matrix:
       do is=1,Ns
          do js=1,Ns
             idum_fluct_dn_dD(is,js)=den(is)*deD(js)
             idum_fluct_dn_dn(is,js)=den(is)*den(js)
             idum_fluct_dD_dD(is,js)=deD(is)*deD(js)
          enddo
       enddo
       forall(is=1:Ns)
          var_dn_dd(is)=idum_fluct_dn_dd(is,is)
          var_dn_dn(is)=idum_fluct_dn_dn(is,is)
          var_dd_dd(is)=idum_fluct_dd_dd(is,is)
       end forall
       cov_dn_dd=cov_dn_dD + idum_fluct_dn_dD/real(list_idum_length,8)
       cov_dn_dn=cov_dn_dn + idum_fluct_dn_dn/real(list_idum_length,8)
       cov_dD_dD=cov_dD_dD + idum_fluct_dD_dD/real(list_idum_length,8)
    enddo
    forall(is=1:Ns)
       var_dn_dd(is)=cov_dn_dd(is,is)
       var_dn_dn(is)=cov_dn_dn(is,is)
       var_dd_dd(is)=cov_dd_dd(is,is)
    end forall
    forall(is=1:Ns)r_local(is) = var_dn_dd(is)/sqrt(var_dn_dn(is)*var_dd_dd(is))
    call splot("covariance_dn_dd.local.data",var_dn_dd)
    call splot("correlation_dn_dd.local.data",r_local)
  end subroutine get_covariance

  !*************************************************************************************
  !*************************************************************************************
  !*************************************************************************************



  subroutine get_correlations()
    character(len=20)              :: idum_dir
    real(8),dimension(Ns)          :: nbar,dbar
    real(8),dimension(Ns)          :: ni,di,den,ded
    integer                        :: i,j,is,js,ird,idum
    real(8),dimension(Ns,Ns)       :: cov_dn_dD,cov_dn_dn,cov_dD_dD
    real(8),dimension(Ns)          :: var_dn_dd,var_dn_dn,var_dd_dd,r_local
    real(8),dimension(Ns,Ns)       :: idum_fluct_dn_dD,idum_fluct_dn_dn,idum_fluct_dD_dD
    cov_dn_dD=0.d0 ; cov_dn_dn=0.d0 ; cov_dD_dD=0.d0
    !loop over all realizations
    do ird=1,list_idum_length
       idum=list_idum(ird)
       write(idum_dir,"(I12)")idum
       idum_dir="idum_"//trim(adjustl(trim(idum_dir)))
       call sread(trim(adjustl(trim(idum_dir)))//"/nVSisite.ipt",ni)
       call sread(trim(adjustl(trim(idum_dir)))//"/deltaVSisite.ipt",Di)
       do is=1,Ns
          do js=1,Ns
             idum_fluct_dn_dD(is,js)=ni(is)*di(js)
             idum_fluct_dn_dn(is,js)=ni(is)*ni(js)
             idum_fluct_dD_dD(is,js)=di(is)*di(js)
          enddo
       enddo
       cov_dn_dd=cov_dn_dD + idum_fluct_dn_dD/real(list_idum_length,8)
       cov_dn_dn=cov_dn_dn + idum_fluct_dn_dn/real(list_idum_length,8)
       cov_dD_dD=cov_dD_dD + idum_fluct_dD_dD/real(list_idum_length,8)
    enddo
    forall(is=1:Ns,js=1:Ns)
       cov_dn_dd(is,js)=cov_dn_dd(is,js)/cov_dn_dd(is,is)
       cov_dn_dn(is,js)=cov_dn_dn(is,js)/cov_dn_dn(is,is)
       cov_dd_dd(is,js)=cov_dd_dd(is,js)/cov_dd_dd(is,is)
    end forall
    call splot("correlation_ninj.data",cov_dn_dn(:,1))
    call splot("correlation_didj.data",cov_dd_dd(:,1))
    call splot("correlation_nidj.data",cov_dn_dd(:,1))
  end subroutine get_correlations



  !****************************************************************************
  !****************************************************************************
  !**************************************************************************** 


end program get_data
