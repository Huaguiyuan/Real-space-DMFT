program get_data
  USE RDMFT_VARS_GLOBAL
  USE TOOLS
  USE MATRIX
  implicit none
  logical                                 :: check,wgf
  integer                                 :: list_idum_length
  integer,allocatable,dimension(:)        :: list_idum
  integer                                 :: i,iread,row,col,is
  character(len=16)                       :: DIR
  real(8),allocatable,dimension(:)        :: rii,sii,zii,cdwnii
  !
  real(8),allocatable,dimension(:,:,:)      :: CRL
  real(8),allocatable,dimension(:,:)      :: cov_nd
  !
  complex(8),allocatable,dimension(:,:,:) :: fg,sigma
  !
  complex(8),allocatable,dimension(:,:) :: afg,asigma
  complex(8),allocatable,dimension(:,:) :: tfg,tsigma
  !
  real(8)                                 :: r


  !READ INPUT FILES:
  !=====================================================================
  call parse_cmd_variable(wgf,"WGF",default=.false.)
  call read_input("inputIPT.in")
  call rdmft_read_input("inputRDMFT.in")
  allocate(wm(L),tau(0:L))
  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  tau(0:)= linspace(0.d0,beta,L+1,mesh=dtau)
  write(*,"(A,I9,A)")"Using ",L," frequencies"


  !     ALLOCATE WORKING ARRAYS:
  !=====================================================================
  Ns    =Nside**2
  wmax  =wmax+Wdis
  allocate(erandom(Ns))
  allocate(H0(Ns,Ns))
  allocate(icol(Ns),irow(Ns))
  allocate(ij2site(Nside,Nside))

  allocate(nii(Ns),dii(Ns),cdwnii(Ns))
  allocate(sii(Ns),rii(Ns),zii(Ns))
  allocate(CRL(3,Ns,Ns))
  allocate(fg(2,Ns,L),sigma(2,Ns,L))
  allocate(afg(2,L),asigma(2,L))
  allocate(tfg(2,L),tsigma(2,L))


  !     BUILD THE LATTICE HAMILTONIAN:
  !=====================================================================
  call get_tb_hamiltonian


  !INQUIRE && READ RANDOM NUMBER SEEDS:
  inquire(file="list_idum",exist=check)
  if(.not.check)call abort("I can not find list_idum! ciao")
  list_idum_length=file_length("list_idum")
  allocate(list_idum(list_idum_length))
  call sread("list_idum",list_idum)


  !CREATE DIRECTORY FOR THE AVERAGED RESULTS:
  DIR="RESULTS" ; call create_data_dir(DIR)

  !START LOOP OVER DISORDER REALIZATIONS: 
  do iread=1,list_idum_length
     idum=list_idum(iread)
     write(name_dir,"(I12)")idum
     name_dir="idum_"//trim(adjustl(trim(name_dir)))

     call read_all(trim(adjustl(trim(name_dir))))


     open(10,file="CorrLenght_n.n.data")
     open(11,file="CorrLenght_delta.delta.data")
     open(11,file="CorrLenght_n.delta.data")
     open(12,file="CorrLenght_cdwn.cdwn.data")
     open(13,file="CorrLenght_zeta.zeta.data")

     do is=1,Ns
        do js=1,Ns
           CRL(1,is,js)=nii(is)*nii(js)
           CRL(2,is,js)=dii(is)*nii(js)
           CRL(3,is,js)=dii(is)*dii(js)
        enddo
     enddo

     do i=1,Nside
        is=ij2site(i,i)
        write(10,*)is,nii(is)
        write(11,*)is,dii(is)
        write(12,*)is,cdwnii(is)
        write(13,*)is,zii(is)
        if(i<Nside)then
           is=ij2site(i,i+1)
           write(10,*)is,nii(is)
           write(11,*)is,dii(is)
           write(12,*)is,cdwnii(is)
           write(13,*)is,zii(is)
        endif
     enddo
     do i=10,13
        close(i)
     enddo

     if(wgf)then
        afg   = afg    + sum(fg,dim=2)/dble(Ns)
        tfg   = tfg    + log(sum(fg,dim=2)/dble(Ns))
        asigma= asigma + sum(sigma,dim=2)/dble(Ns)
        tsigma= tsigma + log(sum(sigma,dim=2)/dble(Ns))
     endif

     call system("mv -vf *.data *.data.gz "//trim(adjustl(trim(name_dir)))//"/ 2>/dev/null")

  enddo
  close(100)



  !PLOT \EE(F,G,Sigma,Self):
  if(wgf)then
     afg=afg/real(list_idum_length,8)
     tfg=tfg/real(list_idum_length,8) ; tfg=exp(tfg)
     asigma=asigma/real(list_idum_length,8)
     tsigma=tsigma/real(list_idum_length,8);tsigma=exp(tsigma)
     call splot("aveG_iw.data",wm,afg(1,1:L))
     call splot("aveF_iw.data",wm,afg(2,1:L))
     call splot("aveSigma_iw.data",wm,asigma(1,1:L))
     call splot("aveSelf_iw.data",wm,asigma(2,1:L))
     call splot("typG_iw.data",wm,tfg(1,1:L))
     call splot("typF_iw.data",wm,tfg(2,1:L))
     call splot("typSigma_iw.data",wm,tsigma(1,1:L))
     call splot("typSelf_iw.data",wm,tsigma(2,1:L))
  endif

  !MOVE EVERYTHING INTO RESULTS/
  call system("mv -vf *.data "//trim(adjustl(trim(DIR)))//"/ 2>/dev/null")

contains

  function distance(is,js) result(dist)
    integer :: is,js
    real(8) :: dist
    dist = abs( (irow(is)-irow(js)) + (icol(is)-icol(js)) )
  end function distance


  subroutine read_all(idum_dir)
    character(len=*) :: idum_dir
    call sread(idum_dir//"/nVSisite.data",nii)
    call sread(idum_dir//"/deltaVSisite.data",dii)
    call sread(idum_dir//"/cdwVSisite.data",cdwnii)
    call sread(idum_dir//"/rhoVSisite.data",rii)
    call sread(idum_dir//"/sigmaVSisite.data",sii)
    call sread(idum_dir//"/zetaVSisite.data",zii)
    if(wgf)then
       call sread(idum_dir//"/LSigma_iw.data",sigma(1,1:Ns,1:L),wm(1:L))
       call sread(idum_dir//"/LSelf_iw.data",sigma(2,1:Ns,1:L),wm(1:L))
       call sread(idum_dir//"/LG_iw.data",fg(1,1:Ns,1:L),wm(1:L))
       call sread(idum_dir//"/LF_iw.data",fg(2,1:Ns,1:L),wm(1:L))
    endif
  end subroutine read_all



  !*************************************************************************************
  !*************************************************************************************
  !*************************************************************************************



  ! subroutine get_covariance(anij_,adij_)
  !   real(8),dimension(Nside,Nside) :: anij_,adij_
  !   character(len=20)              :: idum_dir
  !   real(8),dimension(Ns)          :: nbar,dbar
  !   real(8),dimension(Ns)          :: ni,di,den,ded
  !   integer                        :: i,j,is,js,ird,idum
  !   real(8),dimension(Ns,Ns)       :: cov_dn_dD,cov_dn_dn,cov_dD_dD
  !   real(8),dimension(Ns)          :: var_dn_dd,var_dn_dn,var_dd_dd,r_local
  !   real(8),dimension(Ns,Ns)       :: idum_fluct_dn_dD,idum_fluct_dn_dn,idum_fluct_dD_dD
  !   do is=1,Ns
  !      i=irow(is);j=icol(is)
  !      nbar(is)=anij_(i,j)
  !      dbar(is)=adij_(i,j)
  !   enddo
  !   cov_dn_dD=0.d0 ; cov_dn_dn=0.d0 ; cov_dD_dD=0.d0
  !   !loop over all realizations
  !   do ird=1,list_idum_length
  !      idum=list_idum(ird)
  !      write(idum_dir,"(I12)")idum
  !      idum_dir="idum_"//trim(adjustl(trim(idum_dir)))
  !      call sread(trim(adjustl(trim(idum_dir)))//"/nVSisite.ipt",ni)
  !      call sread(trim(adjustl(trim(idum_dir)))//"/deltaVSisite.ipt",Di)
  !      den  = ni-nbar ; deD=Di-Dbar
  !      !Build-up fluctuations matrix:
  !      do is=1,Ns
  !         do js=1,Ns
  !            idum_fluct_dn_dD(is,js)=den(is)*deD(js)
  !            idum_fluct_dn_dn(is,js)=den(is)*den(js)
  !            idum_fluct_dD_dD(is,js)=deD(is)*deD(js)
  !         enddo
  !      enddo
  !      forall(is=1:Ns)
  !         var_dn_dd(is)=idum_fluct_dn_dd(is,is)
  !         var_dn_dn(is)=idum_fluct_dn_dn(is,is)
  !         var_dd_dd(is)=idum_fluct_dd_dd(is,is)
  !      end forall
  !      cov_dn_dd=cov_dn_dD + idum_fluct_dn_dD/real(list_idum_length,8)
  !      cov_dn_dn=cov_dn_dn + idum_fluct_dn_dn/real(list_idum_length,8)
  !      cov_dD_dD=cov_dD_dD + idum_fluct_dD_dD/real(list_idum_length,8)
  !   enddo
  !   forall(is=1:Ns)
  !      var_dn_dd(is)=cov_dn_dd(is,is)
  !      var_dn_dn(is)=cov_dn_dn(is,is)
  !      var_dd_dd(is)=cov_dd_dd(is,is)
  !   end forall
  !   forall(is=1:Ns)r_local(is) = var_dn_dd(is)/sqrt(var_dn_dn(is)*var_dd_dd(is))
  !   call splot("covariance_dn_dd.local.data",var_dn_dd)
  !   call splot("correlation_dn_dd.local.data",r_local)
  ! end subroutine get_covariance

  !*************************************************************************************
  !*************************************************************************************
  !*************************************************************************************



  ! subroutine get_correlations()
  !   character(len=20)              :: idum_dir
  !   real(8),dimension(Ns)          :: nbar,dbar
  !   real(8),dimension(Ns)          :: ni,di,den,ded
  !   integer                        :: i,j,is,js,ird,idum
  !   real(8),dimension(Ns,Ns)       :: cov_dn_dD,cov_dn_dn,cov_dD_dD
  !   real(8),dimension(Ns)          :: var_dn_dd,var_dn_dn,var_dd_dd,r_local
  !   real(8),dimension(Ns,Ns)       :: idum_fluct_dn_dD,idum_fluct_dn_dn,idum_fluct_dD_dD
  !   cov_dn_dD=0.d0 ; cov_dn_dn=0.d0 ; cov_dD_dD=0.d0
  !   !loop over all realizations
  !   do ird=1,list_idum_length
  !      idum=list_idum(ird)
  !      write(idum_dir,"(I12)")idum
  !      idum_dir="idum_"//trim(adjustl(trim(idum_dir)))
  !      call sread(trim(adjustl(trim(idum_dir)))//"/nVSisite.ipt",ni)
  !      call sread(trim(adjustl(trim(idum_dir)))//"/deltaVSisite.ipt",Di)
  !      do is=1,Ns
  !         do js=1,Ns
  !            idum_fluct_dn_dD(is,js)=ni(is)*di(js)
  !            idum_fluct_dn_dn(is,js)=ni(is)*ni(js)
  !            idum_fluct_dD_dD(is,js)=di(is)*di(js)
  !         enddo
  !      enddo
  !      cov_dn_dd=cov_dn_dD + idum_fluct_dn_dD/real(list_idum_length,8)
  !      cov_dn_dn=cov_dn_dn + idum_fluct_dn_dn/real(list_idum_length,8)
  !      cov_dD_dD=cov_dD_dD + idum_fluct_dD_dD/real(list_idum_length,8)
  !   enddo
  !   forall(is=1:Ns,js=1:Ns)
  !      cov_dn_dd(is,js)=cov_dn_dd(is,js)/cov_dn_dd(is,is)
  !      cov_dn_dn(is,js)=cov_dn_dn(is,js)/cov_dn_dn(is,is)
  !      cov_dd_dd(is,js)=cov_dd_dd(is,js)/cov_dd_dd(is,is)
  !   end forall
  !   call splot("correlation_ninj.data",cov_dn_dn(:,1))
  !   call splot("correlation_didj.data",cov_dd_dd(:,1))
  !   call splot("correlation_nidj.data",cov_dn_dd(:,1))
  ! end subroutine get_correlations



  !****************************************************************************
  !****************************************************************************
  !**************************************************************************** 


end program get_data
