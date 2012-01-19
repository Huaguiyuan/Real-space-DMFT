program get_data
  USE RDMFT_VARS_GLOBAL
  USE DMFT_IPT
  USE CHRONOBAR
  USE VECTORS
  USE LATTICE
  USE IOTOOLS
  USE SLPLOT
  USE SLREAD
  USE TOOLS
  USE GRIDS
  USE MATRIX,    ONLY:mat_inversion_sym,mat_inversion_gj
  USE RANDOM,    ONLY:nrand,drand,init_random_number
  implicit none
  logical :: check,logic
  integer :: list_idum_length
  integer :: is,iread,row,col
  real(8) :: ran
  character(len=16)      :: DIR
  real(8),allocatable    :: nii(:),dii(:)
  real(8),allocatable    :: grid_x(:),grid_y(:)
  real(8),allocatable    :: dij(:,:),nij(:,:),eij(:,:)
  real(8),allocatable    :: dijbar(:,:),nijbar(:,:)
  real(8),allocatable    :: rDOS(:)
  real(8),allocatable    :: arho_a(:,:,:),arho_g(:,:,:)
  real(8),allocatable    :: aarho_a(:,:),aarho_g(:,:)
  complex(8),allocatable,dimension(:,:,:) :: fgsc,fg0sc,sigmasc,gfsc_tmp
  real(8)                                 :: r,n

  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)

  !SET OMP THREADS NUMBER (THIS IS USUALLY 1 AS YOU DO NOT WANT TO HAVE
  !OMP PARALLEL HERE)
  call omp_set_num_threads(OMP_NUM_THREADS)

  !Read input file and allocate necessary memory:
  !distributes realizations among processors
  call read_input("inputIPT.in")
  call rdmft_read_input("inputRDMFT.in")
  Ns   =Nside**2
  n    =0.5d0
  allocate(erandom(Ns))
  allocate(H0(Ns,Ns))
  allocate(icol(Ns),irow(Ns))
  allocate(ij2site(Nside,Nside))
  allocate(fgsc(2,Ns,-L:L))
  allocate(fg0sc(2,Ns,-L:L))
  allocate(sigmasc(2,Ns,-L:L))
  allocate(gfsc_tmp(2,Ns,-L:L))
  allocate(nii(Ns),dii(Ns))
  allocate(nij(Nside,Nside),dij(Nside,Nside),eij(Nside,Nside))
  allocate(grid_x(Nside),grid_y(Nside))
  allocate(nijbar(Nside,Nside),dijbar(Nside,Nside))
  allocate(arho_a(2,Ns,-L:L),arho_g(2,Ns,-L:L))
  allocate(aarho_a(2,-L:L),aarho_g(2,-L:L))
  allocate(rDOS(-L:L))

  nijbar=0.d0
  dijbar=0.d0

  !Create directory for the averaged results:
  DIR="RESULTS" ; call create_data_dir(DIR)

  !Build the Lattice Hamiltonian:
  call get_tb_hamiltonian
  store_size=1024

  !Build grid for the lattice plots:
  do row=1,Nside
     grid_x(row)=dble(row)
     grid_y(row)=dble(row)
  enddo

  !Inquire if list_idum exist && check length:
  inquire(file="list_idum",exist=check)
  if(.not.check)call abort("I can not find list_idum! ciao")
  list_idum_length=file_length("list_idum")

  call init_random_number

  !Start loop over disorder realizations:
  open(100,file="list_idum")
  arho_a=0.d0 ; aarho_a=0.d0
  arho_g=0.d0 ; aarho_g=0.d0
  do iread=1,list_idum_length
     !Read idum:
     read(100,*)idum; print*,"read idum:",idum
     write(name_dir,"(I12)")idum
     name_dir="idum_"//trim(adjustl(trim(name_dir)))

     !Read idum-realized q.ties: [n,delta,Sigma,Self]
     call read_all(trim(adjustl(trim(name_dir))))

     ran=drand()
     logic=.false.;if(ran>0.5d0)logic=.true.
     logic=.true.

     do i=1,Ns
        call splot("LDOS.site.disorder.ipt",wr,-dimag(fgsc(1,i,-L:L))/pi,append=TT)
     enddo
     call system("mv -vf LDOS.site.disorder.ipt "//trim(adjustl(trim(name_dir)))//"/ 2>/dev/null")

     !Get avaraged q.ties:
     aarho_a(1,:) = aarho_a(1,:) + rDOS(:)/dble(list_idum_length)
     aarho_g(1,:) = aarho_g(1,:) + log(rDOS(:))/dble(list_idum_length)

     !site_x, site_y, n/delta/e(site_x,site_y)
     if(logic)then
        call splot("3d_nVSij_"//trim(adjustl(trim(name_dir)))//".ipt",grid_x,grid_y,nij)
        call splot("3d_deltaVSij_"//trim(adjustl(trim(name_dir)))//".ipt",grid_x,grid_y,dij)
        call splot("3d_erandomVSij_"//trim(adjustl(trim(name_dir)))//".ipt",grid_x,grid_y,eij)
        call system("mv -vf plot_3d* 3d_*"//trim(adjustl(trim(name_dir)))//"* "//trim(adjustl(trim(name_dir)))//"/ 2>/dev/null")
     endif
  enddo

  !***********************
  !PLOT REMAINING Q.TIES
  !***********************
  ! aarho_g(1,:) = exp(aarho_g(1,:))
  ! aarho_a(1:2,-L:L)   =sum(arho_a(1:2,1:Ns,-L:L),dim=2)/dble(Ns)
  ! aarho_g(1:2,-L:L)   =sum(arho_g(1:2,1:Ns,-L:L),dim=2)/dble(Ns)

  call splot("avDOSa.ipt",wr,aarho_a(1,-L:L))
  call splot("avDOSg.ipt",wr,exp(aarho_g(1,-L:L)))

  call system("mv -vf *.ipt "//trim(adjustl(trim(DIR)))//"/ 2>/dev/null")
  call MPI_FINALIZE(mpiERR)

contains

  subroutine read_all(idum_dir)
    character(len=*) :: idum_dir
    real(8)          :: w
    call sread(idum_dir//"/nVSisite.ipt",nii)
    call sread(idum_dir//"/deltaVSisite.ipt",dii)
    call sread(idum_dir//"/erandomVSisite.ipt",erandom)
    call sread(idum_dir//"/LSigma.ipt",sigmasc(1,1:Ns,-L:L))
    call sread(idum_dir//"/LSelf.ipt",sigmasc(2,1:Ns,-L:L))
    call sread(idum_dir//"/LG.ipt",fgsc(1,1:Ns,-L:L))
    call sread(idum_dir//"/LF.ipt",fgsc(2,1:Ns,-L:L))
    do row=1,Nside
       do col=1,Nside
          i=ij2site(row,col)
          nij(row,col)=nii(i)
          dij(row,col)=dii(i)
          eij(row,col)=erandom(i)
       enddo
    enddo
    inquire(file=idum_dir//"/DOS.disorder.ipt",exist=check)
    open(20,file=idum_dir//"/DOS.disorder.ipt")
    if(check)then
       do i=-L,L
          read(20,*)w,rDOS(i)
       enddo
    endif
    close(20)
    where(rDOS(:)<0.d0)rDOS=1.d-3
  end subroutine read_all



  include "solve_routines_sc.f90"
end program get_data
