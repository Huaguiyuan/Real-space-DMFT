  !     START MPI:
  !=====================================================================
  call start_mpi()

  !     READ INPUT FILES:
  !=====================================================================
  call read_input("inputIPT.in")
  call rdmft_read_input("inputRDMFT.in")
  !     lm  = int(real(L,8)*beta/pi2)
  !     esp = nint(log(real(lm,8))/log(2.d0))
  !     lm  = 2**esp
  !     L=max(lm,L)
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


  !     CREATE DATA_DIR FOR AND STORE THIS IDUM 
  write(name_dir,"(I12)")idum
  name_dir="idum_"//trim(adjustl(trim(name_dir)))
  if(mpiID==0)then
     call create_data_dir(trim(adjustl(trim(name_dir))))
     open(10,file="list_idum",access='append')
     write(10,*)idum
     close(10)
  endif


  !     BUILD THE LATTICE HAMILTONIAN:
  !=====================================================================
  call get_tb_hamiltonian


  !     BUILD RANDOM ENERGIES:
  !     GET RID OF THE FIRST FEW NON-RANDOM NUMBERS IN NR RAN1 GENERATOR
  !=====================================================================
  do i=1,50
     r=nrand(idum)
  enddo
  do is=1,Ns
     erandom(is)=(2.d0*nrand(idum)-1.d0)*Wdis/2.d0
  enddo

