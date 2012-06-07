  !     START MPI:
  !=====================================================================
  call start_mpi()

  !     READ INPUT FILES:
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
  allocate(etrap(Ns))
  allocate(H0(Ns,Ns))
  allocate(icol(Ns),irow(Ns))
  allocate(ij2site(Nside,Nside))


  !     CREATE DATA_DIR FOR AND STORE THIS IDUM 
  name_dir="data_trap"
  if(mpiID==0)then
     call create_data_dir(trim(adjustl(trim(name_dir))))
  endif


  !     BUILD THE LATTICE HAMILTONIAN:
  !=====================================================================
  call get_tb_hamiltonian


  !     BUILD THE TRAP:
  !=====================================================================
  do is=1,Ns
     etrap(is)= 0.5d0*V0trap*trap_distance_square(is) + a0trap
  enddo
  if(mpiID==0)call splot("fort.53",etrap)