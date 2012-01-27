  !START MPI:
  !=====================================================================
  call start_mpi()

  !READ INPUT FILES:
  !=====================================================================
  call read_input("inputIPT.in")
  call rdmft_read_input("inputRDMFT.in")


  !ALLOCATE WORKING ARRAYS:
  !=====================================================================
  Ns    =Nside**2
  wmax  =wmax+Wdis
  allocate(erandom(Ns))
  allocate(H0(Ns,Ns))
  allocate(icol(Ns),irow(Ns))
  allocate(ij2site(Nside,Nside))


  !CREATE DATA_DIR FOR AND STORE THIS IDUM 
  !=====================================================================
  write(name_dir,"(I12)")idum
  name_dir="idum_"//trim(adjustl(trim(name_dir)))
  if(mpiID==0)then
     call create_data_dir(trim(adjustl(trim(name_dir))))
     open(10,file="list_idum",access='append')
     write(10,*)idum
     close(10)
  endif


  !BUILD THE LATTICE HAMILTONIAN:
  !=====================================================================
  call get_tb_hamiltonian


  !BUILD RANDOM ENERGIES:
  !=====================================================================
  do i=1,50                     !get rid of few spurious random number in NR
     r=nrand(idum)
  enddo
  do is=1,Ns
     erandom(is)=(2.d0*nrand(idum)-1.d0)*Wdis/2.d0
  enddo

