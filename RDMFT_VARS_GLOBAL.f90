!###############################################################
!     PROGRAM  : VARS_GLOBAL
!     TYPE     : Module
!     PURPOSE  : Contains global variables
!     AUTHORS  : Adriano Amaricci
!###############################################################
module RDMFT_VARS_GLOBAL
  !LOAD MODULE FROM LIBRARY:
  USE COMMON_VARS
  USE CHRONOBAR, ONLY:start_timer,stop_timer,eta
  USE IOTOOLS
  USE SLPLOT
  USE SLREAD
  USE MATRIX,    ONLY:mat_inversion_sym,mat_inversion_gj,mat_inversion
  USE RANDOM,    ONLY:nrand
  USE INTEGRATE, ONLY:kronig
  USE TOOLS,     ONLY:fermi,check_convergence
  USE DMFT_IPT
  USE MPI
  USE OMP_LIB
  implicit none

  !Lattice size:
  !=========================================================
  integer,protected :: Nside
  integer           :: Ns

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau
  real(8),dimension(:),allocatable :: wr,t


  !Large matrices for Lattice Hamiltonian/GF
  !=========================================================
  integer,dimension(:),allocatable   :: icol,irow
  integer,dimension(:,:),allocatable :: ij2site
  real(8),dimension(:,:),allocatable :: H0,Id


  !Gloabl variables
  !=========================================================
  real(8) :: Wdis               !Disorder strength
  integer :: idum               !disorder seed
  real(8) :: a0trap,V0trap

  !Other variables:
  !=========================================================
  character(len=20)                  :: name_dir
  logical                            :: pbcflag

  !Random energies
  !=========================================================
  real(8),dimension(:),allocatable   :: erandom,etrap

  !Fix density
  !=========================================================
  real(8) :: nread,nerror
  real(8) :: ndelta,ndelta1
  integer :: nindex,nindex1

  !Namelists:
  !=========================================================
  namelist/disorder/&
       Wdis,     &
       V0trap,   &
       a0trap,   &
       Nside,    &
       idum,     &
       nread,    &
       nerror,   &
       ndelta,   &
       pbcflag,  &
       omp_num_threads



contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine rdmft_read_input(inputFILE)
    character(len=*) :: inputFILE
    integer :: i
    !local variables: default values
    Wdis            = 0.5d0
    Nside           = 10
    a0trap          = 0.d0
    v0trap          = 0.1d0
    nread           = 0.d0
    nerror          = 1.d-4
    ndelta          = 0.1d0
    omp_num_threads =1
    pbcflag         = .true.
    idum            =1234567

    !SET SIZE THRESHOLD FOR FILE ZIPPING:
    store_size=1024

    open(10,file=adjustl(trim(inputFILE)))
    read(10,nml=disorder)
    close(10)
    allocate(help_buffer(23))
    help_buffer=([&
         'NAME',&
         '  xxx_disorder/trap',&
         '',&
         'DESCRIPTION',&
         '  This a layer interface to different codes solving the Real-space Dynamical Mean Field Theory.',&
         '  Each code must be implemented separetely, for convenience some driver routines are in drivers/ dir.',&
         '  The structure is very easy: an thight-binding hamiltonian is generated at fist and the DMFT problem  ',&
         '  is solved in the local wannier basis, solving the impurity problem at each lattice site. Simmetries ',&
         '  can be used to reduce the size of the problem. Lattice sites are related via self-consistency. ',&
         '  The code is MPI parallel. ',&
         '  ',&
         'OPTIONS',&
         ' wdis=[0.5]    -- degree of local disorder.',&
         ' Nside=[10]    -- linear size of the cluster to be solved.',&
         ' a0trap=[0]    -- bottom of the trap. here kept separated from mu.',&
         ' v0trap=[0.1]  -- fix the parabolic shape of the trap.',&
         ' nread=[0.0]   -- density value for chemical potential search.',&
         ' ndelta=[0.1]  -- starting value for chemical potential shift.',&
         ' nerror=[1.d-4]-- max error in adjusting chemical potential. ',&
         ' pbcflag=[T]   -- periodic boundary conditions.',&
         ' idum=[1234567]-- initial seed for the random variable sample.',&
         ' omp_num_threads=[1] -- fix the number of threads in OMP environment.',&
         '  '])
    call parse_cmd_help(help_buffer)
    call parse_cmd_variable(wdis,"WDIS")
    call parse_cmd_variable(v0trap,"V0TRAP")
    call parse_cmd_variable(a0trap,"A0TRAP")
    call parse_cmd_variable(Nside,"NSIDE")
    call parse_cmd_variable(nread,"NREAD")
    call parse_cmd_variable(nerror,"NERROR")
    call parse_cmd_variable(pbcflag,"PBCFLAG")
    call parse_cmd_variable(omp_num_threads,"OMP_NUM_THREADS")
    call parse_cmd_variable(idum,"IDUM")

    !SET OMP THREADS NUMBER
    call omp_set_num_threads(OMP_NUM_THREADS)

    !Print on the screen used vars
    if(mpiID==0)then
       write(*,nml=disorder)
       open(10,file="used."//adjustl(trim(inputFILE)))
       write(10,nml=disorder)
       close(10)
    endif

  end subroutine rdmft_read_input
  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : Build tight-binding Hamiltonian
  !+----------------------------------------------------------------+
  subroutine get_tb_hamiltonian()
    integer :: i,jj,j,k,row,col,link(4)
    H0=0.d0
    do row=0,Nside-1
       do col=0,Nside-1
          i=col+row*Nside+1
          irow(i)=row+1
          icol(i)=col+1
          ij2site(row+1,col+1)=i

          if(pbcflag)then
             !HOPPING w/ PERIODIC BOUNDARY CONDITIONS
             link(1)= row*Nside+1              + mod(col+1,Nside)  ;
             link(3)= row*Nside+1              + (col-1)           ; if((col-1)<0)link(3)=(Nside+(col-1))+row*Nside+1
             link(2)= mod(row+1,Nside)*Nside+1 + col               ; 
             link(4)= (row-1)*Nside+1          + col               ; if((row-1)<0)link(4)=col+(Nside+(row-1))*Nside+1
          else
             !without PBC
             link(1)= row*Nside+1              + col+1   ; if((col+1)==Nside)link(1)=0
             link(3)= row*Nside+1              +(col-1)  ; if((col-1)<0)     link(3)=0
             link(2)= (row+1)*Nside+1 + col              ; if((row+1)==Nside)link(2)=0
             link(4)= (row-1)*Nside+1          + col     ; if((row-1)<0)     link(4)=0
          endif
          do jj=1,4
             if(link(jj)>0)H0(i,link(jj))=ts
          enddo
       enddo
    enddo
  end subroutine get_tb_hamiltonian
  !******************************************************************
  !******************************************************************
  !******************************************************************

end module RDMFT_VARS_GLOBAL
