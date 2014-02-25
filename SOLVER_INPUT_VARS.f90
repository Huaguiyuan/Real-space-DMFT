module SOLVER_INPUT_VARS
  USE COMMON_VARS
  USE PARSE_CMD
  implicit none

  !SIZE OF THE ED PROBLEM
  !Nbath=# of bath sites (per orbital or not depending on bath_type)
  !Norb =# of impurity orbitals
  !Nspin=# spin degeneracy (max 2)
  !=========================================================
  integer                                     :: Nbath
  integer                                     :: Norb,Nspin

  !Global variables
  !=========================================================
  integer              :: L              !a large number fix the size of the problem
  real(8)              :: ts             !hopping amplitude
  integer              :: nloop          !max dmft loop variables
  real(8)              :: u              !local  interaction
  real(8)              :: dt,dtau        !time step
  real(8)              :: fmesh          !freq. step
  real(8)              :: wmin,wmax
  real(8)              :: weight
  real(8)              :: Ust,Jh         !intra-orbitals interactions
  real(8),dimension(3) :: Uloc           !local interactions
  real(8)              :: xmu            !chemical potential
  real(8)              :: deltasc        !breaking symmetry field
  real(8)              :: beta           !inverse temperature
  real(8)              :: eps            !broadening
  real(8)              :: wini,wfin      !
  integer              :: Nsuccess       !
  logical              :: Jhflag         !spin-exchange and pair-hopping flag.
  logical              :: chiflag        !
  logical              :: HFmode         !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)              :: cutoff         !cutoff for spectral summation
  real(8)              :: dmft_error     !dmft convergence threshold
  integer              :: lanc_niter     !Max number of Lanczos iterations
  integer              :: lanc_neigen    !Max number of required eigenvalues per sector
  integer              :: lanc_ngfiter   !Max number of iteration in resolvant tri-diagonalization
  integer              :: lanc_nstates   !Max number of states hold in the finite T calculation
  integer              :: cg_Niter       !Max number of iteration in the fit
  real(8)              :: cg_Ftol        !Tolerance in the cg fit
  integer              :: cg_Weight      !CGfit mode 0=normal,1=1/n weight, 2=1/w weight
  character(len=5)     :: cg_Scheme      !fit scheme: delta (default), weiss for G0^
  logical              :: finiteT        !flag for finite temperature calculation
  character(len=4)     :: ed_method      !flag to set ed method solution: lanc=lanczos method, full=full diagonalization
  character(len=1)     :: ed_type        !flag to set real or complex Ham: d=symmetric H (real), c=hermitian H (cmplx)
  logical              :: ed_supercond   !flag to set ed symmetry type: F=normal (default), T=superc=superconductive
  character(len=7)     :: bath_type      !flag to set bath type: irreducible (1bath/imp), reducible(1bath)
  real(8)              :: nread          !fixed density. if 0.d0 fixed chemical potential calculation.
  real(8)              :: nerr           !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)              :: ndelta         !initial chemical potential step
  integer              :: niter



  !Dimension of the functions:
  !=========================================================
  integer              :: NL,Nw,Nfit,Ltau


  !LOG AND Hamiltonian UNITS
  !=========================================================
  character(len=32)    :: Hfile
  integer              :: LOGfile



  !Namelists:
  !=========================================================
  namelist/solver_vars/&
       beta,ts,xmu,nloop,eps,weight,deltasc,&
       Norb,Nbath,Nspin,U,Uloc,Ust,Jh,& 
       wmax,wini,wfin,      &
       L,NL,Nw,Nfit,Ltau,         &
       nread,nerr,ndelta,       &
       chiflag,Jhflag,cutoff,HFmode,   &
       dmft_error,Nsuccess,      &
       ed_method,ed_type,bath_type,&
       lanc_neigen,lanc_niter,lanc_ngfiter,lanc_nstates,&
       cg_niter,cg_ftol,cg_weight,cg_scheme,  &
       Hfile,LOGfile

contains

  !+----------------------------------------------------------------+
  !PURPOSE  : Read SOLVER input file
  !+----------------------------------------------------------------+
  subroutine solver_read_input(inputFILE)
    character(len=*) :: inputFILE
    integer :: i,iorb
    logical :: control
    !DEFAULT VALUES OF THE PARAMETERS:
    !ModelConf
    Norb       = 1
    Nbath      = 4
    Nspin      = 1
    ts         = 0.5d0
    Uloc       = 0.d0;Uloc(1)=2.d0
    U          = Uloc(1)
    Ust        = 0.d0
    Jh         = 0.d0
    xmu        = 0.d0
    deltasc    = 2.d-2
    beta       = 500.d0
    weight     = 0.9d0
    !Loops
    nloop      = 100
    chiflag    =.true.
    Jhflag     =.false.
    hfmode     =.true.
    !parameters
    L          = 2048
    NL         = L
    Nw         = L
    Ltau       = 1000
    Nfit       = 1000
    eps        = 0.01d0
    nread      = 0.d0
    nerr       = 1.d-4
    ndelta     = 0.1d0
    wmax       = 5.d0
    wini       =-4.d0
    wfin       = 4.d0
    cutoff     = 1.d-9
    dmft_error  = 1.d-5
    nsuccess   = 2
    lanc_niter = 512
    lanc_neigen = 1
    lanc_ngfiter = 100
    lanc_nstates = 1            !set to T=0 calculation
    cg_niter   = 200
    cg_Ftol     = 1.d-9
    cg_weight     = 0
    cg_scheme     = 'delta'
    ed_method    = 'lanc'
    ed_type = 'd'
    ed_supercond = .false.
    bath_type='normal' !hybrid,superc
    !ReadUnits
    Hfile  ="hamiltonian.restart"
    LOGfile=6



    inquire(file=inputFILE,exist=control)    
    if(control)then
       open(10,file=inputFILE,status='old')
       read(10,nml=solver_vars)
       close(10)
    else
       if(mpiID==0)then
          print*,"Can not find INPUT file"
          print*,"Printing a default version in default."//inputFILE
          open(50,file="default."//inputFILE)
          write(50,nml=solver_vars)
          write(50,*)""
       endif
       stop
    endif

    call parse_cmd_variable(Norb,"NORB")    
    call parse_cmd_variable(Nbath,"NBATH")
    call parse_cmd_variable(Nspin,"NSPIN")
    call parse_cmd_variable(beta,"BETA")
    call parse_cmd_variable(xmu,"XMU")
    call parse_cmd_variable(U,"U")
    call parse_cmd_variable(uloc,"ULOC")
    call parse_cmd_variable(uloc(1),"U")
    call parse_cmd_variable(ust,"UST")
    call parse_cmd_variable(Jh,"JH")
    call parse_cmd_variable(nloop,"NLOOP")
    call parse_cmd_variable(deltasc,"DELTASC")
    call parse_cmd_variable(dmft_error,"DMFT_ERROR")
    call parse_cmd_variable(nsuccess,"NSUCCESS")
    call parse_cmd_variable(ts,"TS")
    call parse_cmd_variable(L,"L")
    call parse_cmd_variable(wmax,"WMAX")
    call parse_cmd_variable(weight,"WEIGHT")
    call parse_cmd_variable(L,"NL")
    call parse_cmd_variable(L,"NW")
    call parse_cmd_variable(Ltau,"LTAU")
    call parse_cmd_variable(Nfit,"NFIT")
    call parse_cmd_variable(nread,"NREAD")
    call parse_cmd_variable(nerr,"NERR")
    call parse_cmd_variable(ndelta,"NDELTA")
    call parse_cmd_variable(wini,"WINI")
    call parse_cmd_variable(wfin,"WFIN")
    call parse_cmd_variable(chiflag,"CHIFLAG")
    call parse_cmd_variable(hfmode,"HFMODE")
    call parse_cmd_variable(eps,"EPS")
    call parse_cmd_variable(cutoff,"CUTOFF")
    call parse_cmd_variable(lanc_neigen,"LANC_NEIGEN")
    call parse_cmd_variable(lanc_niter,"LANC_NITER")
    call parse_cmd_variable(lanc_nstates,"LANC_NSTATES")
    call parse_cmd_variable(lanc_ngfiter,"LANC_NGFITER")
    call parse_cmd_variable(cg_niter,"CG_NITER")
    call parse_cmd_variable(cg_scheme,"CG_SCHEME")
    call parse_cmd_variable(cg_ftol,"CG_FTOL")
    call parse_cmd_variable(cg_weight,"CG_WEIGHT")
    call parse_cmd_variable(ed_Type,"ED_TYPE")
    call parse_cmd_variable(ed_Supercond,"ED_SUPERCOND")
    call parse_cmd_variable(ed_Method,"ED_METHOD")
    call parse_cmd_variable(bath_type,"BATH_TYPE")
    call parse_cmd_variable(Hfile,"HFILE")
    call parse_cmd_variable(LOGfile,"LOGFILE")
    Ltau=max(int(beta),Ltau)
    NL=L
    Nw=L
    !


    if(mpiID==0)then
       open(50,file="used."//inputFILE)
       write(50,nml=solver_vars)
       close(50)
       write(LOGfile,*)"SOLVER CONTROL PARAMETERS"
       write(LOGfile,nml=solver_vars)
       write(LOGfile,"(A)")"U_local:"
       write(LOGfile,"(90F12.6,1x)")(Uloc(iorb),iorb=1,Norb)
    endif

    return
  end subroutine solver_read_input



end module SOLVER_INPUT_VARS
