program pade_
  USE COMMON_VARS
  USE TOOLS
  USE IOTOOLS
  USE PADE
  USE RDMFT_VARS_GLOBAL
  implicit none
  integer                :: i,j,is,npade
  !Matsubara:
  complex(8),allocatable,dimension(:,:,:) :: fg
  !Real-axis:
  complex(8),allocatable :: zr(:),gr(:,:)

  !READ INPUT FILES:
  !=====================================================================
  call read_input("inputIPT.in")
  call rdmft_read_input("inputRDMFT.in")
  call parse_cmd_variable(npade,"NPADE",default=10)

  Ns=Nside**2
  Ns=32
  write(name_dir,"(I12)")idum
  name_dir="idum_"//trim(adjustl(trim(name_dir)))


  allocate(wm(L),zr(L))
  allocate(fg(2,Ns,L),gr(2,L))

  wm(:)  = pi/beta*real(2*arange(1,L)-1,8)
  call sread(trim(adjustl(trim(name_dir)))//"/LG_iw.data",fg(1,1:Ns,1:L),wm)
  call sread(trim(adjustl(trim(name_dir)))//"/LF_iw.data",fg(2,1:Ns,1:L),wm)

  zr = linspace(-wmax,wmax,L)+xi*eps
  gr = zero
  call system("rm -fv "//trim(adjustl(trim(name_dir)))//"/L*_realw.pade")
  do is=1,Ns
     gr(1,:) = pade_analytic_continuation(zr,npade,fg(1,is,:),wm(:))
     gr(2,:) = pade_analytic_continuation(zr,npade,fg(2,is,:),wm(:))
     call splot(trim(adjustl(trim(name_dir)))//"/LG_realw.pade",real(zr(:),8),gr(1,:),append=.true.)
     call splot(trim(adjustl(trim(name_dir)))//"/LF_realw.pade",real(zr(:),8),gr(2,:),append=.true.)
  enddo
end program pade_
