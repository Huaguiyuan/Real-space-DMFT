!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : function
!PURPOSE  : 
!+-------------------------------------------------------------------+
function i0_check_convergence_function(Xnew,eps,N1,N2,id,file,index,total,oerr) result(convergence)
  integer,intent(in)       :: Xnew(:)
  real(8),intent(in)       :: eps
  real(8),optional         :: oerr
  integer,intent(in)       :: N1,N2
  integer,optional         :: id,index,total
  integer                  :: id_,index_,total_
  integer                  :: i,j,Msum
  logical                  :: convergence  
  real(8)                  :: error,err
  real(8)                  :: M,S
  integer,save,allocatable :: Xold(:,:)
  integer,save             :: success=0,check=1
  character(len=2)         :: label
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msum=size(Xnew)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msum))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0
     do i=1,Msum
        M=M + abs(Xnew(i)-Xold(index_,i))
        S=S + abs(Xnew(i))
     enddo
     err= M/S
     Xold(index_,:)=Xnew
     include "tools_write_error_file_dim0.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_test_convergence.f90"
     include "tools_print_error_msg_dim0.f90"
     if(present(oerr))oerr=err
     check=check+1
  endif
end function i0_check_convergence_function

function i1_check_convergence_function(Xnew,eps,N1,N2,id,file,index,total,strict,oerr) result(convergence)
  integer,intent(in)              :: Xnew(:,:)
  real(8),intent(in)              :: eps
  real(8),optional         :: oerr
  integer,intent(in)              :: N1,N2
  integer,optional                :: id,index,total
  integer                         :: id_,index_,total_
  integer                         :: i,j,Msize1,Msum
  logical                         :: convergence  
  real(8)                         :: error(2),err
  real(8),dimension(size(Xnew,1)) :: M,S,Verror
  integer,save,allocatable        :: Xold(:,:,:)
  integer,save                    :: success=0,check=1
  character(len=2)         :: label
  logical,optional :: strict
  logical          :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  strict_=.false.;if(present(strict))strict_=strict
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1) ; Msum=size(Xnew,2)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msum))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0 ; Verror=0.d0
     do i=1,Msum
        M=M + abs(Xnew(:,i)-Xold(index_,:,i))
        S=S + abs(Xnew(:,i))
     enddo
     Verror= M/S
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:,:)=Xnew
     include "tools_write_error_file_dim1.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_test_convergence.f90"
     include "tools_print_error_msg_dim1.f90"
     if(present(oerr))oerr=err
     check=check+1
  endif
end function i1_check_convergence_function

function i2_check_convergence_function(Xnew,eps,N1,N2,id,file,index,total,strict,oerr) result(convergence)
  integer,intent(in)                           :: Xnew(:,:,:)
  real(8),intent(in)                           :: eps
  real(8),optional         :: oerr
  integer,intent(in)                           :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2,Msum
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: M,S,Verror
  integer,save,allocatable                     :: Xold(:,:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  logical,optional :: strict
  logical          :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  strict_=.false.;if(present(strict))strict_=strict
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1);Msize2=size(Xnew,2);Msum=size(Xnew,3)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msize2,Msum))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0 ; Verror=0.d0
     do i=1,Msum
        M=M + abs(Xnew(:,:,i)-Xold(index_,:,:,i))
        S=S + abs(Xnew(:,:,i))
     enddo
     Verror= M/S
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:,:,:)=Xnew
     include "tools_write_error_file_dim2.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_test_convergence.f90"
     include "tools_print_error_msg_dim2.f90"
     if(present(oerr))oerr=err
     check=check+1
  endif
end function i2_check_convergence_function


!----------------------------------------------------------------------


function d0_check_convergence_function(Xnew,eps,N1,N2,id,file,index,total,oerr) result(convergence)
  real(8),intent(in)       :: Xnew(:)
  real(8),intent(in)       :: eps
  real(8),optional         :: oerr
  integer,intent(in)       :: N1,N2
  integer,optional         :: id,index,total
  integer                  :: id_,index_,total_
  integer                  :: i,j,Msum
  logical                  :: convergence  
  real(8)                  :: error,err
  real(8)                  :: M,S
  real(8),save,allocatable :: Xold(:,:)
  integer,save             :: success=0,check=1
  character(len=2)         :: label
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msum=size(Xnew)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msum))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0
     do i=1,Msum
        M=M + abs(Xnew(i)-Xold(index_,i))
        S=S + abs(Xnew(i))
     enddo
     err= M/S
     Xold(index_,:)=Xnew
     include "tools_write_error_file_dim0.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_test_convergence.f90"
     include "tools_print_error_msg_dim0.f90"
     if(present(oerr))oerr=err
     check=check+1
  endif
end function d0_check_convergence_function

function d1_check_convergence_function(Xnew,eps,N1,N2,id,file,index,total,strict,oerr) result(convergence)
  real(8),intent(in)              :: Xnew(:,:)
  real(8),intent(in)              :: eps
  real(8),optional         :: oerr
  integer,intent(in)              :: N1,N2
  integer,optional                :: id,index,total
  integer                         :: id_,index_,total_
  integer                         :: i,j,Msize1,Msum
  logical                         :: convergence  
  real(8)                         :: error(2),err
  real(8),dimension(size(Xnew,1)) :: M,S,Verror
  real(8),save,allocatable        :: Xold(:,:,:)
  integer,save                    :: success=0,check=1
  character(len=2)         :: label
  logical,optional :: strict
  logical          :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  strict_=.false.;if(present(strict))strict_=strict
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1) ; Msum=size(Xnew,2)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msum))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0 ; Verror=0.d0
     do i=1,Msum
        M=M + abs(Xnew(:,i)-Xold(index_,:,i))
        S=S + abs(Xnew(:,i))
     enddo
     Verror= M/S
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:,:)=Xnew
     include "tools_write_error_file_dim1.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_test_convergence.f90"
     include "tools_print_error_msg_dim1.f90"
     if(present(oerr))oerr=err
     check=check+1
  endif
end function d1_check_convergence_function

function d2_check_convergence_function(Xnew,eps,N1,N2,id,file,index,total,strict,oerr) result(convergence)
  real(8),intent(in)                           :: Xnew(:,:,:)
  real(8),intent(in)                           :: eps
  real(8),optional         :: oerr
  integer,intent(in)                           :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2,Msum
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: M,S,Verror
  real(8),save,allocatable                     :: Xold(:,:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  logical,optional :: strict
  logical          :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  strict_=.false.;if(present(strict))strict_=strict
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1);Msize2=size(Xnew,2);Msum=size(Xnew,3)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msize2,Msum))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0 ; Verror=0.d0
     do i=1,Msum
        M=M + abs(Xnew(:,:,i)-Xold(index_,:,:,i))
        S=S + abs(Xnew(:,:,i))
     enddo
     Verror= M/S
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:,:,:)=Xnew
     include "tools_write_error_file_dim2.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_test_convergence.f90"
     include "tools_print_error_msg_dim2.f90"
     if(present(oerr))oerr=err
     check=check+1
  endif
end function d2_check_convergence_function


!----------------------------------------------------------------------


function z0_check_convergence_function(Xnew,eps,N1,N2,id,file,index,total,oerr) result(convergence)
  complex(8),intent(in)       :: Xnew(:)
  real(8),intent(in)          :: eps
  real(8),optional         :: oerr
  integer,intent(in)          :: N1,N2
  integer,optional            :: id,index,total
  integer                     :: id_,index_,total_
  integer                     :: i,j,Msum
  logical                     :: convergence  
  real(8)                     :: error,err
  real(8)                     :: M,S
  complex(8),save,allocatable :: Xold(:,:)
  integer,save                :: success=0,check=1
  character(len=2)         :: label
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msum=size(Xnew)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msum))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0
     do i=1,Msum
        M=M + abs(Xnew(i)-Xold(index_,i))
        S=S + abs(Xnew(i))
     enddo
     err= M/S
     Xold(index_,:)=Xnew
     include "tools_write_error_file_dim0.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_test_convergence.f90"
     include "tools_print_error_msg_dim0.f90"
     if(present(oerr))oerr=err
     check=check+1
  endif
end function z0_check_convergence_function

function z1_check_convergence_function(Xnew,eps,N1,N2,id,file,index,total,strict,oerr) result(convergence)
  complex(8),intent(in)           :: Xnew(:,:)
  real(8),intent(in)              :: eps
  real(8),optional         :: oerr
  integer,intent(in)              :: N1,N2
  integer,optional                :: id,index,total
  integer                         :: id_,index_,total_
  integer                         :: i,j,Msize1,Msum
  logical                         :: convergence  
  real(8)                         :: error(2),err
  real(8),dimension(size(Xnew,1)) :: M,S,Verror
  complex(8),save,allocatable     :: Xold(:,:,:)
  integer,save                    :: success=0,check=1
  character(len=2)         :: label
  logical,optional :: strict
  logical          :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  strict_=.false.;if(present(strict))strict_=strict
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1) ; Msum=size(Xnew,2)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msum))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0 ; Verror=0.d0
     do i=1,Msum
        M=M + abs(Xnew(:,i)-Xold(index_,:,i))
        S=S + abs(Xnew(:,i))
     enddo
     Verror= M/S
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:,:)=Xnew
     include "tools_write_error_file_dim1.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_test_convergence.f90"
     include "tools_print_error_msg_dim1.f90"
     if(present(oerr))oerr=err
     check=check+1
  endif
end function z1_check_convergence_function

function z2_check_convergence_function(Xnew,eps,N1,N2,id,file,index,total,strict,oerr) result(convergence)
  complex(8),intent(in)                        :: Xnew(:,:,:)
  real(8),intent(in)                           :: eps
  real(8),optional         :: oerr
  integer,intent(in)                           :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2,Msum
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: M,S,Verror
  complex(8),save,allocatable                  :: Xold(:,:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  logical,optional :: strict
  logical          :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  strict_=.false.;if(present(strict))strict_=strict
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1);Msize2=size(Xnew,2);Msum=size(Xnew,3)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msize2,Msum))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0 ; Verror=0.d0
     do i=1,Msum
        M=M + abs(Xnew(:,:,i)-Xold(index_,:,:,i))
        S=S + abs(Xnew(:,:,i))
     enddo
     Verror= M/S
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:,:,:)=Xnew
     include "tools_write_error_file_dim2.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_test_convergence.f90"
     include "tools_print_error_msg_dim2.f90"
     if(present(oerr))oerr=err
     check=check+1
  endif
end function z2_check_convergence_function
