module parallelized_cg_method
implicit none
private
public parallelized_inv1
contains

subroutine parallelized_inv1(Eig,x,restart)
use parameters, only : ldim
use working_area, only : wwk
implicit none
complex(KIND(0d0)), intent(inout) :: x(:)
double precision, intent(inout) :: Eig 
logical, intent(in) :: restart
call parallelized_inv1z(restart,Eig, x, wwk(:,1), wwk(:,2), wwk(:,3), wwk(:,4))
end subroutine parallelized_inv1
!
!****************** inverse iteration ***********************
subroutine parallelized_inv1z(restart,Eig,x,b,p,r,y)
use parameters, only : ldim, max_ldim, progress, itr_counter,max_comp_time
use parallelized_diag, only : calc_vec_inprdct, parallelized_set_initial_vec_p
use openmp_lib
use parallel_lib
use def_output, only : interval_dump_read,interval_dump
implicit none
logical, intent(in) :: restart 
double precision, intent(inout) :: Eig
complex(KIND(0d0)), intent(inout) :: x(:),b(:),r(:),y(:),p(:)
!!
logical :: loop=.false.
integer :: i, itr, iterat, start
double precision :: xnorm, xb, eigen(1)
complex(KIND(0d0)) :: xcnorm, cxb
!!
if(itr_counter(7)<2)then
  if(.not.restart)then
    do i=1,ldim
      b(i)=dcmplx(0d0,0d0)
    enddo
    call parallelized_set_initial_vec_p(b)
  else
    do i=1,ldim
      b(i)=x(i)
    enddo
  endif
else
  call interval_dump_read(position=progress, data_size=max_ldim, se=eigen, &
                          vec_size=max_ldim, vec0=x,vec1=b,vec2=r,vec3=p)
  Eig=eigen(1)
endif
!!
start=itr_counter(5)
do itr=start,itr_counter(6)
  call barrier
  call parallelized_cg1(Eig,x,b,p,r,y,iterat)
  itr_counter(7)=1
  call calc_vec_inprdct(ldim,x,x,xcnorm)
  xnorm=1.d0/dsqrt(DREAL(xcnorm))
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(xnorm)
!$OMP do 
    do i=1,ldim
      x(i)=x(i)*xnorm
    enddo
!$OMP END DO
!$OMP END PARALLEL
!  if(iterat>ldim)then
!          print *,' #(W10)# Iterat in cg1 exceeds ldim or 5000'
!          print *,'         Approximate eigenvector returned'
!          print *,'         Itration number in inv1 is',itr
!          progress(3)=.false.
!          call barrier
!          return
!  endif
!!
  call calc_vec_inprdct(ldim,x,b,cxb)
!!
  if(lme==0) write(*,*) itr,dabs((CDABS(cxb))-1d0)
  if(dabs(CDABS(cxb)-1.0d0)<1.0d-13) loop=.true.
  call mpi_check_logical(loop)
  if(loop)then
          if(lme==0) write(*,100)itr
 100      format('       number of iterations in inv1 :',i5)
    progress(3)=.false.
    call barrier
    return
  endif
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED)
!$OMP do 
  do i=1,ldim
    b(i)=x(i)
  enddo
!$OMP END DO
!$OMP END PARALLEL
itr_counter(5)=itr_counter(5)+1
enddo
write(23,*)' #(W11)# inv1 did not converge'
end subroutine parallelized_inv1z
!
!************** solution of linear equations -- cg method ************
!
subroutine parallelized_cg1(Eig,x,b,p,r,y,itr)
use parameters, only : ldim, max_ldim, ibond, itr_counter, max_comp_time, progress
use parallelized_diag, only : p_mltply_v3, calc_vec_inprdct
use parallel_lib
use openmp_lib
use timer
use def_output, only : interval_dump
implicit none
double precision, intent(in) :: Eig
complex(KIND(0d0)), intent(inout) :: b(:),x(:),r(:),y(:),p(:)
integer, intent(inout) :: itr
complex(KIND(0d0)) :: crp, cyp, calpha, bcnorm, rcnorm, rcnorm2, cdmm
double precision :: eigen(1), eperbd, beta
integer :: i, start
logical :: loop=.false.
!*** initialization
call barrier
eigen(1)=Eig
call calc_vec_inprdct(ldim,b,b,bcnorm)
if(itr_counter(7)<2)then
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED)
!$OMP do 
  do i=1,ldim
    r(i)=b(i)
    p(i)=b(i)
    x(i)=dcmplx(0d0,0d0)
  enddo
!$OMP END DO
!$OMP END PARALLEL
endif
!!
    eperbd=Eig/dble(ibond)
!!
!*** iteration
start=itr_counter(7)
do itr=start,itr_counter(8)
  call barrier
  if(tell_me_time()>max_comp_time)then
    call interval_dump(position=progress, data_size=max_ldim, se=eigen, &
                          vec_size=max_ldim, vec0=x,vec1=b,vec2=r,vec3=p)
  endif  
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED)
!$OMP do 
  do i=1,ldim
    y(i)=dcmplx(0d0,0d0)
  enddo
!$OMP END DO
!$OMP END PARALLEL
  call p_mltply_v3(p,y,cdmm,eperbd) 
  crp=dcmplx(0d0,0d0)
  cyp=dcmplx(0d0,0d0)
  call calc_vec_inprdct(ldim,p,r,crp)
  call calc_vec_inprdct(ldim,p,y,cyp)
!!
  calpha=crp/cyp
  rcnorm=dcmplx(0d0,0d0)
  rcnorm2=dcmplx(0d0,0d0)
  call calc_vec_inprdct(ldim,r,r,rcnorm)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(calpha) 
!$OMP do 
  do i=1,ldim
    x(i)=x(i)+calpha*p(i)
    r(i)=r(i)-calpha*y(i)
  enddo
!$OMP END DO
!$OMP END PARALLEL
  call calc_vec_inprdct(ldim,r,r,rcnorm2)
!!!!!!!!!!
  beta=DREAL(rcnorm2)/DREAL(rcnorm)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(beta)
!$OMP do 
  do i=1,ldim
    p(i)=r(i)+beta*p(i)
  enddo
!$OMP END DO
!$OMP END PARALLEL
  itr_counter(7)=itr_counter(7)+1
!!  if(mod(itr,5)==0) write(*,'(A,I6,A,1pe16.6e3,e16.6e3)') 'ITR :',itr,' ',rcnorm2
  if(mod(itr,5).ne.0) cycle
  if( dsqrt(CDABS(rcnorm2))<1.0d-9*dsqrt(CDABS(bcnorm)) ) loop=.true.
  call barrier
  call mpi_check_logical(loop)
  if(loop)then
          if(lme==0)write(*,150) itr
 150      format('       number of iterations in cg1     :',i5)
          call barrier
          return
  end if
enddo
write(23,*)' #(Wxx)# cg1 did not converge'
end subroutine parallelized_cg1
!
end module parallelized_cg_method
