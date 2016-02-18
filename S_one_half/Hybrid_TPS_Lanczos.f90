module parallelized_lanczos_method
implicit none
private
public parallelized_lnc1, parallelized_lncv1
!!
contains
subroutine parallelized_lnc1(nvec2,itr,ene,x,restart)
use parameters, only : n, ldim, max_ldim, ibond
use working_area, only : wwk
use data_mem, only : ipair
use basic, only : datack
use parallelized_diag, only : parallelized_set_initial_vec_p
use parallel_lib
implicit none
logical, intent(in) :: restart
integer, intent(in) :: nvec2
integer, intent(inout) :: itr
complex(kind(0d0)), intent(inout) :: x(:)
double precision :: ene(:)
!!
if(nvec2<0)then
  print *,' #(E06)# Wrong value given to nvec in lnc1'
  print *,'         Only the eigenvalues are calculated'
  stop
endif
call datack(ipair,ibond,n)
!!
!*** initialization
!$OMP PARALLEL
!$OMP WORKSHARE
wwk(1:max_ldim,1:2)=dcmplx(0d0,0d0)
!$OMP END WORKSHARE
!$OMP END PARALLEL
if(.not.restart)then
  call parallelized_set_initial_vec_p(wwk(:,1)) 
else
!$OMP PARALLEL
!$OMP WORKSHARE
  wwk(1:ldim,1)=x(1:ldim)
!$OMP END WORKSHARE
!$OMP END PARALLEL
endif

call parallelized_lnc1z(nvec2,ene,itr,wwk(:,1),wwk(:,2))
end subroutine parallelized_lnc1
!
!************ eigenvalues by the Lanczos method
subroutine parallelized_lnc1z(nvec2,ene,itr,v1,v0)
use parameters, only : ldim, max_ldim, max_itr, itr_counter,progress,max_comp_time
use working_area, only : alpha, beta, coef
use parallelized_diag, only : p_mltply_v3, calc_beta, calc_alpha_beta
use parallel_lib
use basic
use openmp_lib
use def_output, only : interval_dump_read,interval_dump
use timer
implicit none
integer, intent(in) :: nvec2
double precision, intent(inout) :: ene(:)
complex(KIND(0d0)), intent(inout) :: v1(:),v0(:)
integer, intent(inout) :: itr
integer :: i, j, start
double precision :: alpha1, beta1, ebefor, eps
complex(KIND(0d0)) :: temp1, temp2
logical :: loop=.false.
integer :: t_max_itr
double precision,allocatable :: wk(:,:)
integer,allocatable :: iwk(:)
t_max_itr=itr_counter(2)
allocate(iwk(t_max_itr))
allocate(wk(t_max_itr,5))
!! initial setting
!$OMP PARALLEL
!$OMP WORKSHARE
wk(1:t_max_itr,1:5)=0d0
iwk(1:t_max_itr)=0
!$OMP END WORKSHARE
!$OMP END PARALLEL
!!
if(itr_counter(1)<3)then
!*** alpha(1) and beta(1)
  call calc_alpha_beta(ldim,alpha1,beta1,v1,v0)
  alpha(1)=alpha1
  beta(1) =beta1
  if(lme==0)write(*,'(A,F18.6)') 'ITR :      1',tell_me_time()
!!
else
  call interval_dump_read(position=progress, data_size=max_ldim, &
                          vec_size=t_max_itr, sa=alpha, sb=beta ,sc=coef, se=ene,&
                          vec1=v0,vec2=v1)
       alpha1=alpha(itr_counter(1)-1)
       beta1=beta(itr_counter(1)-1)
endif
!!
!*** iteration
start=itr_counter(1)
do i=start,itr_counter(2)
if(lme==0)  write(*,'(A,I8,A,F12.4)') 'Itr :',i,'Time: ',tell_me_time()
!  if(mod(i,100)==0) write(*,'(A,F12.4)') 'Time: ',tell_me_time()
  if(tell_me_time()>max_comp_time)then
    call interval_dump(position=progress, data_size=max_ldim, &
                       vec_size=t_max_itr, sa=alpha, sb=beta, sc=coef, se=ene,&
                       vec1=v0, vec2=v1)
  endif
!!
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(ldim,alpha1,beta1) &
!$OMP& PRIVATE(j,temp1,temp2) 
!$OMP do 
  do j=1,ldim
    temp1=v1(j)
    temp2=(v0(j)-alpha1*v1(j))/beta1
    v0(j)=-beta1*temp1
    v1(j)=temp2
  enddo
!$OMP END DO
!$OMP END PARALLEL
!!
  call calc_alpha_beta(ldim,alpha1,beta1,v1,v0)
  alpha(i)=alpha1
  beta(i) =beta1
!!write(*,*) 'ITR :',i,alpha1,beta1
!!
!!
  if(beta(i).lt.0.5d-30)then
            print *,' #(E07)# Tridiagonalization unsuccessful in lnc1'
            print *,'         Beta(i) is too small at i=  ',i
            stop
  end if
!!
!*** convergence check
  if( (i>4*nvec2).and.(mod(i,5)==0) )then
    call bisec(alpha,beta,i,ene,nvec2,eps)
    if(dabs((ebefor-ene(nvec2))/ene(nvec2))<1d-13) loop=.true.
    write(*,*) 'GS ENE:',ene(1)
    call mpi_check_logical(loop)
    if(loop)then
      if(nvec2>0)then 
        call vec12(alpha,beta,coef,ene,i,nvec2, &
             wk(:,1),wk(:,2),wk(:,3),wk(:,4),wk(:,5),iwk)
      end if
      itr=i
      progress(1)=.false.
      itr_counter(4)=itr
      deallocate(wk,iwk)
      return
    end if
    ebefor=ene(nvec2)
  end if
  if(i==4*nvec2)then
    eps=1d-25
    call bisec(alpha,beta,4*nvec2,ene,nvec2,eps)
    ebefor=ene(nvec2)
  end if
  itr_counter(1)=itr_counter(1)+1
enddo
!
write(*,'(A,I4,A)')' #(W07)# lnc1 did not converge within', itr-1,' steps'
itr=max(max_itr-1,ldim)
deallocate(wk,iwk)
!
end subroutine parallelized_lnc1z
!
!************ eigenvector by the Lanczos method *************
subroutine parallelized_lncv1(nvec2,ene,x,itr,restart)
use parameters, only : ldim
use working_area, only : wwk
implicit none
integer, intent(in) :: nvec2, itr
complex(KIND(0d0)), intent(inout) :: x(:) 
double precision, intent(inout) :: ene(:)
logical, intent(in) :: restart
!
if(nvec2<0)then
          print *,'#(W08)# nvec given to lncv1 out of range'
          return
end if
call parallelized_lncv1z(nvec2,ene,x,wwk(:,1),wwk(:,2),itr,restart)
end subroutine parallelized_lncv1
!
!************ eigenvector by the Lanczos method
subroutine parallelized_lncv1z(nvec2,ene,x,v1,v0,itr,restart)
use parameters, only : ldim, max_ldim, itr_counter, progress,max_comp_time
use working_area, only : alpha, beta, coef
use parallelized_diag, only : parallelized_set_initial_vec_p, p_mltply_v3, calc_vec_inprdct
use openmp_lib
use def_output, only : interval_dump_read,interval_dump
use timer
implicit none 
!
integer, intent(in) :: nvec2, itr
complex(KIND(0d0)), intent(inout) :: x(:)
complex(KIND(0d0)), intent(inout) :: v1(:), v0(:)
double precision, intent(inout) :: ene(:)
logical, intent(in) :: restart
integer :: i, j, k, start
double precision :: prdct, beta1, alpha1, dnorm, dd
complex(KIND(0d0)) :: temp1, temp2, dcnorm, cprdct
integer :: t_max_itr
t_max_itr=itr_counter(4)
!!
!*** initialization
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(ldim) &
!$OMP& PRIVATE(i)
!$OMP do
do i=1,ldim
  v0(i)=dcmplx(0d0,0d0)
  v1(i)=dcmplx(0d0,0d0)
enddo
!$OMP END DO
!$OMP END PARALLEL
!!!!!!
!!
if(itr_counter(3)<3)then
  do k=1,1
    if(.not.restart)then
!$OMP PARALLEL&
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(k,ldim) &
!$OMP& PRIVATE(i)
!$OMP do
      do i=1,ldim
        x(i)=dcmplx(0d0,0d0)
      enddo
!$OMP END DO
!$OMP END PARALLEL
      call parallelized_set_initial_vec_p(x)
    endif
    dd = coef(1,k)
!$OMP PARALLEL&
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(ldim,k,dd) &
!$OMP& PRIVATE(i)
!$OMP do
    do i=1,ldim
      x(i)=x(i)*dd
    enddo
!$OMP END DO
!$OMP END PARALLEL
  enddo
!*** alpha(1) and beta(1)
  call p_mltply_v3(v1,v0,cprdct,0d0)
  alpha1=alpha(1)
  beta1=beta(1)
  do k=1,1  !do k=1,nvec2
    dd = coef(2,k)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(alpha1,beta1,k,ldim,dd) &
!$OMP& PRIVATE(j) 
!$OMP do 
    do j=1,ldim
      x(j)=x(j)+dd*(v0(j)-alpha1*v1(j))/beta1
    enddo
!$OMP END DO
!$OMP END PARALLEL
  enddo
else
  call interval_dump_read(position=progress, data_size=max_ldim, &
                          vec_size=t_max_itr, sa=alpha, sb=beta ,sc=coef, se=ene, &
                          vec0=x, vec1=v0, vec2=v1)
  alpha1=alpha(itr_counter(3)-1)
  beta1=beta(itr_counter(3)-1)
endif
!!
!*** iteration
start=itr_counter(3)
do i=start,itr_counter(4)
!  if(mod(i,5)==0) write(*,'(A,I8)') 'Itr :',i
  write(*,'(A,I8)') 'Itr :',i
  if(tell_me_time()>max_comp_time)then
    call interval_dump(position=progress, data_size=max_ldim, &
                       vec_size=t_max_itr, sa=alpha, sb=beta, sc=coef, se=ene, &
                       vec0=x, vec1=v0, vec2=v1)
  endif  
!!
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(i,ldim,alpha1,beta1) &
!$OMP& PRIVATE(temp1,temp2,j) 
!$OMP do 
  do j=1,ldim
    temp1=v1(j)
    temp2=(v0(j)-alpha1*v1(j))/beta1
    v0(j)=-beta1*temp1
    v1(j)=temp2
  enddo
!$OMP END DO
!$OMP END PARALLEL
  call p_mltply_v3(v1,v0,cprdct,0d0)
  alpha1=alpha(i)
  beta1=beta(i)
!!
  do k=1,1
    dd = coef(i+1,k)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(alpha1,beta1,i,k,dd,ldim) &
!$OMP& PRIVATE(j)
!$OMP do 
    do j=1,ldim
      x(j)=x(j)+dd*(v0(j)-alpha1*v1(j))/beta1
    enddo
!$OMP END DO
!$OMP END PARALLEL
  enddo
  itr_counter(3)=itr_counter(3)+1
enddo
!
!*** normalization
do k=1,1
  call calc_vec_inprdct(ldim,x,x,dcnorm)
  dnorm=dsqrt(DREAL(dcnorm))
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(dnorm,ldim,k) &
!$OMP& PRIVATE(j)
!$OMP do 
  do j=1,ldim
    x(j)=x(j)/dnorm
  enddo
!$OMP END DO
!$OMP END PARALLEL
enddo
progress(2)=.false.
end subroutine parallelized_lncv1z
!!
end module parallelized_lanczos_method
!**********************************************************************************************************************
