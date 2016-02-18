module DynamicStractureFactor
implicit none
private
public calc_DSF, calc_extended_DSF, intalv_zz, intalv_xy

contains
subroutine clear_dsf_mem
use working_area, only : salpha, sbeta
implicit none
!$OMP PARALLEL
!$OMP WORKSHARE
salpha(:)=0d0
sbeta(:)=0d0
!$OMP END WORKSHARE
!$OMP END PARALLEL
end subroutine clear_dsf_mem

subroutine intalv_zz(v,x,q,cnorm)
use parameters, only : n, ldim, max_ldim, b_bit
use data_mem, only : site_vec, list1
use parallel_lib, only : lme, barrier
use parallelized_diag, only : calc_vec_inprdct
use openmp_lib
implicit none
double precision, intent(in) :: q(:)
complex(KIND(0d0)), intent(in) :: x(:)
complex(KIND(0d0)), intent(inout) :: v(:), cnorm
double precision :: qnorm, norm2
integer :: i, j, is, istate, ibit, bit_shift, t_lme
logical :: bp
!! Sz(Q)
!$OMP PARALLEL 
!$OMP do
do i =1,max_ldim
  v(i)=dcmplx(0d0,0d0)
enddo
!$OMP END DO
!$OMP END PARALLEL 
!!
t_lme=lme
do j=0, n-1
  qnorm = q(1) * site_vec(1,j+1) + q(2) * site_vec(2,j+1)
  bit_shift=j
  bp = .false.
  if(j>b_bit-1)then
    bit_shift=j-(b_bit)
    bp = .true.
  endif
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,istate,is) &
!$OMP& FIRSTPRIVATE(qnorm,t_lme,bit_shift,ibit,bp)
!$OMP do 
  do i=1, ldim
           istate=list1(i)
           is=iand(ishft(istate,-bit_shift),1)
    if(bp) is=iand(ishft(t_lme, -bit_shift),1)
    v(i)=v(i)+(dble(is)-0.5d0)*dcmplx(dcos(qnorm)*DREAL(x(i))-dsin(qnorm)*DIMAG(x(i)),&
                                      dsin(qnorm)*DREAL(x(i))+dcos(qnorm)*DIMAG(x(i)))
  enddo
!$OMP END DO
!$OMP END PARALLEL
enddo
!!
!******************* normarize *****************************
!!
cnorm=dcmplx(0d0,0d0)
call calc_vec_inprdct(ldim,v,v,cnorm)
if(CDABS(cnorm)<1d-25)then
  write(*,*) 'Warning : vec norm is too small.',cnorm
  cnorm=dcmplx(1d-25,0d0)
!!  stop
endif
norm2=1d0/dsqrt(DREAL(cnorm))
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(ldim,norm2)
!$OMP do 
do i=1, ldim
  v(i)=v(i)*norm2
enddo
!$OMP END DO
!$OMP END PARALLEL
cnorm=cnorm/dble(n)
write(*,*) 'sq',cnorm

end subroutine intalv_zz

subroutine intalv_xy(v,x,q,cnorm)
use parameters, only : n, ldim, max_ldim, max_tdim, pdim, b_bit, total_sz
use data_mem, only : site_vec, list3, listt
use parallel_lib, only : lme, barrier, wait, exchange_data
use parallelized_diag, only : calc_vec_inprdct, get_mask, alloc_listt, get_listt, dealloc_listt
use openmp_lib
implicit none
double precision, intent(in) :: q(:)
complex(KIND(0d0)), intent(in) :: x(:)
complex(KIND(0d0)), intent(inout) :: v(:), cnorm

double precision :: qnorm, norm2
integer :: i, j
integer :: is, it, it0, it1, it2, lipm
complex(KIND(0d0)), allocatable :: r(:)
integer :: target_cpu, ibit
integer :: irght2, ilft2, ihf2
integer, parameter :: ipm=1
call get_mask(irght2,ilft2,ihf2)
!$OMP PARALLEL
!$OMP DO
do i=1,max_ldim
  v(i)=dcmplx(0d0,0d0)
enddo
!$OMP END DO
!$OMP END PARALLEL
!ipm = -1 -> S-
!ipm = +1 -> S+
call get_listt(lme,total_sz,lme,pdim)
do j=0, b_bit-1
  qnorm = q(1) * site_vec(1,j+1) + q(2) * site_vec(2,j+1)
  lipm=ipm
  ibit=ishft(1,j)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,is,it0,it1,it2,it)&
!$OMP& FIRSTPRIVATE(j,qnorm,lipm,ibit,irght2,ilft2,ihf2)
!$OMP do 
  do i=1, pdim
    is=ishft(iand(listt(i),ibit),-j)
    if(is==1-lipm)then
      it0=listt(i) + (1-2*is) * ishft(1,j)
      it1=iand(it0,irght2)
      it2=ishft(iand(it0,ilft2),-ihf2)
      it=list3(1,it1,0)+list3(2,it2,0)
      v(it)=v(it)+dcmplx(dcos(qnorm)*DREAL(x(i))-dsin(qnorm)*DIMAG(x(i)),&
                         dsin(qnorm)*DREAL(x(i))+dcos(qnorm)*DIMAG(x(i)))
    endif
  enddo
!$OMP END DO
!$OMP END PARALLEL
enddo
call dealloc_listt
!!
!!
allocate(r(1:max_tdim))
do j=b_bit, n-1
  qnorm = q(1) * site_vec(1,j+1) + q(2) * site_vec(2,j+1)
  is=iand(ishft(lme,-(j-b_bit)),1)
  target_cpu = lme + (1-2*is) * ishft(1,(j-b_bit))
  call get_listt(target_cpu,total_sz,lme,pdim)
  call exchange_data(target_cpu,max_tdim,j,x,r)
  if(is==ipm)then
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,it0,it1,it2,it)&
!$OMP& FIRSTPRIVATE(j,qnorm,pdim,ilft2,irght2,ihf2)
!$OMP do 
    do i=1, pdim
      it0=listt(i)
      it1=iand(it0,irght2)
      it2=ishft(iand(it0,ilft2),-ihf2)
      it=list3(1,it1,0)+list3(2,it2,0)
      v(it)=v(it)+dcmplx(dcos(qnorm)*DREAL(r(i))-dsin(qnorm)*DIMAG(r(i)),&
                         dsin(qnorm)*DREAL(r(i))+dcos(qnorm)*DIMAG(r(i)))
    enddo
!$OMP END DO
!$OMP END PARALLEL
  endif
  call dealloc_listt
  call barrier
enddo
deallocate(r)
!!
!******************* normarize *****************************
!!
call calc_vec_inprdct(ldim,v,v,cnorm)
if(CDABS(cnorm)<1d-25)then
  write(*,*) 'Warning : vec norm is too small.'
  cnorm=dcmplx(1d-25,0d0)
!  stop
endif
norm2=1d0/dsqrt(DREAL(cnorm))
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(norm2)
!$OMP do 
do i=1, ldim
  v(i)=v(i)*norm2
enddo
!$OMP END DO
!$OMP END PARALLEL
cnorm=cnorm/dble(n)
write(*,*) 'sq',cnorm
end subroutine intalv_xy
!!
subroutine parallelized_dsflnc1z(iq,icf,cnorm,v1,v0)
use parameters, only : max_ldim, ldim, s_maxitr, itr_counter, max_comp_time, progress
use working_area, only : salpha, sbeta
use parallelized_diag, only : calc_alpha_beta
use parallel_lib
use openmp_lib
use def_output, only : interval_dump
use timer
implicit none
integer, intent(in) :: iq
complex(KIND(0d0)), intent(inout) :: v1(:),v0(:)
complex(KIND(0d0)), intent(in) :: cnorm
integer, intent(inout) :: icf
integer :: i, j, start, t_s_maxitr
double precision :: alpha1, beta1
complex(KIND(0d0)) :: temp0, temp1
!!
t_s_maxitr=itr_counter(12)
if(itr_counter(11)==2)then
!$OMP PARALLEL 
!$OMP WORKSHARE
  v0(:)=dcmplx(0d0,0d0)
!$OMP END WORKSHARE
!$OMP END PARALLEL
!!
  sbeta(1) =cnorm
  call calc_alpha_beta(ldim,alpha1,beta1,v1,v0)
  salpha(1)=alpha1
  sbeta(2) =beta1
else
  alpha1=salpha(itr_counter(11)-1)
  beta1=sbeta(itr_counter(11))
endif
write(*,*) 'ITR:', itr_counter(11)
write(*,*) 'DSF ITR :  start  :',salpha(itr_counter(11)),sbeta(itr_counter(11))
!!
!*** iteration
icf=1
start=itr_counter(11)
do i=start,itr_counter(12)-1
  if(tell_me_time()>max_comp_time)then
    call interval_dump(position=progress, data_size=max_ldim, &
                       vec_size=t_s_maxitr, sa=salpha, sb=sbeta, &
                       vec1=v1, vec2=v0)
  endif
  if(mod(i,5)==0.and.lme==0) write(*,*) 'DSF ITR :',i
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j,temp1,temp0) &
!$OMP& FIRSTPRIVATE(alpha1,beta1)
!$OMP do 
  do j=1,ldim
    temp1=v1(j)
    temp0=(v0(j)-alpha1*v1(j))/beta1
    v0(j)=-beta1*temp1
    v1(j)=temp0
  enddo
!$OMP END DO
!$OMP END PARALLEL
  call calc_alpha_beta(ldim,alpha1,beta1,v1,v0)
  salpha(i)=alpha1
  sbeta(i+1) =beta1
!  if(mod(i,10)==0) write(*,*) 'DSF ITR :',i,salpha(i),sbeta(i)
  itr_counter(11)=itr_counter(11)+1
enddo
icf=i-1
end subroutine parallelized_dsflnc1z

subroutine exec_renbun(i_range,s_maxitr,split,ene,epsln,w,conv,dsf)
use working_area, only : salpha, sbeta
implicit none
integer, intent(in) :: i_range, s_maxitr
double precision, intent(out) :: w(:), dsf(:)
integer, intent(out) :: conv(:)
double precision, intent(in) :: ene, epsln, split
double precision :: xi, omega
integer :: j, i_conv, m1
double precision :: yy,xr,c,cr,ci,d,dr,di,dl,dlr,dli,f1,f2,xr1,xi1
double precision :: K_split, K_epsln
yy=ene
K_epsln=epsln
K_split=split
!$OMP PARALLEL &
!$OMP DEFAULT(PRIVATE)&
!$OMP SHARED(i_range,s_maxitr,salpha,sbeta,w,conv,dsf) &
!$OMP FIRSTPRIVATE(K_split,yy,K_epsln) 
!$OMP do SCHEDULE(STATIC,1)
do j=1,i_range
  omega=dble(j-1)*K_split
  w(j)=omega
  i_conv=0
  xr=1.0d-10
  xi=1.0d-10
  cr=xr
  ci=xi
  dr=0.0d0
  di=0.0d0
  do m1=1,s_maxitr-1
    f1=dble(omega+yy-salpha(m1))
    f2=dble((-1.0d0)*(sbeta(m1)**2.0d0))
    xr1=xr
    xi1=xi
    dr=dr*f2
    dr=dr+f1
    di=di*f2
    di=di+K_epsln
    d=dr**2+di**2
    if(d==0d0)then
      dr=0.5d-28
      di=0.5d-28
    endif
    c=cr**2+ci**2
    cr=cr*f2/c
    cr=cr+f1
    ci=-ci*f2/c
    ci=ci+K_epsln
    c=cr**2+ci**2
    if(c==0d0)then
      cr=0.5d-28
      ci=0.5d-28
    endif
    d=dr**2+di**2
    dr=dr/d
    di=-di/d
    dlr=cr*dr-ci*di
    dli=ci*dr+cr*di
    xr=xr1*dlr-xi1*dli
    xi=xr1*dli+xi1*dlr
    dl=(dlr-1d0)**2+dli**2
    if (dl<0.5d-28)then
      i_conv=m1
      exit
    endif
  enddo 
  conv(j)=i_conv
  dsf(j)=xi
enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine exec_renbun


subroutine nomalization_DSF(dsf,sq)
use parameters, only : i_range, split
implicit none
double precision, intent(inout) :: dsf(:)
complex(KIND(0d0)) :: sq
double precision :: sum
integer :: i
sum=0.0d0
do i=2,i_range
  sum=sum+0.5d0*split*(dsf(i)+dsf(i-1))
enddo
do i=1,i_range
  dsf(i)=DREAL(sq)/sum*dsf(i)
enddo
end subroutine nomalization_DSF

subroutine dump_alpha_beta(iq,indx,sq,comp,v1,v2)
use parallel_lib, only : get_me
use working_area, only : salpha, sbeta
implicit none
integer, intent(in) :: iq, indx
complex(KIND(0d0)), intent(in) :: sq
complex(KIND(0d0)), intent(in) :: v1(:), v2(:)
character(LEN=2), intent(in) :: comp
character(LEN=128) :: FILENAME, IDENTIFY
integer :: j, pme
pme=get_me()
write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'SalphaSbeta_q',iq,'me',pme,trim(comp)
FILENAME=trim(IDENTIFY)
j=100+pme+1
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
write(j) sq
write(j) salpha
write(j) sbeta
close(j)
write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'vec_q',iq,'me',pme,trim(comp)
FILENAME=trim(IDENTIFY)
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
write(j) indx
write(j) v1
write(j) v2
close(j)
end subroutine dump_alpha_beta


subroutine read_alpha_beta2(iq,indx,sq,comp,v1,v2)
use parallel_lib, only : get_me
use working_area, only : salpha, sbeta
implicit none
integer, intent(in) :: iq
integer, intent(inout) ::  indx
complex(KIND(0d0)), intent(inout) :: sq
complex(KIND(0d0)), intent(inout) :: v1(:), v2(:)
character(LEN=2), intent(in) :: comp
character(LEN=128) :: FILENAME, IDENTIFY
integer :: j, pme

double precision, allocatable :: ta(:), tb(:)

pme=get_me()
j=100+pme+1
write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'vec_q',iq,'me',pme,trim(comp)
FILENAME=trim(IDENTIFY)
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
read(j) indx
read(j) v1
read(j) v2
close(j)

allocate(ta(indx+1),tb(indx+1))
write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'SalphaSbeta_q',iq,'me',pme,trim(comp)
FILENAME=trim(IDENTIFY)
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
read(j) sq
read(j) ta
read(j) tb
close(j)
salpha(:)=0d0; sbeta(:)=0d0
do j =1 ,indx+1
  salpha(j)=ta(j); sbeta(j)=tb(j)
enddo
deallocate(ta,tb)
end subroutine read_alpha_beta2


subroutine read_dump_alpha_beta(iq,sq,comp)
use parameters
use working_area, only : salpha, sbeta
use parallel_lib, only : get_me
implicit none
integer, intent(in) :: iq
complex(KIND(0d0)), intent(inout) :: sq
character(LEN=2),intent(in) :: comp
character(LEN=128) :: FILENAME, IDENTIFY
integer :: j, pme
pme=get_me()
write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'SalphaSbeta_q',iq,'me',pme,trim(comp)
FILENAME=trim(IDENTIFY)
j=100+pme+1
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
read(j) sq
read(j) salpha
read(j) sbeta
close(j)
write(*,*) 'read data :',sq,salpha(1),sbeta(1)
write(*,*) 'read data :',   salpha(20),sbeta(20)
end subroutine read_dump_alpha_beta


subroutine calc_DSF(Iam,ene,x,i_color)
use parameters
use working_area, only : wwk, salpha, sbeta
use data_mem, only : q_list
use def_output
implicit none
integer, intent(in) :: Iam, i_color
double precision, intent(in) :: ene
complex(KIND(0d0)), intent(inout) :: x(:)
integer :: i, j, start, idsfitr
complex(KIND(0d0)) :: norm
double precision :: q(2)
double precision :: dsf(1:i_range), w(1:i_range)
integer :: conv(1:i_range)
integer :: t_s_maxitr
t_s_maxitr=itr_counter(12)
!!
if(Iam==0) write(*,*) 'start : calc_DSF'
call clear_dsf_mem
call open_res_files(Iam,dsf_comp,i_color)
start=itr_counter(9)
do i=start,itr_counter(10)
    q(1:2)=q_list(1:2,i)
    write(*,'(A,I3,2F12.6)') ' q = ',i,q(1),q(2)


  if(itr_counter(11)==2)then
    if(dsf_comp=='zz')then
      call intalv_zz(wwk(:,1),x,q,norm)
    elseif(dsf_comp=='xx')then
      call intalv_xy(wwk(:,1),x,q,norm)
    else
      write(*,*) 'Error: DSF comp.'
      stop
    endif
  else
    call interval_dump_read(position=progress, data_size=max_ldim, &
                       vec_size=t_s_maxitr, sa=salpha, sb=sbeta,  &
                       vec1=wwk(:,1), vec2=wwk(:,2))
    norm=sbeta(1)
  endif


  if(Iam==0)then
  open(23,file='sq.txt',status='unknown',access='append')
    write(23,'(A,I3,2F12.6,F25.18)') ' q = ',i,q(1),q(2),DREAL(norm)
  close(23)
  endif

  if(progress(4))then
    call parallelized_dsflnc1z(i,idsfitr,norm,wwk(:,1),wwk(:,2))
    call dump_alpha_beta(i,idsfitr,norm,dsf_comp,wwk(:,1),wwk(:,2))
  else
    call read_dump_alpha_beta(i,norm,dsf_comp)
  endif
!!
  call exec_renbun(i_range,s_maxitr,split,ene,delt,w,conv,dsf)
!!
  call nomalization_DSF(dsf,norm)
  call write_dsf_res(Iam,i_color,dsf_comp,i_range,norm,q,w,dsf,conv)
  call clear_dsf_mem
  itr_counter(11)=2
  itr_counter(9)=itr_counter(9)+1
enddo
call close_res_files(Iam,i_color)
progress(4)=.false.
end subroutine calc_DSF


subroutine calc_extended_DSF(Iam,ene,x,i_color)
use parameters
use working_area, only : wwk, salpha, sbeta
use data_mem, only : q_list
use def_output
implicit none
integer, intent(in) :: Iam, i_color
double precision, intent(in) :: ene
complex(KIND(0d0)), intent(inout) :: x(:)
integer :: i, j, start, idsfitr, icn
complex(KIND(0d0)) :: norm
double precision :: q(2)
double precision :: dsf(1:i_range), w(1:i_range)
integer :: conv(1:i_range)
!!
call clear_dsf_mem
call open_res_files(Iam,dsf_comp,i_color)
start=itr_counter(9)
do i=start,itr_counter(10)
  q(1:2)=q_list(1:2,i)
  write(*,'(A,I3,2F12.6)') ' q = ',i,q(1),q(2)

  call read_alpha_beta2(i,icn,norm,dsf_comp,wwk(:,1),wwk(:,2))
  itr_counter(11)=icn+1 

  call parallelized_dsflnc1z(i,idsfitr,norm,wwk(:,1),wwk(:,2))
  call dump_alpha_beta(i,idsfitr,norm,dsf_comp,wwk(:,1),wwk(:,2))
!!
  call exec_renbun(i_range,s_maxitr,split,ene,delt,w,conv,dsf)
!!
  call nomalization_DSF(dsf,norm)
  call write_dsf_res(Iam,i_color,dsf_comp,i_range,norm,q,w,dsf,conv)
  call clear_dsf_mem
  itr_counter(11)=2
  itr_counter(9)=itr_counter(9)+1
enddo
call close_res_files(Iam,i_color)
progress(4)=.false.
end subroutine calc_extended_DSF

end module DynamicStractureFactor
