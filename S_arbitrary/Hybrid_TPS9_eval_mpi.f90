module my_mpi_lib
implicit none
include 'mpif.h'
integer, public, save :: ierr, np, me

contains
subroutine my_mpi_start()
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NP,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,ME,ierr)
end subroutine my_mpi_start

subroutine my_mpi_end()
implicit none
call barrier
call MPI_FINALIZE(ierr)
end subroutine my_mpi_end

subroutine barrier
implicit none
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine barrier

subroutine bcast_logical_data(l_dat)
implicit none
integer, intent(inout) :: l_dat
call MPI_BCAST(l_dat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
end subroutine bcast_logical_data

subroutine bcast_int_data(num_data1,num_data2,seed_num,tmin,tmax,num_div)
implicit none
integer, intent(inout) :: num_data1, num_data2, seed_num, num_div
double precision, intent(inout) :: tmin, tmax
integer :: idt(4) 
double precision :: dt(2)
idt(1)=num_data1; idt(2)=num_data2; idt(3)=seed_num; idt(4)=num_div
call MPI_BCAST(idt,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
num_data1=idt(1);num_data2=idt(2);seed_num=idt(3);num_div=idt(4)

dt(1)=tmin;dt(2)=tmax
call MPI_BCAST(dt,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
tmin=dt(1);tmax=dt(2)
end subroutine bcast_int_data


subroutine bcast_double_data(num_data1,num_data2,seed_num,alpha,beta,coef)
implicit none
integer, intent(in) :: num_data1, num_data2, seed_num
double precision :: alpha(num_data1,seed_num),beta(num_data1,seed_num),coef(num_data1,num_data2,seed_num)
double precision, allocatable :: tmp(:)
integer :: i1,i2,i3,i4
allocate(tmp(num_data1*seed_num))
do i2=1,seed_num
do i1=1,num_data1
  i3= (i2-1)*num_data1 + i1
  tmp(i3)=alpha(i1,i2)
enddo
enddo
call MPI_BCAST(tmp,num_data1*seed_num,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
do i2=1,seed_num
do i1=1,num_data1
  i3= (i2-1)*num_data1 + i1
  alpha(i1,i2)=tmp(i3)
enddo
enddo
deallocate(tmp)
call barrier
allocate(tmp(num_data1*seed_num))
do i2=1,seed_num
do i1=1,num_data1
  i3= (i2-1)*num_data1 + i1
  tmp(i3)=beta(i1,i2)
enddo
enddo
call MPI_BCAST(tmp,num_data1*seed_num,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
do i2=1,seed_num
do i1=1,num_data1
  i3= (i2-1)*num_data1 + i1
  beta(i1,i2)=tmp(i3)
enddo
enddo
deallocate(tmp)
call barrier
allocate(tmp(num_data1*num_data2*seed_num))
do i3=1,seed_num
do i2=1,num_data2
do i1=1,num_data1
  i4= (i3-1)*num_data1*num_data2 + (i2-1)*num_data1 + i1
  tmp(i4)=coef(i1,i2,i3)
enddo
enddo
enddo
call MPI_BCAST(tmp,num_data1*num_data2*seed_num,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
do i3=1,seed_num
do i2=1,num_data2
do i1=1,num_data1
  i4= (i3-1)*num_data1*num_data2 + (i2-1)*num_data1 + i1
  coef(i1,i2,i3)=tmp(i4)
enddo
enddo
enddo
deallocate(tmp)
call barrier
end subroutine bcast_double_data


subroutine get_all_data(data_num,dat1,dat2)
implicit none
double precision, intent(inout) :: dat1(0:),dat2(0:)
integer, intent(in) :: data_num
integer :: ierr1
call MPI_ALLREDUCE(dat1,dat2,data_num,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr1)
end subroutine get_all_data

end module my_mpi_lib

module evaluate
implicit none
double precision, allocatable :: alpha(:,:),beta(:,:),coef(:,:,:)
double precision, allocatable :: fac1(:), fac2(:), dat1(:), dat2(:)
character(LEN=2) :: comp='tp'
integer :: seed_num, num_data1, num_data2
integer :: pme=0
double precision :: t_min, t_max
integer :: num_div


contains
subroutine get_args(t1,t2,t_num)
implicit none
character(LEN=256) :: ct_min, ct_max, cnum_div
integer :: num
double precision, intent(inout) :: t1, t2
integer :: t_num
num=iargc()
if(num==3)then
  call getarg(1,ct_min)
  read(ct_min,*) t1
  call getarg(2,ct_max)
  read(ct_max,*) t2
  call getarg(3,cnum_div)
  read(cnum_div,*) t_num
else
  t1=1d-2
  t2=1d0
  t_num=1000
endif
end subroutine get_args

subroutine read_alpha_beta_coef
implicit none
character(LEN=128) :: FILENAME, IDENTIFY
integer :: j, k, m, indx1
double precision, allocatable :: ta(:),tb(:),tc(:,:)
j=100+pme+1
indx1=1
write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'SalphaSbeta_me',pme,'seed',indx1,trim(comp)
FILENAME=trim(IDENTIFY)
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
read(j) num_data1
read(j) num_data2
close(j)
!!
allocate(alpha(num_data1,seed_num),beta(num_data1,seed_num),coef(num_data1,num_data2,seed_num))
allocate(ta(num_data1),tb(num_data1),tc(num_data1,num_data2))
allocate(fac1(num_data1),fac2(num_data1),dat1(num_data1),dat2(num_data1))
do indx1=1,seed_num
  write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'SalphaSbeta_me',pme,'seed',indx1,trim(comp)
  FILENAME=trim(IDENTIFY)
  open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
  read(j) num_data1
  read(j) num_data2
  read(j) ta
  read(j) tb
  read(j) tc
  close(j)
  do k =1 ,num_data1
    alpha(k,indx1)=ta(k); beta(k,indx1)=tb(k)
  enddo
  do m= 1 ,num_data2
  do k =1 ,num_data1
    coef(k,m,indx1)=tc(k,m)
  enddo
  enddo
enddo
deallocate(ta,tb,tc)
end subroutine read_alpha_beta_coef

!!
subroutine alloc_mem
implicit none
if(.not.allocated(alpha)) allocate(alpha(num_data1,seed_num))
if(.not.allocated(beta)) allocate(beta(num_data1,seed_num))
if(.not.allocated(coef)) allocate(coef(num_data1,num_data2,seed_num))
if(.not.allocated(fac1)) allocate(fac1(num_data1))
if(.not.allocated(fac2)) allocate(fac2(num_data1))
if(.not.allocated(dat1)) allocate(dat1(num_data1))
if(.not.allocated(dat2)) allocate(dat2(num_data1))
alpha(:,:)=0; beta(:,:)=0;coef(:,:,:)=0
end subroutine alloc_mem

subroutine dealloc_mem
implicit none
if(allocated(alpha)) deallocate(alpha)
if(allocated(beta)) deallocate(beta)
if(allocated(coef)) deallocate(coef)
if(allocated(fac1)) deallocate(fac1)
if(allocated(fac2)) deallocate(fac2)
if(allocated(dat1)) deallocate(dat1)
if(allocated(dat2)) deallocate(dat2)
end subroutine dealloc_mem


subroutine check_dump_data(w_logic)
implicit none
integer, intent(inout) :: w_logic
character(LEN=1024) :: COMM
integer :: j
w_logic=0
write(COMM,'(a,I5.5,a,a,a)') 'ls | grep SalphaSbeta_me',pme,'seed | grep ',trim(comp),' | wc -l > list_fort'
call system(COMM)
!!
open(10,file='list_fort', status='old', err=100)
read(10,*,end=25,err=100) seed_num
w_logic=1
return
!!
25 continue
100 continue
write(*,*) 'error : Data files do not exist.'
close(10)
return
end subroutine check_dump_data
!!
!!
subroutine calc_fac(bt,seed)
implicit none
double precision, intent(in) :: bt
integer, intent(in) :: seed
double precision :: y, z, maxfac
integer  :: m, i
!!!!!
!!!!!
y=0d0
maxfac=-1d200
fac1(:)=0d0; fac2(:)=0d0 ! val : log f
do i=1,2*(num_data1-1)
  m= ishft((i+1),-1) + 1
  z=dlog( bt/dble(2*i)*beta(m,seed))
  y = y + z
  fac1(m)=y
  maxfac=dmax1(y,maxfac)
enddo
do i=1,num_data1-1
  fac1(i)=dexp(fac1(i)-maxfac)
  fac2(i)=fac1(i)*bt*beta(i+1,seed) / dble(2*i-1)
enddo
fac1(num_data1)=dexp(fac1(num_data1)-maxfac)
fac2(num_data1)=0d0
write(*,*) 'R maxfac',maxfac
!!!!!
!!!!!
!y=0d0
!maxfac=-1d200
!do i=1,2*num_data1
!  m= ishft((i-1),-1) + 1
!  z=dlog( bt/dble(2*i)*beta(m,seed))
!  y = y + z
!  fac1(m)=y
!  maxfac=dmax1(y,maxfac)
!enddo
!fac1(1)=dexp(fac1(1)-maxfac)
!fac2(1)=0d0
!do i=2,num_data1
!!!  write(*,*) i, fac1(i)
!  fac1(i)=dexp(fac1(i)-maxfac)
!  fac2(i)=fac1(i-1)*bt / dble(2*i-1)
!enddo
!write(*,*) 'R maxfac',maxfac
!!!!!
end subroutine calc_fac
!!
!!
subroutine eval_q(seed,bt,dat1,dat2,res)
implicit none
integer, intent(in) :: seed
double precision, intent(in) :: bt
double precision, intent(in) :: dat1(:), dat2(:)
double precision, intent(inout) :: res
integer :: i ,j 
integer :: sorted_list(num_data1-1)
double precision :: x, y, z
double precision :: temp(num_data1-1,2)

  call calc_fac(bt,seed)
  z=0d0
  do i=1,num_data1-1
    temp(i,1)=fac1(i)*dat1(i) 
    temp(i,2)=fac2(i)*dat2(i)
    z = z + fac1(i) + fac2(i)
  enddo
  call lsort(temp(:,1),sorted_list)
  x=0d0
  do i=1,num_data1-1
    x = x + temp(sorted_list(i),1)
  enddo
  call lsort(temp(:,2),sorted_list)
  y=0d0
  do i=1,num_data1-1
    y = y + temp(sorted_list(i),2)
  enddo
  res= (x+y)/z
end subroutine eval_q
!!
subroutine lsort(data_list,sorted_list)
implicit none
double precision, intent(in)    :: data_list(:)
integer,          intent(inout) :: sorted_list(:)
double precision, allocatable   :: tmp_dat(:)
integer, allocatable :: itmp_dat(:)
integer :: i, j, itmp, ub, lb
double precision :: tmp
ub=ubound(data_list,1)
lb=lbound(data_list,1)
allocate(tmp_dat(lb:ub),itmp_dat(lb:ub))

do i=lb,ub
  itmp_dat(i)=i
  tmp_dat(i)=data_list(i)
enddo
do j=lb,ub-1
do i=j,ub
  if(tmp_dat(i)>tmp_dat(j))then
    itmp=itmp_dat(i)
    itmp_dat(i)=itmp_dat(j)
    itmp_dat(j)=itmp
    tmp=tmp_dat(i)
    tmp_dat(i)=tmp_dat(j)
    tmp_dat(j)=tmp
  endif
enddo
enddo
do i=lb,ub
  sorted_list(i)=itmp_dat(i)
enddo
deallocate(tmp_dat,itmp_dat)
end subroutine lsort

subroutine calc_ave(dat,res)
implicit none
double precision, intent(inout) :: res(2)
double precision, intent(in) :: dat(:)
integer :: i
double precision :: ave, var
ave=0d0
var=0d0
do i =1,seed_num
  ave=ave+dat(i)
  var=var+dat(i)*dat(i)
enddo
res(1)=ave/dble(seed_num)
res(2)=var/dble(seed_num) - res(1)*res(1)
end subroutine calc_ave

end module evaluate
!!
program main
use my_mpi_lib
use evaluate
implicit none
integer :: seed
double precision :: bt, res(2)
double precision, allocatable :: tmp(:)
integer :: w_logic
integer :: i1, i2
double precision :: t0, log_t, log_t_min, log_t_max, delta, t_delta
double precision :: dat(0:8)
double precision :: d_e, spec, var_spec, spec_d1, spec_d2, del_c
double precision, allocatable :: tres(:), fres(:)

w_logic=0
call get_args(t_min,t_max,num_div)
call my_mpi_start()

if(me==0) call check_dump_data(w_logic)
call barrier
call bcast_logical_data(w_logic)

if(w_logic==1)then
  call barrier
  if(me==0) call read_alpha_beta_coef
  call bcast_int_data(num_data1,num_data2,seed_num,t_min,t_max,num_div)
  if(me/=0) call alloc_mem
  call barrier
  call bcast_double_data(num_data1,num_data2,seed_num,alpha,beta,coef)

  allocate(tmp(seed_num),tres(0:4*(num_div+1)-1),fres(0:4*(num_div+1)-1))
  tres(:)=0
  log_t_max = dlog(t_max)
  log_t_min = dlog(t_min)
  delta=(log_t_max - log_t_min)/dble(num_div)

  t_delta=1d-4
  if ( t_min<t_delta ) t_delta = t_min * 0.9

  do i1 = me, num_div, np
    log_t = log_t_min + dble(i1) * delta
    t0 = dexp(log_t)
    do i2 = 0,2
      if(i2==0) bt = 1.d0/(t0 - t_delta)
      if(i2==1) bt = 1.0d0/t0
      if(i2==2) bt = 1.d0/(t0 + t_delta)
      do seed=1,seed_num
        call calc_fac(bt,seed)
        dat1(:)=coef(:,4,seed); dat2(:)=coef(:,5,seed)
!!  coef(1,4,:)=0d0; coef(1,5,:)=0d0
        call eval_q(seed,bt,dat1,dat2,tmp(seed))
      enddo
      call calc_ave(tmp,res)
      dat(i2*3+0)=bt; dat(i2*3+1)=res(1); dat(i2*3+2)=res(2)
!!  write(*,*) bt,res(1),res(2)
    enddo

    d_e = 0.5d0*(dat(7) - dat(1))
    spec = d_e/t_delta
    var_spec = 0.5d0*(dat(2)+dat(8))
    spec_d1 = (dat(4)-dat(1))/t_delta
    spec_d2 = (dat(7)-dat(4))/t_delta
    del_c   = (spec_d2 - spec_d1)/t_delta
  
    tres(4*i1)=1d0/dat(3); tres(4*i1+1)=spec; tres(4*i1+2)=var_spec; tres(4*i1+3)=del_c
!    write(*,'(4E20.6)') 1d0/dat(3),spec,var_spec,del_c


  enddo

  deallocate(tmp)
endif
call barrier
call get_all_data(4*(num_div+1),tres,fres)

if(me==0)then
  do i1=0, num_div
    write(*,'(4E20.6)') fres(4*i1),fres(4*i1+1),fres(4*i1+2),fres(4*i1+3)
  enddo
endif

call barrier
call dealloc_mem
deallocate(tres,fres)
call my_mpi_end
end




