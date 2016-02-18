module evaluate
implicit none
double precision, allocatable :: alpha(:,:),beta(:,:),coef(:,:,:)
double precision, allocatable :: fac1(:), fac2(:), dat1(:), dat2(:)
character(LEN=2) :: comp='tp'
integer :: seed_num, num_data1, num_data2
integer :: pme=0

contains
subroutine get_args(bt)
implicit none
character(LEN=256) :: cbt
integer :: num
double precision, intent(out) :: bt
num=iargc()
if(num==1)then
  call getarg(1,cbt)
  read(cbt,*) bt 
else
  bt=1d0
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
logical, intent(inout) :: w_logic
character(LEN=1024) :: COMM
integer :: j
w_logic=.false.
write(COMM,'(a,I5.5,a,a,a)') 'ls | grep SalphaSbeta_me',pme,'seed | grep ',trim(comp),' | wc -l > list_fort'
call system(COMM)
!!
open(10,file='list_fort', status='old', err=100)
read(10,*,end=25,err=100) seed_num
w_logic=.true.
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
use evaluate
implicit none
integer :: seed
double precision :: bt, res(2)
double precision, allocatable :: tmp(:)
logical :: w_logic
w_logic=.false.
call get_args(bt)
call check_dump_data(w_logic)
if(w_logic)then
call read_alpha_beta_coef

allocate(tmp(seed_num))
do seed=1,seed_num
  call calc_fac(bt,seed)
  dat1(:)=coef(:,4,seed); dat2(:)=coef(:,5,seed)
!  coef(1,4,:)=0d0; coef(1,5,:)=0d0
  call eval_q(seed,bt,dat1,dat2,tmp(seed))
enddo
  call calc_ave(tmp,res)
write(*,*) bt,res(1),res(2)
deallocate(tmp)
endif
call dealloc_mem

end




