module parallel_lib
!use mpi
implicit none
include 'mpif.h'
private

logical, private :: debug=.true.
integer, private :: LOCAL_MPI_COMM_WORLD
integer, private :: NP, ME, o_np
integer, private :: COLOR, KEY
integer, public, save :: ierr, lnp, lme
integer, allocatable, private, save :: list(:)
integer, allocatable, private, save :: stat_send(:,:),stat_recv(:,:)
integer, allocatable, private, save :: request_send(:), request_recv(:)
integer, allocatable, private, save :: stat_send2(:,:),stat_recv2(:,:)
integer, allocatable, private, save :: request_send2(:), request_recv2(:)

integer, public :: num_lbits_p
integer, public :: num_hbits_p
integer, public :: num_bbits_p
integer, allocatable, public :: cpu_pair_list(:,:,:)
integer, allocatable, public :: lbits_pair(:,:)
integer, allocatable, public :: bbits_pair(:,:)
integer, allocatable, public :: hbits_pair(:,:)
integer, allocatable, public :: ibond_type_p_l(:)
integer, allocatable, public :: ibond_type_p_b(:)
integer, allocatable, public :: ibond_type_p_h(:)
integer, allocatable, public :: exchange_list_b(:),exchange_list_h(:,:)
integer, allocatable, public :: matrix_element(:,:)
logical, allocatable, public :: exlog(:,:)
double precision, allocatable, public :: bondwt_p_l(:)
double precision, allocatable, public :: bondwt_p_b(:)
double precision, allocatable, public :: bondwt_p_h(:)
!! for mltply v3
integer, public :: system_mem_size
integer, public :: size_of_eff_mtx
integer, public :: num_int_bond
integer, public :: num_tot_ex_cpus
integer, allocatable, public :: cpu_list(:)
integer, allocatable, public :: int_type(:)
integer, allocatable, public :: site_pair(:,:)
integer, allocatable, public :: matrix_element2(:)
integer, allocatable, public :: tot_ex_cpus(:)
logical, allocatable, public :: exlog2(:)
double precision, allocatable, public :: bondwt_p2(:)
integer, allocatable, public :: ibond_type_p2(:)
integer, public :: max_local_loop, ex_loop_start
double precision, allocatable :: buff(:)
!!

public set_parallel_bits, my_mpi_start, exchange_data, &
       allreduce_rdata, allreduce_data, error_stop, mpi_end, barrier, &
       dealloc_mem_parallel, mpi_get_maxint, wait, make_cpu_list2,mpi_check_logical, &
       mpi_mywait, mpi_test_isendrecv, mpi_isendrecv, mpi_probe_data, &
       get_buffer,release_buffer, get_np, get_me, get_color, my_mpi_split, bcast_sim_data, &
       check_system_memsize


contains
subroutine check_system_memsize
implicit none
character(LEN=1024) :: command
system_mem_size=1024*1024*1024
if(ME==0)then
!  command='free -k -o | head -n 2 | tail -n 1 | awk '//"'"//'{print $2}'//"'"//' > temp'
!  call system(command)
!  open(10,file='temp',status='old')
!  read(10,*) system_mem_size
!  close(10)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine check_system_memsize

subroutine set_parallel_bits
use parameters, only : n, ldim, p_bit, b_bit
implicit none
p_bit=nint(log10(dble(lnp))/log10(2d0))
if(2**p_bit.ne.lnp)then
  write(*,*) 'error : wrong cpu number',p_bit,lnp
  stop
endif
b_bit=n-p_bit
ldim=2**b_bit
end subroutine set_parallel_bits

function get_np()
implicit none
integer :: get_np
  get_np=NP
end function get_np

function get_me()
implicit none
integer :: get_me
  get_me=ME
end function get_me

function get_color()
implicit none
integer :: get_color
  get_color=COLOR
end function get_color

subroutine make_cpu_list2
use parameters, only : n, ibond, b_bit, p_bit
use data_mem, only : ipair, bondwt, ibond_type
use openmp_lib
implicit none
integer :: i, j, k, isite1, isite2, start, num_b, num_h, l1, l2, istat,  ibase, num_l
integer, allocatable :: temp_cpu_list(:),temp_int_type(:),temp_matrix_element2(:),temp_site_pair(:,:)
integer, allocatable :: sorted_list(:),temp_ibond_type_p2(:)
logical, allocatable :: temp_exlog2(:)
double precision, allocatable :: temp_bondwt_p2(:)
integer :: pre_cpu

allocate(lbits_pair(2,ibond), bbits_pair(2,ibond), hbits_pair(2,ibond))
allocate(ibond_type_p_l(ibond), ibond_type_p_b(ibond), ibond_type_p_h(ibond))
allocate(bondwt_p_l(ibond), bondwt_p_b(ibond), bondwt_p_h(ibond))
allocate(list(0:p_bit-1))
start=lme
do i=0,p_bit-1
  list(i)=iand(start,1)
  start=start/2
enddo

num_lbits_p=0
num_hbits_p=0
num_bbits_p=0
lbits_pair=-1
bbits_pair=-1
hbits_pair=-1
ibond_type_p_l=-1
ibond_type_p_h=-1
ibond_type_p_b=-1
bondwt_p_l=0d0
bondwt_p_h=0d0
bondwt_p_b=0d0

do i=1,ibond
  isite1=ipair(2*i-1)-1
  isite2=ipair(2*i  )-1
  if((isite1<b_bit).and.(isite2<b_bit))then
    num_lbits_p=num_lbits_p+1
    lbits_pair(1,num_lbits_p)=min(isite1,isite2)
    lbits_pair(2,num_lbits_p)=max(isite1,isite2)
    bondwt_p_l(num_lbits_p)=bondwt(i)
    ibond_type_p_l(num_lbits_p)=ibond_type(i)

  elseif(((isite1<b_bit).and.(isite2>=b_bit)).or.((isite1>=b_bit).and.(isite2<b_bit)))then
    num_bbits_p=num_bbits_p+1  
    bbits_pair(1,num_bbits_p)=min(isite1,isite2)
    bbits_pair(2,num_bbits_p)=max(isite1,isite2) - b_bit
    bondwt_p_b(num_bbits_p)=bondwt(i)
    ibond_type_p_b(num_bbits_p)=ibond_type(i)

  elseif((isite1>=b_bit).and.(isite2>=b_bit))then
    num_hbits_p=num_hbits_p+1  
    hbits_pair(1,num_hbits_p)=min(isite1,isite2) - b_bit
    hbits_pair(2,num_hbits_p)=max(isite1,isite2) - b_bit
    bondwt_p_h(num_hbits_p)=bondwt(i)
    ibond_type_p_h(num_hbits_p)=ibond_type(i)

  else
    write(*,*) 'error: #=', i,isite1,isite2
  endif
enddo
!!
!!
num_l=max(num_lbits_p,1)
num_b=max(num_bbits_p,1)
num_h=max(num_hbits_p,1)
size_of_eff_mtx=4
num_int_bond=size_of_eff_mtx*ibond
allocate(temp_cpu_list(1:num_int_bond),temp_int_type(1:num_int_bond),sorted_list(1:num_int_bond), &
         temp_matrix_element2(1:num_int_bond),temp_exlog2(1:num_int_bond),temp_site_pair(2,1:num_int_bond), &
         temp_bondwt_p2(1:num_int_bond),temp_ibond_type_p2(1:num_int_bond))
temp_exlog2(:)=.false.


do i=1,num_lbits_p
  temp_cpu_list(4*i-3)=lme
  temp_cpu_list(4*i-2)=lme
  temp_cpu_list(4*i-1)=lme
  temp_cpu_list(4*i-0)=lme
  temp_int_type(4*i-3)=1
  temp_int_type(4*i-2)=1
  temp_int_type(4*i-1)=1
  temp_int_type(4*i-0)=1
  temp_matrix_element2(4*i-3)=0
  temp_matrix_element2(4*i-2)=1
  temp_matrix_element2(4*i-1)=2
  temp_matrix_element2(4*i-0)=3
  temp_exlog2(4*i-3)=.true.
  temp_exlog2(4*i-2)=.true.
  temp_exlog2(4*i-1)=.true.
  temp_exlog2(4*i-0)=.true.
  temp_site_pair(1,4*i-3)=lbits_pair(1,i);temp_site_pair(2,4*i-3)=lbits_pair(2,i)
  temp_site_pair(1,4*i-2)=lbits_pair(1,i);temp_site_pair(2,4*i-2)=lbits_pair(2,i)
  temp_site_pair(1,4*i-1)=lbits_pair(1,i);temp_site_pair(2,4*i-1)=lbits_pair(2,i)
  temp_site_pair(1,4*i-0)=lbits_pair(1,i);temp_site_pair(2,4*i-0)=lbits_pair(2,i)
  temp_ibond_type_p2(4*i-3)=ibond_type_p_l(i)
  temp_ibond_type_p2(4*i-2)=ibond_type_p_l(i)
  temp_ibond_type_p2(4*i-1)=ibond_type_p_l(i)
  temp_ibond_type_p2(4*i-0)=ibond_type_p_l(i)
  temp_bondwt_p2(4*i-3)=bondwt_p_l(i)
  temp_bondwt_p2(4*i-2)=bondwt_p_l(i)
  temp_bondwt_p2(4*i-1)=bondwt_p_l(i)
  temp_bondwt_p2(4*i-0)=bondwt_p_l(i)
enddo
do i=1,num_bbits_p
  j=i+num_lbits_p
  isite1=bbits_pair(1,i)
  isite2=bbits_pair(2,i)
  temp_cpu_list(4*j-3)=lme
  temp_cpu_list(4*j-2)=lme
  temp_cpu_list(4*j-1)=lme+(1-2*list(isite2))*2**isite2
  temp_cpu_list(4*j-0)=lme+(1-2*list(isite2))*2**isite2
  temp_int_type(4*j-3)=2
  temp_int_type(4*j-2)=2
  temp_int_type(4*j-1)=2
  temp_int_type(4*j-0)=2
  temp_matrix_element2(4*j-3)=0
  temp_matrix_element2(4*j-2)=1
  temp_matrix_element2(4*j-1)=2
  temp_matrix_element2(4*j-0)=3
  temp_exlog2(4*j-3)=.true.
  temp_exlog2(4*j-2)=.true.
  temp_exlog2(4*j-1)=.true.
  temp_exlog2(4*j-0)=.true.
  temp_site_pair(1,4*j-3)=bbits_pair(1,i);temp_site_pair(2,4*j-3)=bbits_pair(2,i)
  temp_site_pair(1,4*j-2)=bbits_pair(1,i);temp_site_pair(2,4*j-2)=bbits_pair(2,i)
  temp_site_pair(1,4*j-1)=bbits_pair(1,i);temp_site_pair(2,4*j-1)=bbits_pair(2,i)
  temp_site_pair(1,4*j-0)=bbits_pair(1,i);temp_site_pair(2,4*j-0)=bbits_pair(2,i)
  temp_ibond_type_p2(4*j-3)=ibond_type_p_b(i)
  temp_ibond_type_p2(4*j-2)=ibond_type_p_b(i)
  temp_ibond_type_p2(4*j-1)=ibond_type_p_b(i)
  temp_ibond_type_p2(4*j-0)=ibond_type_p_b(i)
  temp_bondwt_p2(4*j-3)=bondwt_p_b(i)
  temp_bondwt_p2(4*j-2)=bondwt_p_b(i)
  temp_bondwt_p2(4*j-1)=bondwt_p_b(i)
  temp_bondwt_p2(4*j-0)=bondwt_p_b(i)
enddo
do i=1,num_hbits_p
  k=i+num_lbits_p+num_bbits_p
  isite1=hbits_pair(1,i)
  isite2=hbits_pair(2,i)
  istat=2*list(isite2)+list(isite1)
  ibase=lme -list(isite1)*2**isite1 -list(isite2)*2**isite2
  temp_site_pair(1,4*k-3)=hbits_pair(1,i);temp_site_pair(2,4*k-3)=hbits_pair(2,i)
  temp_site_pair(1,4*k-2)=hbits_pair(1,i);temp_site_pair(2,4*k-2)=hbits_pair(2,i)
  temp_site_pair(1,4*k-1)=hbits_pair(1,i);temp_site_pair(2,4*k-1)=hbits_pair(2,i)
  temp_site_pair(1,4*k-0)=hbits_pair(1,i);temp_site_pair(2,4*k-0)=hbits_pair(2,i)
  temp_ibond_type_p2(4*k-3)=ibond_type_p_h(i)
  temp_ibond_type_p2(4*k-2)=ibond_type_p_h(i)
  temp_ibond_type_p2(4*k-1)=ibond_type_p_h(i)
  temp_ibond_type_p2(4*k-0)=ibond_type_p_h(i)
  temp_bondwt_p2(4*k-3)=bondwt_p_h(i)
  temp_bondwt_p2(4*k-2)=bondwt_p_h(i)
  temp_bondwt_p2(4*k-1)=bondwt_p_h(i)
  temp_bondwt_p2(4*k-0)=bondwt_p_h(i)
   if(istat==0)then
  temp_matrix_element2(4*k-3)=0
  temp_matrix_element2(4*k-2)=1
  temp_matrix_element2(4*k-1)=2
  temp_matrix_element2(4*k-0)=3
   elseif(istat==1)then
  temp_matrix_element2(4*k-3)=3
  temp_matrix_element2(4*k-2)=0
  temp_matrix_element2(4*k-1)=1
  temp_matrix_element2(4*k-0)=2
   elseif(istat==2)then
  temp_matrix_element2(4*k-3)=2
  temp_matrix_element2(4*k-2)=1
  temp_matrix_element2(4*k-1)=0
  temp_matrix_element2(4*k-0)=3
   elseif(istat==3)then
  temp_matrix_element2(4*k-3)=3
  temp_matrix_element2(4*k-2)=2
  temp_matrix_element2(4*k-1)=1
  temp_matrix_element2(4*k-0)=0
   endif
   do j=0,3
     l1=iand(temp_matrix_element2(4*k-(3-j)),1)
     l2=iand(ishft(temp_matrix_element2(4*k-(3-j)),-1),1)
!!
     temp_exlog2(4*k-(3-j))=.true.
     temp_cpu_list(4*k-(3-j))=ibase+l1*2**isite1+l2*2**isite2
     temp_int_type(4*k-(3-j))=3
   enddo
enddo

allocate(cpu_list(1:num_int_bond),int_type(1:num_int_bond),matrix_element2(1:num_int_bond), &
         exlog2(1:num_int_bond), site_pair(2,1:num_int_bond),ibond_type_p2(1:num_int_bond), &
         bondwt_p2(1:num_int_bond))
call sort_list(num_int_bond,temp_cpu_list,sorted_list)
!!write(*,*) 'num_int_bond',num_int_bond
!!!!!!
k=0
pre_cpu=lme
max_local_loop=0
do i=1,num_int_bond
  cpu_list(i)=temp_cpu_list(sorted_list(i))
  int_type(i)=temp_int_type(sorted_list(i))
  matrix_element2(i)=temp_matrix_element2(sorted_list(i))
  exlog2(i)=temp_exlog2(sorted_list(i))
  if(cpu_list(i)<0) exlog2(i)=.false.
  site_pair(1,i)=temp_site_pair(1,sorted_list(i))
  site_pair(2,i)=temp_site_pair(2,sorted_list(i))
  bondwt_p2(i)=temp_bondwt_p2(sorted_list(i))
  ibond_type_p2(i)=temp_ibond_type_p2(sorted_list(i))
  if(cpu_list(i)==lme) max_local_loop=max_local_loop+1
  if((pre_cpu.ne.cpu_list(i)).and.(cpu_list(i))>=0)then 
    k=k+1
    pre_cpu = cpu_list(i)
  endif
enddo
deallocate(temp_cpu_list,temp_int_type,temp_matrix_element2,temp_exlog2,temp_site_pair,temp_bondwt_p2,temp_ibond_type_p2)
num_tot_ex_cpus=k
ex_loop_start=max_local_loop + 1
!!
write(*,*) 'Total number of interacting bonds:',num_int_bond
write(*,*) 'Total number of cpus for exchange data:',num_tot_ex_cpus
write(*,*) 'Maximum number for inter-node interactions:',max_local_loop
! num_int_bond - max_local_loop
!!
!!call get_thread_num(threads)
!!
allocate(tot_ex_cpus(0:num_tot_ex_cpus))
k=0
tot_ex_cpus(0)=lme
pre_cpu=tot_ex_cpus(0)
do i=1,num_int_bond
  if((pre_cpu.ne.cpu_list(i)).and.(cpu_list(i))>=0)then
    k=k+1
    tot_ex_cpus(k) = cpu_list(i)
    pre_cpu = tot_ex_cpus(k)
  endif
enddo
call barrier
!!!!!!
end subroutine make_cpu_list2


subroutine sort_list(data_size,data_list,sorted_list)
implicit none
integer, intent(in) :: data_size
integer, intent(in) :: data_list(1:data_size)
integer, intent(inout) :: sorted_list(1:data_size)
logical :: loop
logical, allocatable :: log1(:),log2(:)
integer :: i, start, cnt, s
allocate(log1(1:data_size),log2(1:data_size))
do i=1,data_size
  sorted_list(i)=i
  log1(i)=.true.
  log2(i)=.true.
enddo

start=1
cnt=1
do
  do i=1,data_size
    if(log2(i))then
      start=i
      s=data_list(start)
      exit
    endif
  enddo
  log1(:)=log2(:)
  do i=1,data_size
    if(log2(i))then
      if(data_list(i).ne.s) log1(i)=.false.
    endif
  enddo

  do i=1,data_size
    if(log1(i))then
      log2(i)=.false.
      sorted_list(cnt)=i
      cnt=cnt+1
    endif
  enddo

  loop=.true.
  do i=1,data_size
    loop=loop.and.(.not.log2(i))
  enddo
  if(loop) exit
  if(cnt>data_size) exit
enddo
deallocate(log1,log2)  
end subroutine sort_list


subroutine my_mpi_start()
implicit none
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NP,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,ME,ierr)
allocate(stat_send(MPI_STATUS_SIZE,NP), stat_recv(MPI_STATUS_SIZE,NP))
allocate(request_send(0:NP-1), request_recv(0:NP-1))
allocate(stat_send2(MPI_STATUS_SIZE,NP), stat_recv2(MPI_STATUS_SIZE,NP))
allocate(request_send2(0:NP-1), request_recv2(0:NP-1))
end subroutine my_mpi_start


subroutine bcast_sim_data(num_th,old_np,sim_tm,sn,sl,qs,qe,pm,typ,rnd)
implicit none
real(KIND(0d0)),intent(inout) :: pm
integer, intent(inout) :: num_th, old_np, qs, qe, sn, sl, rnd
real, intent(inout) :: sim_tm
integer :: idt(8), dat_num
character(LEN=2), intent(inout) :: typ
dat_num=ubound(idt,1)-lbound(idt,1)
idt(1)=num_th
idt(2)=old_np
idt(3)=qs
idt(4)=qe
idt(5)=sn
idt(6)=sl
idt(7)=system_mem_size
idt(8)=rnd
call MPI_BCAST(idt,dat_num,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(sim_tm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(pm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(typ,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
num_th=idt(1)
old_np=idt(2)
qs=idt(3)
qe=idt(4)
sn=idt(5)
sl=idt(6)
system_mem_size=idt(7)
rnd=idt(8)
end subroutine bcast_sim_data

subroutine my_mpi_split(i_color,w_parallel)
implicit none
integer, intent(in) :: i_color
logical, intent(in) :: w_parallel
KEY=0
COLOR=i_color
!if(w_parallel)then
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,COLOR,KEY,LOCAL_MPI_COMM_WORLD,ierr)
  call MPI_COMM_SIZE(LOCAL_MPI_COMM_WORLD,lnp,ierr)
  call MPI_COMM_RANK(LOCAL_MPI_COMM_WORLD,lme,ierr)
!else
!  lme=ME
!  lnp=NP
!  MPI_COMM_WORLD=LOCAL_MPI_COMM_WORLD
!endif
end subroutine my_mpi_split


subroutine mpi_get_maxint(x_in,x_out)
implicit none
integer :: x_in, x_out
call MPI_ALLREDUCE(x_in,x_out,1,MPI_INTEGER,MPI_MAX,LOCAL_MPI_COMM_WORLD,ierr)
end subroutine mpi_get_maxint

subroutine wait(itag)
implicit none
integer, intent(in) :: itag
integer :: istat
call MPI_WAIT(itag,istat,ierr)
end subroutine wait

subroutine get_buffer(i)
implicit none
integer :: i
allocate(buff(8*10816624))
call MPI_BUFFER_ATTACH(buff,8*10816624,ierr)
end subroutine get_buffer

subroutine release_buffer
implicit none
call MPI_BUFFER_DETACH(buff,8*10816624,ierr)
if(allocated(buff))deallocate(buff)
end subroutine release_buffer

subroutine mpi_isendrecv(target_cpu,data_num,x,y,itags,itagr)
implicit none
integer, intent(in) :: data_num, target_cpu, itagr, itags
complex(KIND(0d0)), intent(inout) :: x(:)
complex(KIND(0d0)), intent(inout) :: y(:)
call MPI_ISEND(x,data_num,MPI_DOUBLE_COMPLEX,target_cpu,itags,LOCAL_MPI_COMM_WORLD,request_send(ME),ierr)
call MPI_IRECV(y,data_num,MPI_DOUBLE_COMPLEX,target_cpu,itagr,LOCAL_MPI_COMM_WORLD,request_recv(ME),ierr)
end subroutine mpi_isendrecv

subroutine mpi_isendrecv_c(target_cpu,data_num,xr,yr,xi,yi,itags1,itagr1,itags2,itagr2)
implicit none
integer, intent(in) :: data_num, target_cpu, itagr1, itags1, itagr2, itags2
double precision, intent(inout) :: xr(:),xi(:)
double precision, intent(inout) :: yr(:),yi(:)
call MPI_ISEND(xr,data_num,MPI_DOUBLE_PRECISION,target_cpu,itags1,LOCAL_MPI_COMM_WORLD,request_send(ME), ierr)
call MPI_ISEND(xi,data_num,MPI_DOUBLE_PRECISION,target_cpu,itags2,LOCAL_MPI_COMM_WORLD,request_send2(ME),ierr)
call MPI_IRECV(yr,data_num,MPI_DOUBLE_PRECISION,target_cpu,itagr1,LOCAL_MPI_COMM_WORLD,request_recv(ME), ierr)
call MPI_IRECV(yi,data_num,MPI_DOUBLE_PRECISION,target_cpu,itagr2,LOCAL_MPI_COMM_WORLD,request_recv2(ME),ierr)
end subroutine mpi_isendrecv_c

subroutine mpi_test_isendrecv(send_logical,recv_logical)
implicit none
logical, intent(out) :: send_logical,recv_logical
call MPI_TEST(request_send(ME),send_logical,stat_send,ierr)
call MPI_TEST(request_recv(ME),recv_logical,stat_recv,ierr)
end subroutine mpi_test_isendrecv

subroutine mpi_test_isendrecv_c(send_logical1,recv_logical1,send_logical2,recv_logical2)
implicit none
logical, intent(out) :: send_logical1,recv_logical1,send_logical2,recv_logical2
call MPI_TEST(request_send(ME), send_logical1,stat_send, ierr)
call MPI_TEST(request_recv(ME), recv_logical1,stat_recv, ierr)
call MPI_TEST(request_send2(ME),send_logical2,stat_send2,ierr)
call MPI_TEST(request_recv2(ME),recv_logical2,stat_recv2,ierr)
end subroutine mpi_test_isendrecv_c

subroutine mpi_probe_data(target_cpu,itags,itagr,send_logical,recv_logical)
implicit none
integer, intent(in) :: target_cpu, itags, itagr
logical, intent(out) :: send_logical,recv_logical
!call MPI_IPROBE(Iam,itags,LOCAL_MPI_COMM_WORLD,send_logical,stat_send,ierr)
send_logical=.true.
call MPI_IPROBE(target_cpu,itagr,LOCAL_MPI_COMM_WORLD,recv_logical,stat_recv,ierr)
end subroutine mpi_probe_data


subroutine mpi_mywait
implicit none
call MPI_WAIT(request_send(ME),stat_send,ierr)
call MPI_WAIT(request_recv(ME),stat_recv,ierr)
!call MPI_WAITALL(1, request_recv(lme), stat_recv, ierr)
!call MPI_WAITALL(1, request_send(lme), stat_send, ierr)
end subroutine mpi_mywait

subroutine mpi_mywait_c
implicit none
call MPI_WAIT(request_send(ME), stat_send, ierr)
call MPI_WAIT(request_recv(ME), stat_recv, ierr)
call MPI_WAIT(request_send2(ME),stat_send2,ierr)
call MPI_WAIT(request_recv2(ME),stat_recv2,ierr)
!call MPI_WAITALL(1, request_recv(lme),  stat_recv,  ierr)
!call MPI_WAITALL(1, request_send(lme),  stat_send,  ierr)
!call MPI_WAITALL(1, request_recv2(lme), stat_recv2, ierr)
!call MPI_WAITALL(1, request_send2(lme), stat_send2, ierr)
end subroutine mpi_mywait_c

subroutine mpi_check_logical(loop)
implicit none
logical, intent(inout) :: loop
logical :: l_in, l_out
l_in = loop
call MPI_ALLREDUCE(l_in,l_out,1,MPI_LOGICAL,MPI_LAND,LOCAL_MPI_COMM_WORLD,ierr)
loop = l_out
end subroutine mpi_check_logical

subroutine allreduce_rdata(data_num,x_in,x_out)
implicit none
integer, intent(in) :: data_num
integer :: ierr1
double precision, intent(inout) :: x_in, x_out
call MPI_ALLREDUCE(x_in,x_out,data_num,MPI_DOUBLE_PRECISION,MPI_SUM,LOCAL_MPI_COMM_WORLD,ierr1)
end subroutine allreduce_rdata

subroutine allreduce_data(data_num,x_in,x_out)
implicit none
integer, intent(in) :: data_num
integer :: ierr1
complex(KIND(0d0)), intent(inout) :: x_in, x_out
call MPI_ALLREDUCE(x_in,x_out,data_num,MPI_DOUBLE_COMPLEX,MPI_SUM,LOCAL_MPI_COMM_WORLD,ierr1)
end subroutine allreduce_data

subroutine exchange_data(target_cpu,num_data,itag0,x,y)
implicit none
integer :: itag0
integer, intent(in) :: target_cpu,num_data
complex(KIND(0d0)), intent(in) :: x(:)
complex(KIND(0d0)), intent(inout) :: y(:)
integer :: ireq1, ireq2
integer :: ierr1, ierr2
call MPI_SENDRECV(x,num_data,MPI_DOUBLE_COMPLEX,target_cpu,itag0, &
                  y,num_data,MPI_DOUBLE_COMPLEX,target_cpu,itag0, &
                  LOCAL_MPI_COMM_WORLD,stat_send,ierr1)
end subroutine exchange_data


subroutine error_stop(char)
implicit none
character(LEN=256), intent(in) :: char
integer :: errcode
write(*,*) 'rank: ',ME
write(*,*) 'message: ',char
call MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
end subroutine error_stop

subroutine barrier
implicit none
call MPI_BARRIER(LOCAL_MPI_COMM_WORLD,ierr)
end subroutine barrier

subroutine all_barrier
implicit none
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine all_barrier



subroutine dealloc_mem_parallel
implicit none
if(allocated(cpu_pair_list))deallocate(cpu_pair_list)
if(allocated(lbits_pair))deallocate(lbits_pair)
if(allocated(bbits_pair))deallocate(bbits_pair)
if(allocated(hbits_pair))deallocate(hbits_pair)
if(allocated(ibond_type_p_h))deallocate(ibond_type_p_h)
if(allocated(ibond_type_p_l))deallocate(ibond_type_p_l)
if(allocated(ibond_type_p_b))deallocate(ibond_type_p_b)
if(allocated(list))deallocate(list)
if(allocated(exchange_list_b))deallocate(exchange_list_b)
if(allocated(exchange_list_h))deallocate(exchange_list_h)
if(allocated(bondwt_p_b))deallocate(bondwt_p_b)
if(allocated(bondwt_p_h))deallocate(bondwt_p_h)
if(allocated(matrix_element))deallocate(matrix_element)
if(allocated(cpu_list))deallocate(cpu_list)
if(allocated(int_type))deallocate(int_type)
if(allocated(matrix_element2))deallocate(matrix_element2)
if(allocated(exlog2))deallocate(exlog2)
if(allocated(site_pair))deallocate(site_pair)
if(allocated(stat_send))deallocate(stat_send)
if(allocated(stat_recv))deallocate(stat_recv)
if(allocated(request_send))deallocate(request_send)
if(allocated(request_recv))deallocate(request_recv)
if(allocated(stat_send2))deallocate(stat_send2)
if(allocated(stat_recv2))deallocate(stat_recv2)
if(allocated(request_send2))deallocate(request_send2)
if(allocated(request_recv2))deallocate(request_recv2)
end subroutine dealloc_mem_parallel

subroutine mpi_end
implicit none
call all_barrier
call MPI_FINALIZE(ierr)
end subroutine mpi_end
end module parallel_lib
