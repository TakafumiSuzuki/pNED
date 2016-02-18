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
integer, allocatable, public :: real_cpu_list(:)
integer, allocatable, public :: lbits_pair(:,:)
integer, allocatable, public :: bbits_pair(:,:)
integer, allocatable, public :: hbits_pair(:,:)
integer, allocatable, public :: ibond_type_p_l(:)
integer, allocatable, public :: ibond_type_p_b(:)
integer, allocatable, public :: ibond_type_p_h(:)
integer, allocatable, public :: matrix_element(:,:)
integer, allocatable, public :: label_list(:)
logical, allocatable, public :: exlog(:,:)
double precision, allocatable, public :: bondwt_p_l(:)
double precision, allocatable, public :: bondwt_p_b(:)
double precision, allocatable, public :: bondwt_p_h(:)
!! for mltply v3
integer, public :: size_of_eff_mtx
integer, public :: num_int_bond
integer, public :: num_tot_ex_cpus
integer, allocatable, public :: cpu_list_send(:), cpu_list_recv(:)
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
       get_buffer,release_buffer, get_np, get_me, get_color, my_mpi_split, bcast_sim_data


contains
subroutine set_parallel_bits
use parameters, only : n, ldim, p_bit, b_bit, q_spin, size_of_local_matrix
implicit none
p_bit=nint(log10(dble(lnp))/log10(dble(q_spin)))
if((q_spin)**p_bit.ne.lnp)then
  write(*,*) 'error : wrong cpu number',p_bit,lnp
  stop
endif
b_bit=n-p_bit
ldim=(q_spin)**b_bit
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
use parameters, only : n, ibond, b_bit, p_bit, q_spin, size_of_local_matrix
use data_mem, only : ipair, bondwt, ibond_type, plocal_Ham
use openmp_lib
implicit none
integer :: i, j, k, isite1, isite2, start, num_b, num_h, l1, l2, istat,  ibase, num_l
integer :: pre_cpu, fshift, t_bit0, t_bit1, t_bit2, t_bit3, tind
integer :: cme, c_bbit, send_c, recv_c, send_c1, recv_c1, send_c2, recv_c2
integer :: p, num_hp, a1,a2
integer, allocatable :: temp_label_list(:),temp_ipair(:),sorted_list(:)

integer, allocatable :: temp_cpu_list_recv(:),temp_cpu_list_send(:),temp_real_cpu_list(:)
integer, allocatable :: temp_int_type(:),temp_matrix_element2(:),temp_site_pair(:,:)
integer, allocatable :: temp_ibond_type_p2(:)
integer, allocatable :: temp_list(:)
logical, allocatable :: temp_exlog2(:)
logical, allocatable :: l_ex_cpus(:)
double precision, allocatable :: temp_bondwt_p2(:)


allocate(lbits_pair(2,ibond), bbits_pair(2,ibond), hbits_pair(2,ibond))
allocate(ibond_type_p_l(ibond), ibond_type_p_b(ibond), ibond_type_p_h(ibond))
allocate(bondwt_p_l(ibond), bondwt_p_b(ibond), bondwt_p_h(ibond))
allocate(list(0:p_bit-1))

start=lme
do i=0,p_bit-1
  list(i)=mod(start,q_spin)
  start=start/q_spin
enddo
size_of_local_matrix=q_spin**2

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

!!
allocate(temp_list(ibond),sorted_list(ibond),temp_ipair(2*ibond), &
         temp_bondwt_p2(1:ibond),temp_ibond_type_p2(1:ibond))
do i=1,ibond
  isite1=ipair(2*i-1)
  isite2=ipair(2*i  )
  temp_ipair(2*i-1)=min(isite1,isite2)
  temp_ipair(2*i  )=max(isite1,isite2)
  temp_list(i)=temp_ipair(2*i)
  temp_bondwt_p2(i)=bondwt(i)
  temp_ibond_type_p2(i)=ibond_type(i)
enddo
!!
call sort_list(ibond,temp_list,sorted_list)
!!
do i=1,ibond
  ipair(2*i-1)=temp_ipair(2*sorted_list(i)-1)
  ipair(2*i  )=temp_ipair(2*sorted_list(i)  )
  bondwt(i)=temp_bondwt_p2(sorted_list(i))
  ibond_type(i)=temp_ibond_type_p2(sorted_list(i))
enddo
deallocate(temp_list,temp_ipair,sorted_list,temp_bondwt_p2,temp_ibond_type_p2)
!!
!!
!!
num_hp=0
do i=1,ibond
  isite1=ipair(2*i-1)-1
  isite2=ipair(2*i  )-1
  if((isite1<b_bit).and.(isite2<b_bit))then
    num_lbits_p=num_lbits_p+1
    lbits_pair(1,num_lbits_p)=isite1
    lbits_pair(2,num_lbits_p)=isite2
    bondwt_p_l(num_lbits_p)=bondwt(i)
    ibond_type_p_l(num_lbits_p)=ibond_type(i)
!!    write(*,'(3I4)') isite1,isite2,ibond_type(i)
  elseif(((isite1<b_bit).and.(isite2>=b_bit)).or.((isite1>=b_bit).and.(isite2<b_bit)))then
    num_bbits_p=num_bbits_p+1  
    bbits_pair(1,num_bbits_p)=isite1
    bbits_pair(2,num_bbits_p)=isite2 - b_bit
    bondwt_p_b(num_bbits_p)=bondwt(i)
    ibond_type_p_b(num_bbits_p)=ibond_type(i)

  elseif((isite1>=b_bit).and.(isite2>=b_bit))then
    num_hbits_p=num_hbits_p+1  
    hbits_pair(1,num_hbits_p)=isite1 - b_bit
    hbits_pair(2,num_hbits_p)=isite2 - b_bit
    bondwt_p_h(num_hbits_p)=bondwt(i)
    ibond_type_p_h(num_hbits_p)=ibond_type(i)
    do k=0,size_of_local_matrix-1
      if(plocal_Ham(k,ibond_type(i))) num_hp=num_hp+1
    enddo
  else
    write(*,*) 'error: #=', i,isite1,isite2
  endif
enddo
!!
!!
num_l=max(num_lbits_p,1)
num_b=max(num_bbits_p,1)
num_h=max(num_hbits_p,1)
num_int_bond=size_of_local_matrix*(num_lbits_p+num_bbits_p)+num_hp
fshift=0
!!
!!
allocate(temp_cpu_list_recv(1:num_int_bond),&
         temp_cpu_list_send(1:num_int_bond),&
         temp_real_cpu_list(1:num_int_bond), &
         temp_int_type(1:num_int_bond),&
         temp_matrix_element2(1:num_int_bond),&
         temp_exlog2(1:num_int_bond),temp_site_pair(2,1:num_int_bond), &
         temp_bondwt_p2(1:num_int_bond),temp_ibond_type_p2(1:num_int_bond), &
         sorted_list(1:num_int_bond))
temp_exlog2(:)=.false.
!!
!!
tind=1
do i=0,num_lbits_p-1
  do k=0,size_of_local_matrix-1
    temp_cpu_list_send(tind)=lme
    temp_cpu_list_recv(tind)=lme
    temp_real_cpu_list(tind)=3*(lme*lnp+lme)
    temp_int_type(tind)=1
    temp_matrix_element2(tind)=mod(k,size_of_local_matrix)
    temp_exlog2(tind)=.true.
    temp_site_pair(1,tind)=lbits_pair(1,i+1)
    temp_site_pair(2,tind)=lbits_pair(2,i+1)
    temp_bondwt_p2(tind)=bondwt_p_l(i+1)
    temp_ibond_type_p2(tind)=ibond_type_p_l(i+1)
!!    write(*,'(4I6)') temp_site_pair(1,tind),temp_site_pair(2,tind),temp_ibond_type_p2(tind)
    tind=tind+1
  enddo
enddo
!!
do j=0,num_bbits_p-1
  i=j+num_lbits_p
  isite2=bbits_pair(2,j+1)
  cme = lme - list(isite2)*(q_spin**isite2)
  t_bit0=list(isite2)
!!  write(*,'(A,3I6)') 'PMe,isite2,bit',lme,isite2,list(isite2)
  do k=0,size_of_local_matrix-1
    t_bit2=mod(k+fshift,size_of_local_matrix)/q_spin
    temp_cpu_list_send(tind)=cme+mod(q_spin-1 -t_bit2 +t_bit0, q_spin)*(q_spin**isite2)
    temp_cpu_list_recv(tind)=cme+mod(1        +t_bit2 +t_bit0, q_spin)*(q_spin**isite2)
    t_bit3=1
    if(temp_cpu_list_recv(tind)*lnp+temp_cpu_list_send(tind)==lme*lnp+lme) t_bit3=0
    temp_real_cpu_list(tind)=3*(temp_cpu_list_recv(tind)*lnp+temp_cpu_list_send(tind))+t_bit3
    temp_int_type(tind)=2

    temp_matrix_element2(tind)=mod(mod(list(isite2)+1,q_spin)*q_spin+k,size_of_local_matrix)

    temp_exlog2(tind)=.true.
    temp_site_pair(1,tind)=bbits_pair(1,j+1)
    temp_site_pair(2,tind)=bbits_pair(2,j+1)
    temp_bondwt_p2(tind)=bondwt_p_b(j+1)
    temp_ibond_type_p2(tind)=ibond_type_p_b(j+1)
!    write(*,'(A,3I6,4I5,2I5)')'PMe,i,k,s,r,R,m',lme,i,tind,temp_cpu_list_send(tind),temp_cpu_list_recv(tind),&
!          temp_real_cpu_list(tind),temp_matrix_element2(tind),temp_site_pair(1,tind),temp_site_pair(2,tind)
    tind=tind+1
  enddo
enddo
!!
do j=0,num_hbits_p-1
  i=j+num_lbits_p+num_bbits_p
  isite1=hbits_pair(1,j+1)
  isite2=hbits_pair(2,j+1)
  cme =lme - list(isite1)*(q_spin**isite1) - list(isite2)*(q_spin**isite2)
  t_bit0 = list(isite2)*q_spin + list(isite1)
  p=ibond_type_p_h(j+1)
  do k=0,size_of_local_matrix-1
    if(plocal_Ham(k,p))then
      ibase  =mod(-k +t_bit0+size_of_local_matrix, size_of_local_matrix)

      send_c =mod( k +t_bit0+size_of_local_matrix, size_of_local_matrix)
      recv_c =mod(-k +t_bit0+size_of_local_matrix, size_of_local_matrix)

      send_c1=mod(send_c,q_spin); send_c2=send_c/q_spin
      recv_c1=mod(recv_c,q_spin); recv_c2=recv_c/q_spin
      temp_cpu_list_send(tind)=cme +send_c1*q_spin**isite1 +send_c2*q_spin**isite2
      temp_cpu_list_recv(tind)=cme +recv_c1*q_spin**isite1 +recv_c2*q_spin**isite2
      t_bit3=2
      if(temp_cpu_list_recv(tind)*lnp+temp_cpu_list_send(tind)==lme*lnp+lme) t_bit3=0
      temp_real_cpu_list(tind)=3*(temp_cpu_list_recv(tind)*lnp+temp_cpu_list_send(tind))+t_bit3
      temp_int_type(tind)=3

      temp_matrix_element2(tind)=ibase

      temp_exlog2(tind)=.true.
      temp_site_pair(1,tind)=hbits_pair(1,j+1)
      temp_site_pair(2,tind)=hbits_pair(2,j+1)
      temp_bondwt_p2(tind)=bondwt_p_h(j+1)
      temp_ibond_type_p2(tind)=ibond_type_p_h(j+1)
!    write(*,'(A,3I6,4I5,2I5)')'Me,i,k,s,r,R,m',lme,i,tind,temp_cpu_list_send(tind),temp_cpu_list_recv(tind),&
!          temp_real_cpu_list(tind),temp_matrix_element2(tind),temp_site_pair(1,tind),temp_site_pair(2,tind)
      tind=tind+1
    endif
  enddo
enddo
if(tind-1.ne.num_int_bond)then
  write(*,*) 'bond error'
  stop
endif
!!
!!
call sort_list2(num_int_bond,temp_real_cpu_list,sorted_list)
call MPI_BCAST(sorted_list,num_int_bond,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!
!!
allocate(cpu_list_send(1:num_int_bond),cpu_list_recv(1:num_int_bond),real_cpu_list(1:num_int_bond), &
         int_type(1:num_int_bond),matrix_element2(1:num_int_bond), &
         exlog2(1:num_int_bond), site_pair(2,1:num_int_bond),ibond_type_p2(1:num_int_bond), &
         bondwt_p2(1:num_int_bond))
!!
!do i=1,num_int_bond
!    write(*,'(A,2I6,2I5,2I6)')'Me,i,k,s,r,R',lme,i,temp_cpu_list_send(i),temp_cpu_list_recv(i),&
!         temp_real_cpu_list(i),temp_matrix_element2(i)
!enddo
!!!!!!
!!!!!!
allocate(l_ex_cpus(0:lnp-1))
l_ex_cpus(:)=.true.
l_ex_cpus(lme)=.false.
k=0
pre_cpu=lme*lnp+lme
max_local_loop=0
do i=1,num_int_bond
  cpu_list_recv(i)=temp_cpu_list_recv(sorted_list(i))
  cpu_list_send(i)=temp_cpu_list_send(sorted_list(i))
  real_cpu_list(i)=temp_real_cpu_list(sorted_list(i))/3
  int_type(i)=temp_int_type(sorted_list(i))
  matrix_element2(i)=temp_matrix_element2(sorted_list(i))
  exlog2(i)=temp_exlog2(sorted_list(i))
  if(real_cpu_list(i)<0) exlog2(i)=.false.
  site_pair(1,i)=temp_site_pair(1,sorted_list(i))
  site_pair(2,i)=temp_site_pair(2,sorted_list(i))
  bondwt_p2(i)=temp_bondwt_p2(sorted_list(i))
  ibond_type_p2(i)=temp_ibond_type_p2(sorted_list(i))
!!  write(*,'(4I6)') site_pair(1,i),site_pair(2,i),ibond_type_p2(i),matrix_element2(i)
!!
  if(real_cpu_list(i)==lme*lnp+lme) max_local_loop=max_local_loop+1
!  if((pre_cpu.ne.real_cpu_list(i)).and.(real_cpu_list(i)>=0))then
!    pre_cpu = real_cpu_list(i)
    if(l_ex_cpus(cpu_list_recv(i)))then 
      k=k+1
      l_ex_cpus(cpu_list_recv(i))=.false.
    endif
!  endif
enddo
deallocate(temp_cpu_list_recv,temp_cpu_list_send,temp_real_cpu_list,temp_int_type,sorted_list, &
           temp_matrix_element2,temp_exlog2,temp_site_pair,temp_bondwt_p2,temp_ibond_type_p2)
num_tot_ex_cpus=k
ex_loop_start=max_local_loop + 1
!!
!!
!do i=1,num_int_bond
!    write(*,'(A,2I6,2I5,2I6)')'Me,i,k,s,r,R',lme,i,cpu_list_send(i),cpu_list_recv(i),&
!         real_cpu_list(i),matrix_element2(i)
!enddo
!!
!!
! num_int_bond - max_local_loop
!!
allocate(tot_ex_cpus(0:num_tot_ex_cpus))
k=0
l_ex_cpus(:)=.true.
l_ex_cpus(lme)=.false.
tot_ex_cpus(:)=lme
pre_cpu=lme*lnp+lme
do i=1,num_int_bond
!  if((pre_cpu.ne.real_cpu_list(i)).and.(real_cpu_list(i)>=0))then
!    pre_cpu = real_cpu_list(i)
    if(l_ex_cpus(cpu_list_recv(i)))then
      k=k+1
      tot_ex_cpus(k) = cpu_list_recv(i)
      l_ex_cpus(cpu_list_recv(i))=.false.
    endif
!  endif
!  if(i>=ex_loop_start)  write(*,'(A,2I6,2I5,2I4,2I4,I6)')'Me,i,s,r,R,m',lme,i,cpu_list_send(i),cpu_list_recv(i),&
!          real_cpu_list(i),matrix_element2(i),site_pair(1,i),site_pair(2,i),int_type(i)
enddo
deallocate(l_ex_cpus)
!!
allocate(label_list(1:num_int_bond))
!!
!!
label_list(:)=0
do i=ex_loop_start,num_int_bond
  if(cpu_list_recv(i).ne.lme)then
    do k=1,num_tot_ex_cpus
      if(cpu_list_recv(i)==tot_ex_cpus(k))then
        label_list(i)=k
      endif
    enddo
  endif
!  write(*,'(A,4I6)') 'list i me label',i,lme,label_list(i),tot_ex_cpus(label_list(i))
enddo
!!!!
write(*,*) 'Total number of interacting bonds:',num_int_bond
write(*,*) 'Total number of cpus for exchange data:',num_tot_ex_cpus
write(*,*) 'Maximum number for inter-node interactions:',max_local_loop
!!!!
!!
end subroutine make_cpu_list2


subroutine sort_list(data_size,data_list,sorted_list)
implicit none
integer, intent(in) :: data_size
integer, intent(in) :: data_list(1:data_size)
integer, intent(inout) :: sorted_list(1:data_size)
integer, allocatable :: temp(:)
logical :: loop
integer :: i, j, k
allocate(temp(data_size))
do i=1,data_size
  sorted_list(i)=i
  temp(i)=data_list(i)
enddo

do i=1,data_size-1
  do j=i+1,data_size
    if (temp(i)>temp(j)) then
      k=temp(i)
      temp(i)=temp(j)
      temp(j)=k
      k=sorted_list(i)
      sorted_list(i)=sorted_list(j)
      sorted_list(j)=k
    endif
  enddo
enddo
deallocate(temp)  
end subroutine sort_list

subroutine sort_list2(data_size,data_list,sorted_list)
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
end subroutine sort_list2


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


subroutine bcast_sim_data(num_th,old_np,sim_tm,sn,sl,qs,qe,pm,typ,sp,nrnd)
implicit none
real(KIND(0d0)),intent(inout) :: pm
integer, intent(inout) :: num_th, old_np, qs, qe, sn, sl, sp,nrnd
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
idt(7)=sp
idt(8)=nrnd
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
sp=idt(7)
nrnd=idt(8)
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
integer, intent(in) :: data_num, target_cpu(2), itagr, itags
complex(KIND(0d0)), intent(inout) :: x(:)
complex(KIND(0d0)), intent(inout) :: y(:)
call MPI_ISEND(x,data_num,MPI_DOUBLE_COMPLEX,target_cpu(1),itags,LOCAL_MPI_COMM_WORLD,request_send(ME),ierr)
call MPI_IRECV(y,data_num,MPI_DOUBLE_COMPLEX,target_cpu(2),itagr,LOCAL_MPI_COMM_WORLD,request_recv(ME),ierr)
end subroutine mpi_isendrecv

subroutine mpi_isendrecv_c(target_cpu,data_num,xr,yr,xi,yi,itags1,itagr1,itags2,itagr2)
implicit none
integer, intent(in) :: data_num, target_cpu(2), itagr1, itags1, itagr2, itags2
double precision, intent(inout) :: xr(:),xi(:)
double precision, intent(inout) :: yr(:),yi(:)
call MPI_ISEND(xr,data_num,MPI_DOUBLE_PRECISION,target_cpu(1),itags1,LOCAL_MPI_COMM_WORLD,request_send(ME), ierr)
call MPI_ISEND(xi,data_num,MPI_DOUBLE_PRECISION,target_cpu(1),itags2,LOCAL_MPI_COMM_WORLD,request_send2(ME),ierr)
call MPI_IRECV(yr,data_num,MPI_DOUBLE_PRECISION,target_cpu(2),itagr1,LOCAL_MPI_COMM_WORLD,request_recv(ME), ierr)
call MPI_IRECV(yi,data_num,MPI_DOUBLE_PRECISION,target_cpu(2),itagr2,LOCAL_MPI_COMM_WORLD,request_recv2(ME),ierr)
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
complex(KIND(0d0)), intent(inout) :: x(:)
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
if(allocated(real_cpu_list))deallocate(real_cpu_list)
if(allocated(lbits_pair))deallocate(lbits_pair)
if(allocated(bbits_pair))deallocate(bbits_pair)
if(allocated(hbits_pair))deallocate(hbits_pair)
if(allocated(ibond_type_p_h))deallocate(ibond_type_p_h)
if(allocated(ibond_type_p_l))deallocate(ibond_type_p_l)
if(allocated(ibond_type_p_b))deallocate(ibond_type_p_b)
if(allocated(list))deallocate(list)
if(allocated(bondwt_p_b))deallocate(bondwt_p_b)
if(allocated(bondwt_p_h))deallocate(bondwt_p_h)
if(allocated(matrix_element))deallocate(matrix_element)
if(allocated(real_cpu_list))deallocate(real_cpu_list)
if(allocated(tot_ex_cpus))deallocate(tot_ex_cpus)
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
if(allocated(label_list))deallocate(label_list)
end subroutine dealloc_mem_parallel

subroutine mpi_end
implicit none
call all_barrier
call MPI_FINALIZE(ierr)
end subroutine mpi_end
end module parallel_lib
