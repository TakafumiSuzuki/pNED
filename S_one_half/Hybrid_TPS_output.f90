module def_output
implicit none
private
character(LEN=128), save :: DUMP_NAME

public readwrite_vec_data, write_dsf_res, open_res_files, close_res_files, &
       interval_dump_read,interval_dump,check_interval_dump, set_dump_name



contains
subroutine readwrite_vec_data(ind,type,ll,e,x)
use parallel_lib, only : barrier
implicit none
integer, intent(in) :: ll, ind
logical, intent(in) :: type
double precision, intent(inout) :: e
complex(KIND(0d0)), intent(inout) :: x(1:ll)
character(LEN=32) :: IDENTIFY
integer :: j
write(IDENTIFY,'(A,I5.5)') 'me',ind
j=ind+100
!!
open(j,file="basis"//trim(IDENTIFY),form='unformatted',access='sequential',status='unknown')
if(type)then
  write(*,*) 'dumped vec data'
  write(j) e
  write(j) x
else
  write(*,*) 'reading vec data'
  read(j) e
  read(j) x
  call barrier
endif
close(j)
end subroutine readwrite_vec_data
!!
!!
subroutine set_dump_name
use parallel_lib, only : get_me
implicit none
integer :: pme
character(LEN=128) :: IDENTIFY
pme=get_me()
write(IDENTIFY,'(a,i5.5)') 'vecdump_me',pme
DUMP_NAME=trim(IDENTIFY)
end subroutine set_dump_name
!!
!!
subroutine interval_dump_read(position,data_size,vec_size,sa,sb,sc,se,vec0,vec1,vec2,vec3,vec4)
use parameters
use parallel_lib, only : get_me, get_np, barrier
implicit none
logical, intent(inout) :: position(:)
integer, intent(inout) :: data_size, vec_size
double precision, intent(inout), optional :: sa(:)
double precision, intent(inout), optional :: sb(:)
double precision, intent(inout), optional :: sc(:,:)
double precision, intent(inout), optional :: se(:)
complex(KIND(0d0)), intent(inout), optional :: vec0(:)
complex(KIND(0d0)), intent(inout) :: vec1(:)
complex(KIND(0d0)), intent(inout) :: vec2(:)
complex(KIND(0d0)), intent(inout), optional :: vec3(:)
complex(KIND(0d0)), intent(inout), optional :: vec4(:)
character(LEN=128) :: IDENTIFY
integer :: i, j, k, m, pme, pnp
pme=get_me()
pnp=get_np()
j=100+pme+1
do i=0,pnp-1
if(i==pme)then
write(*,*) 'start: reading files',i
open(j,file=trim(DUMP_NAME),form='unformatted',access='sequential',status='old')
read(j) data_size
read(j) vec_size
if(present(sa)) read(j) sa
if(present(sb)) read(j) sb
if(present(sc)) read(j) sc
if(present(se)) read(j) se
if(present(vec0)) read(j) vec0
read(j) vec1
read(j) vec2
if(present(vec3)) read(j) vec3
if(present(vec4)) read(j) vec4
close(j)
write(*,*) 'end: reading files',i
endif
call barrier
enddo
!!
end subroutine interval_dump_read
!!
subroutine check_interval_dump(position)
use parameters
use parallel_lib, only : get_me, get_np
implicit none
logical, intent(out) :: position(:)
character(LEN=128) :: EXTENTION
integer :: i, j, k, pme, pnp
logical :: exec=.true.
integer :: t_np, t_max_itr, t_s_maxitr, t_n, t_ibond, t_nvec, t_total_sz  
character(LEN=2) :: t_dsf_comp
double precision :: t_kitaev, t_xy, t_xxz
pme=get_me()
pnp=get_np()
EXTENTION='.info'
j=100+pme+1
open(j,file=trim(DUMP_NAME)//trim(EXTENTION),form='unformatted',access='sequential',status='old', err=50)
read(j)  DUMP_NAME
read(j)  t_np
read(j)  position
read(j)  t_max_itr
read(j)  t_dsf_comp
read(j)  t_s_maxitr
read(j)  t_n
read(j)  t_ibond
read(j)  t_nvec
read(j)  t_total_sz
!!
write(*,*) 'Report of progress (current itr numbers)'
do i=1,12
  read(j) k
  itr_counter(i) = max(k,itr_counter(i))
  write(*,*) itr_counter(i)
enddo
close(j)
!if(t_np.ne.pnp) exec=exec .and. .false.
!if(t_max_itr.ne.max_itr) exec=exec .and. .false.
if(t_dsf_comp.ne.t_dsf_comp) exec=exec .and. .false.
!if(t_s_maxitr.ne.s_maxitr) exec=exec .and. .false.
if(t_n.ne.n) exec=exec .and. .false.
if(t_ibond.ne.ibond) exec=exec .and. .false.
!!if(t_nvec.ne.nvec) exec=exec .and. .false.
if(t_total_sz.ne.total_sz) exec=exec .and. .false.
if(.not. exec)then
  write(*,*) 'Error : The given parameter set is different.'
  write(*,'(A,3I5)') 'N:np,max_itr,s_maxitr', t_np,t_max_itr,t_s_maxitr
  write(*,'(A,4I5)') 'N:n ,ibond,nvec,total_sz', t_n,t_ibond,t_nvec,t_total_sz
  write(*,'(A,3I5)') 'O:max_itr,s_maxitr', max_itr,s_maxitr
  write(*,'(A,4I5)') 'O:n ,ibond,nvec,total_sz', n,ibond,nvec,total_sz
  stop
endif
 50 write(*,*) 'Read info files...'
end subroutine check_interval_dump
!!
!!
subroutine interval_dump(position,data_size,vec_size,sa,sb,sc,se,vec0,vec1,vec2,vec3,vec4)
use parameters
use data_mem
use working_area
use parallel_lib
use timer
implicit none
logical, intent(in) :: position(:)
integer, intent(in) :: data_size, vec_size
double precision, intent(in), optional :: sa(:)
double precision, intent(in), optional :: sb(:)
double precision, intent(in), optional :: sc(:,:)
double precision, intent(in), optional :: se(:)
complex(KIND(0d0)), intent(in), optional :: vec0(:)
complex(KIND(0d0)), intent(in) :: vec1(:)
complex(KIND(0d0)), intent(in) :: vec2(:)
complex(KIND(0d0)), intent(in), optional :: vec3(:)
complex(KIND(0d0)), intent(in), optional :: vec4(:)
character(LEN=128) :: EXTENTION, FILENAME
integer :: i, j, pme, pnp
pme=get_me()
pnp=get_np()
EXTENTION='.info'
FILENAME=trim(DUMP_NAME)//trim(EXTENTION)
write(*,*) FILENAME
j=100+pme+1
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
write(j)  DUMP_NAME
write(j)  lnp
write(j)  position
write(j)  max_itr
write(j)  dsf_comp
write(j)  s_maxitr
write(j)  n
write(j)  ibond
write(j)  nvec
write(j)  total_sz
do i=1,12
  write(j)  itr_counter(i)
enddo
close(j)
!!
open(j,file=trim(DUMP_NAME),form='unformatted',access='sequential',status='unknown')
write(j) data_size
write(j) vec_size
if(present(sa)) write(j) sa
if(present(sb)) write(j) sb
if(present(sc)) write(j) sc
if(present(se)) write(j) se
if(present(vec0)) write(j) vec0
write(j) vec1
write(j) vec2
if(present(vec3)) write(j) vec3
if(present(vec4)) write(j) vec4
close(j)
call barrier
!!
call dealloc_mem
call dealloc_mem_2
call dealloc_mem_parallel
call timer_stop(10)
call make_timer_log(pme)
call mpi_end
stop
end subroutine interval_dump
!!
!!
!!
!!
!!
subroutine open_res_files(ind,comp,i_color)
implicit none
integer, intent(in) :: ind,i_color
character(LEN=2), intent(in) :: comp
character(LEN=64) :: IDENTIFY
integer :: j
if(ind==0)then
  j=ind+100*i_color
  write(IDENTIFY,'(A,I5.5,A,I5.5)') 'cl',i_color,'lme',ind
  if(comp=='zz')then
    open(j,file="dsf_res_zz.dat"//trim(IDENTIFY),status='unknown')
  elseif(comp=='xx')then
    open(j,file="dsf_res_xx.dat"//trim(IDENTIFY),status='unknown')
  else
    open(j,file="dsf_comp_error.dat"//trim(IDENTIFY),status='unknown')
  endif
endif
end subroutine open_res_files
!!
subroutine close_res_files(ind,i_color)
implicit none
integer, intent(in) :: ind, i_color
integer :: j
if(ind==0)then
  j=ind+100*i_color
  close(j)
endif
end subroutine close_res_files
!!
subroutine write_dsf_res(ind,i_color,comp,ll,sq,q,w,dsf,conv)
implicit none
complex(KIND(0d0)), intent(in) :: sq
double precision, intent(in) :: q(:), dsf(1:ll), w(1:ll)
integer, intent(in) :: ind, i_color, ll, conv(1:ll)
character(LEN=2), intent(in) :: comp
integer :: j, k
if(ind==0)then
  j=ind+100*i_color
  write(j,'(A,F25.18)') '#        ', DREAL(sq)
  do k=1,ll
    write(j,'(2F18.10,2F25.16,I5)') q(1),q(2),w(k),dsf(k),conv(k)
  enddo
  write(j,*) ''
  write(j,*) ''
  write(j,*) ''
  write(j,*) ''
endif
end subroutine write_dsf_res
end module def_output
