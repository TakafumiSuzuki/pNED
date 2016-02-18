module timer
!!!!
!!
!!
!!!!
implicit none
private
integer :: start(8)
integer :: checkpoint(0:8,100)
integer :: checkpoint_c(0:8,100)
integer :: v(8)
real :: y(8)
real :: mds(0:12)=(/0e0,31e0,28e0,31e0,30e0,31e0,30e0,31e0,31e0,30e0,31e0,30e0,31e0/)
!!                                             hour  min  sec msec  
real :: tunit(8)=(/3.1104e7,2.592e6,8.64e4,0e0,3.6e3,60e0,1e0,1e-3/)
logical :: t_stat(100)

public init_timer, timer_start, timer_stop, tell_me_time, make_timer_log

contains
real function get_secs(x)
integer :: i
real :: x(8)
real :: secs
secs=0.0
do i=8,5,-1
  secs=secs+x(i)*tunit(i)
enddo
get_secs=secs
end function get_secs

function tell_me_time()
real :: tell_me_time
real :: x(8), s1, s2
call date_and_time(values=v)
x(3:8)=float(v(3:8)-start(3:8))
s1 = ( x(3) + mds(v(2)-start(2)) )*tunit(3) 
s2 = get_secs(x)
tell_me_time = s1+s2
end function tell_me_time

subroutine init_timer
implicit none
real :: sumd, tp(0:11)
integer :: i, j
checkpoint(:,:)=0
checkpoint_c(:,:)=0
v(1:8)=0
t_stat(:)=.false.
call date_and_time(values=v)
start(1:8)=v(1:8)
!!
if((mod(start(1),4)==0).and.start(2)<2) mds(2)=mds(2)+1.d0
sumd=0.0
j=start(2)
do i=0,11
 tp(i)=sumd
 sumd=sumd+mds(j)
 j=j+1
 if(j>12) j=1
enddo
mds(0:11)=tp(0:11)
mds(12)=365e0
end subroutine init_timer

subroutine timer_start(i)
implicit none
integer, intent(in) :: i
if(.not.t_stat(i)) write(*,*) 'timer start (mem cleared):',i
  call date_and_time(values=v)
  checkpoint_c(5:8,i)=v(5:8)
  t_stat(i)=.true.
end subroutine timer_start

subroutine timer_stop(i)
implicit none
integer, intent(in) :: i
if(t_stat(i))then
  call date_and_time(values=v)
  checkpoint(5:8,i)=checkpoint(5:8,i)+(v(5:8)-checkpoint_c(5:8,i))
else
  write(*,*) 'timer error (timer_start is not called):',i
  t_stat(i)=.false.
endif
end subroutine timer_stop

subroutine make_timer_log(Iam)
implicit none
integer, intent(in) :: Iam
integer :: i
real :: t1
character(LEN=5) :: IDENTIFY
write(IDENTIFY,'(I5.5)') Iam
open(Iam+500,file='timer_log'//trim(IDENTIFY)//'.txt',status='unknown',access='append')
do i=1,15
    y(1:8)=float(checkpoint(1:8,i))
    t1=get_secs(y)
    write(Iam+500,'(A,I3.3,A,F18.4)') 'timer ',i,' :',t1
enddo
close(Iam+500)
end subroutine make_timer_log

end module timer
