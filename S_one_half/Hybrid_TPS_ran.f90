module rand_num1
implicit none
integer :: seedsize
integer,allocatable :: seed(:)

contains
subroutine int_ran(indx)
implicit none
integer, intent(in) :: indx
integer :: j
character(LEN=32) :: IDENTIFY
call random_seed(size=seedsize)
allocate(seed(seedsize))
!call random_seed(get=seed)
seed(:)=indx
call random_seed(put=seed)
write(IDENTIFY,'(A,I5.5)') 'me',indx
j=indx+100
open(j,file="seed"//trim(IDENTIFY),form='unformatted',access='sequential',status='unknown')
write(j) seedsize
write(j) seed
close(j)
end subroutine int_ran
!!
double precision function grnd()
double precision :: x(1)
call random_number(x)
grnd=x(1)
end function grnd
!!
end module rand_num1
!!
!!
!!
!!
!!
module rand_num2
      implicit none
      integer, private, parameter :: N     =  624
      integer, private, parameter :: N1    =  N+1
      integer, private, parameter :: M     =  397
      integer, private, parameter :: MATA  = -1727483681
      integer, private, parameter :: UMASK = -2147483648
      integer, private, parameter :: LMASK =  2147483647
      integer, private, parameter :: TMASKB= -1658038656
      integer, private, parameter :: TMASKC= -272236544
      integer :: mti
      integer :: mt(0:N-1)
      integer :: mag01(0:1)
      data mti/N1/
      data mag01/0, MATA/
!!
      contains
      subroutine int_ran(ixx)
        implicit none
        integer,intent(in) :: ixx 
        call sgrnd(ixx)
      end subroutine int_ran

      subroutine sgrnd(seed)
        implicit none
        integer,intent(in) :: seed
        mt(0)= iand(seed,-1)
        do mti=1,N-1
          mt(mti) = iand(69069 * mt(mti-1),-1)
        enddo
      end subroutine sgrnd
!!
      integer function int_grnd()
      implicit none
      integer :: y,kk
      if(mti.ge.N) then
        if(mti.eq.N+1) then
          call sgrnd(4357)
        endif
        do kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        do kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
        enddo
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
!      if(y.lt.0) then
!        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
!      else
!        grnd=dble(y)/(2.0d0**32-1.0d0)
!      endif
      int_grnd=y
      end function
!!
      function TSHFTU(y)
      implicit none
      integer,intent(in) :: y
      integer :: TSHFTU
           TSHFTU=ishft(y,-11)
      end function TSHFTU
!!
      function TSHFTS(y)
      implicit none
      integer,intent(in) :: y
      integer :: TSHFTS
           TSHFTS=ishft(y,7)
      end function TSHFTS
!!
      function TSHFTT(y)
      implicit none
      integer,intent(in) :: y
      integer :: TSHFTT
           TSHFTT=ishft(y,15)
      end function TSHFTT
!!
      function TSHFTL(y)
      implicit none
      integer,intent(in) :: y
      integer :: TSHFTL
           TSHFTL=ishft(y,-18)
      end function TSHFTL
!!
      double precision function grnd()
      implicit none
      integer :: int_g
      int_g=int_grnd()
      if(int_g<0)then
        grnd=(dble(int_g)+2.d0**32)/(2.d0**32-1.d0)
      else
        grnd=dble(int_g)/(2.d0**32-1.d0)
      endif
      end function
end module rand_num2
!!
!!
!!
!!
!!
!include 'mkl_vsl.fi'
!module rand_num3
!use MKL_VSL_TYPE
!use MKL_VSL
!implicit none
!TYPE(VSL_STREAM_STATE), save :: stream
!integer :: error_code
!integer :: method=VSL_METHOD_DUNIFORM_STD
!integer :: brng=VSL_BRNG_MT19937
!
!contains
!subroutine int_ran
!error_code = vslnewstream( stream, brng, 4357 )
!end subroutine int_ran
!
!double precision function grnd()
!double precision :: x(1)
!error_code = vdrnguniform(method,stream,1,x,0d0,1d0)
!grnd = x(1)
!end function grnd
!
!end module rand_num3
