!!
module check_results
implicit none
private
public check1
!!
contains
subroutine check1(x)
use parameters, only : n, ldim
use parallel_lib, only : lme, get_me
use parallelized_diag, only : p_mltply_v3, calc_vec_inprdct
use working_area, only : wwk
implicit none
complex(KIND(0d0)), intent(inout) :: x(:)
complex(KIND(0d0)) :: cprd, cdnorm, cdmm
double precision :: Hexpec
integer :: i, pme
character(LEN=64) :: ACC_FNAME
pme=get_me()
write(ACC_FNAME,'(A,I5.5,A)') 'acc_log_me',pme,'.txt'
!!
cdnorm=0d0
call calc_vec_inprdct(ldim,x,x,cdnorm)
!!
if(CDABS(cdnorm).lt.1.d-30)then
         print *,' #(W09)# Null vector given to check1'
         return
endif
!!
!$OMP PARALLEL
!$OMP WORKSHARE
wwk(:,1)=0d0
!$OMP END WORKSHARE
!$OMP END PARALLEL
call p_mltply_v3(x,wwk(:,1),cdmm,0d0)
!
cprd=0.0d0
call calc_vec_inprdct(ldim,wwk(:,1),x,cprd)

Hexpec=DREAL(cprd)
open(23,file=trim(ACC_FNAME),status='unknown',access='append')
write(23,*)
write(23,200)
 200  format(' ---------------------------- Information from check1')
write(23,210) cprd
 210  format(' <x*H*x> =(',1pe16.8,',',1pe16.8,')')
write(23,220)
 220  format(' H*x(j)/x(j) (j=min(ldim/3,13),ldim,max(1,ldim/20))')
write(23,230) (wwk(i,1)/x(i),i=min(ldim/3,13),ldim,max(1,ldim/20))
 230  format(4('(',1pe16.6e3,',',1pe16.6e3,')'))
write(23,240)
 240  format(' ---------------------------------------------------')
close(23)
end subroutine check1
end module check_results
