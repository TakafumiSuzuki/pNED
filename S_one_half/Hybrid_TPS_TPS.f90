module parallelized_TPS_method
implicit none
private
public parallelized_TPS, init_res_dat, dealloc_mem_TPS
integer, parameter :: num_res_dat=30
integer, parameter :: ipm=1
integer :: local_tdim_l, sq_num
double precision, allocatable :: sq(:)
complex(KIND(0d0)), allocatable :: rt(:)
complex(KIND(0d0)), allocatable :: rk(:)
complex(KIND(0d0)), allocatable :: rr(:)
double precision :: cz(2)
!!
!!
contains
subroutine init_res_dat(Iam,tot_sz)
use parameters, only : itr_counter, max_ldim
use working_area, only : coef
use data_mem, only : q_list_symmetric
implicit none
integer, intent(in) :: tot_sz, Iam
integer :: i,lb,ub
if(allocated(coef))then
  deallocate(coef)
endif
allocate(coef(itr_counter(2),num_res_dat))
!!
sq_num=0
lb=lbound(q_list_symmetric,1)
ub=ubound(q_list_symmetric,1)
do i=lb, ub
  if(q_list_symmetric(i)>0) sq_num=sq_num+1
enddo
!!
if(.not.allocated(sq)) allocate(sq(4*sq_num))
!!
end subroutine init_res_dat

subroutine dealloc_mem_TPS
implicit none
if(allocated(sq)) deallocate(sq)
end subroutine dealloc_mem_TPS


subroutine parallelized_TPS(nvec2,itr,ene,restart)
use parameters, only : n, ldim, max_ldim, ibond, itr_counter
use working_area, only : wwk
use data_mem, only : ipair
use basic, only : datack
use parallelized_diag, only : parallelized_set_initial_vec_p
use parallel_lib
implicit none
logical, intent(in) :: restart
integer, intent(in) :: nvec2
integer, intent(inout) :: itr
double precision :: ene(:)
integer :: start, i
!!
if(nvec2<0)then
  print *,' #(E06)# Wrong value given to nvec in lnc1'
  print *,'         Only the eigenvalues are calculated'
  stop
endif
call datack(ipair,ibond,n)
!!
write(*,*) 'start TPS'
!!
start=itr_counter(9)
do i=start,itr_counter(10)
!*** initialization
  !$OMP PARALLEL
  !$OMP WORKSHARE
    wwk(1:max_ldim,1:4)=dcmplx(0d0,0d0)
  !$OMP END WORKSHARE
  !$OMP END PARALLEL
  if(.not.restart) call parallelized_set_initial_vec_p(wwk(:,1)) 
!!
  call parallelized_TPSz(i,nvec2,ene,itr,wwk(:,1),wwk(:,2),wwk(:,3),wwk(:,4))
  itr_counter(9)=itr_counter(9) + 1
enddo
end subroutine parallelized_TPS
!!
subroutine eval_phys_val(ict,alp,bet,l_max,ldim,vg,v1,v0,vt)
use working_area, only : coef
implicit none
integer, intent(in) :: ict, ldim
double precision, intent(in) :: alp, bet, l_max
complex(KIND(0d0)), intent(in) :: vg(:), v1(:), v0(:)
complex(KIND(0d0)), intent(inout) :: vt(:)
double precision :: res
integer :: offset, i
coef(ict,1)=dble(2*ict)/alp
coef(ict,2)=l_max
!    internal energy 
coef(ict,3)=bet                                  ! k-1 <i|h|i>
coef(ict,4)=l_max - alp                          ! normlization factor of |k>
!
call calc_internal_energy(ldim,l_max,res,vg,v1,v0)
coef(ict,5)=res

call calc_nn_correlation(ldim,vg,v1,vt)
coef(ict,6)=cz(1)
coef(ict,7)=cz(2)

offset=7
call calc_StaticStractureFactor(ldim,vg,v1,vt)
do i =1,sq_num
    coef(ict,offset+i)=sq(i)
enddo
!!
end subroutine eval_phys_val
!!
!!
subroutine calc_internal_energy(ldim,l_max,res,va,vb,vc)
use parallelized_diag
implicit none
integer, intent(in) :: ldim
double precision, intent(in) :: l_max
complex(KIND(0d0)), intent(in) :: va(:),vb(:),vc(:)
double precision, intent(inout) :: res
complex(KIND(0d0)) :: val1, val2
call calc_vec_inprdct(ldim,va,vc,val1)
call calc_vec_inprdct(ldim,va,vb,val2)
res=DREAL(l_max*val2-val1)
end subroutine calc_internal_energy
!!
!!
subroutine calc_StaticStractureFactor(ldim,vg,v1,vt)
use parallelized_diag
use data_mem, only : q_list, q_list_symmetric
use DynamicStractureFactor, only : intalv_zz, intalv_xy
implicit none
integer, intent(in) :: ldim
complex(KIND(0d0)), intent(in) :: v1(:), vg(:)
complex(KIND(0d0)), intent(inout) :: vt(:)
complex(KIND(0d0)) :: cnorm, val1,val2
double precision :: q(2)
integer :: i
do i=1,sq_num
  q(1)=q_list(1,q_list_symmetric(i)); q(2)=q_list(2,q_list_symmetric(i))
  call calc_sq_zz(vt,v1,q)
  call calc_vec_inprdct(ldim,v1,vt,val1)
  call calc_vec_inprdct(ldim,vg,vt,val2)
  sq(4*i-3)=DREAL(val1)
  sq(4*i-2)=DREAL(val2)
  call calc_sq_pp(vt,v1,q)
  call calc_vec_inprdct(ldim,v1,vt,val1)
  call calc_vec_inprdct(ldim,vg,vt,val2)
  sq(4*i-1)=DREAL(val1)
  sq(4*i  )=DREAL(val2)
enddo
contains
subroutine calc_sq_zz(v,x,q)
use parameters, only : n, ldim, max_ldim, b_bit, p_bit
use data_mem, only : site_vec, list1
use parallel_lib, only : lme, barrier
use parallelized_diag, only : calc_vec_inprdct
use openmp_lib
implicit none
double precision, intent(in) :: q(:)
complex(KIND(0d0)), intent(in) :: x(:)
complex(KIND(0d0)), intent(inout) :: v(:)
double precision :: qnorm, norm2
integer :: i, j, is, is0, istate, ibit, bit_shift, t_lme
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
!$OMP& PRIVATE(i,istate,is,is0) &
!$OMP& FIRSTPRIVATE(qnorm,t_lme,bit_shift,ibit,bp)
!$OMP do 
  do i=1, ldim
           istate=list1(i)
           is0=iand(istate,1)
           is=iand(ishft(istate,-bit_shift),1)
    if(bp) is=iand(ishft(t_lme, -bit_shift),1)
    v(i)=v(i)+(dble(is)-0.5d0)*(dble(is0)-0.5d0)*dcmplx(dcos(qnorm),dsin(qnorm))*x(i)
  enddo
!$OMP END DO
!$OMP END PARALLEL
enddo
end subroutine calc_sq_zz
!!
subroutine calc_sq_pp(v,x,q)
use parameters, only : n, ldim, max_ldim, b_bit, p_bit
use data_mem, only : site_vec, list1
use parallel_lib, only : lme, barrier
use parallelized_diag, only : calc_vec_inprdct
use openmp_lib
implicit none
double precision, intent(in) :: q(:)
complex(KIND(0d0)), intent(in) :: x(:)
complex(KIND(0d0)), intent(inout) :: v(:)
double precision :: qnorm, norm2
integer :: i, j, is, is0, istate, ibit, bit_shift, t_lme
logical :: bp
double precision :: my_xor(0:1,0:1)=0d0
my_xor(0,1)=0.25d0; my_xor(1,0)=0.25d0
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
!$OMP& PRIVATE(i,istate,is,is0) &
!$OMP& FIRSTPRIVATE(qnorm,t_lme,bit_shift,ibit,bp)
!$OMP do 
  do i=1, ldim
           istate=list1(i)
           is0=iand(istate,1)
           is=iand(ishft(istate,-bit_shift),1)
    if(bp) is=iand(ishft(t_lme, -bit_shift),1)
    v(i)=v(i)+my_xor(is,is0)*dcmplx(dcos(qnorm),dsin(qnorm))*x(i)
  enddo
!$OMP END DO
!$OMP END PARALLEL
enddo
end subroutine calc_sq_pp
!!
end subroutine calc_StaticStractureFactor
!!
!!
subroutine calc_nn_correlation(ldim,vg,v1,vt)
use parallelized_diag
implicit none
integer, intent(in) :: ldim
complex(KIND(0d0)), intent(in) :: v1(:), vg(:)
complex(KIND(0d0)), intent(inout) :: vt(:)
complex(KIND(0d0)) :: val1,val2
!!
call calc_zz_bond_corr(vt,v1)
call calc_vec_inprdct(ldim,v1,vt,val1)
call calc_vec_inprdct(ldim,vg,vt,val2)
cz(1)=DREAL(val1)
cz(2)=DREAL(val2)
!!
!!
contains 
subroutine calc_zz_bond_corr(v,x)
use parameters, only : n, max_ldim, ldim, b_bit, p_bit, ibond
use parallel_lib, only : lme
use data_mem, only : ipair, ibond_type, list1
implicit none
complex(KIND(0d0)), intent(in) :: x(:)
complex(KIND(0d0)), intent(inout) :: v(:)
integer :: i, j, t_lme
integer :: istate, bit_shift1, bit_shift2, bp1, bp2, is1, is2, is3, is4
t_lme = lme
!$OMP PARALLEL
!$OMP do
do i =1,max_ldim
  v(i)=dcmplx(0d0,0d0)
enddo
!$OMP END DO
!$OMP END PARALLEL 
!!
do i = 1, ibond
  if(ibond_type(i) /= 3) cycle
  is1=ipair(2*i-1)
  is2=ipair(2*i  )
  bit_shift1=is1
  bp1 = .false.
  if(is1>b_bit-1)then
    bit_shift1=is1-(b_bit)
    bp1 = .true.
  endif
  bit_shift2=is2
  bp2 = .false.
  if(is2>b_bit-1)then
    bit_shift2=is2-(b_bit)
    bp2 = .true.
  endif
!! zz bond correlation
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j,istate,is3,is4) &
!$OMP& FIRSTPRIVATE(t_lme,bit_shift1,bp1,bit_shift2,bp2)
!$OMP do 
  do j=1, ldim
            istate=list1(j)
            is3=iand(ishft(istate,-bit_shift1),1)
    if(bp1) is3=iand(ishft(t_lme, -bit_shift1),1)
            is4=iand(ishft(istate,-bit_shift2),1)
    if(bp2) is4=iand(ishft(t_lme, -bit_shift2),1)
    v(j)=v(j)+(dble(is3)-0.5d0)*(dble(is4)-0.5d0)*x(j)
  enddo
!$OMP END DO
!$OMP END PARALLEL
enddo
end subroutine calc_zz_bond_corr
end subroutine calc_nn_correlation
!!
!************ eigenvalues by the Lanczos method
subroutine parallelized_TPSz(ind,nvec2,ene,itr,v1,v0,vg,vt)
use parameters, only : ldim, max_ldim, max_itr, itr_counter,progress,max_comp_time, localHam_max,dsf_comp
use working_area, only : alpha, beta, coef
use parallelized_diag, only : p_mltply_v3, calc_beta, calc_alpha_beta2
use parallel_lib
use openmp_lib
use def_output, only : interval_dump_read,interval_dump
use timer
implicit none
integer, intent(in) :: nvec2, ind
double precision, intent(inout) :: ene(:)
complex(KIND(0d0)), intent(inout) ::vg(:),v1(:),v0(:), vt(:)
integer, intent(inout) :: itr
integer :: i, j, start
double precision :: alpha1, beta1, gamma1, ebefor, eps, xnorm
complex(KIND(0d0)) :: temp1, temp2
logical :: loop=.false.
integer :: t_max_itr
t_max_itr=itr_counter(2)
!! initial setting
!!
!*** calc_alpha_beta
!
! |k> = (L-H)|k-1> : v1 = H * v0 / ||v0||
! alpha = <k-1|(L-H)|k-1> : u_(k-1)
! beta = <k|k> : v1*v1
!!
if(itr_counter(1)<3)then
  call calc_alpha_beta2(ldim,alpha1,beta1,v1,v0)
  alpha(1)=localHam_max - alpha1 ! The zero-th order
  beta(1) =beta1  !
  if(lme==0) write(*,'(A,2e25.6)') 'ITR :     1',2d0/alpha1,localHam_max-alpha1
  call eval_phys_val(1,alpha1,beta1,localHam_max,ldim,vg,v1,v0,vt)
else
  call interval_dump_read(position=progress, data_size=max_ldim, &
                          vec_size=t_max_itr, sa=alpha, sb=beta ,sc=coef, se=ene,&
                          vec1=v0,vec2=v1, vec3=vg)
       alpha1=alpha(itr_counter(1)-1)
       beta1=beta(itr_counter(1)-1)
endif
!!
!*** iteration
start=itr_counter(1)
do i=start,itr_counter(2)
  if(tell_me_time()>max_comp_time)then
      call interval_dump(position=progress, data_size=max_ldim, &
                         vec_size=t_max_itr, sa=alpha, sb=beta, sc=coef, se=ene,&
                         vec1=v0, vec2=v1, vec3=vg)
  endif
!!
!!
  xnorm=1d0/beta1
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(ldim) &
!$OMP& PRIVATE(j) 
!$OMP do 
  do j=1,ldim
    vg(j)=v1(j)
    v1(j)=xnorm*v0(j)
    v0(j)=dcmplx(0d0,0d0)
  enddo
!$OMP END DO
!$OMP END PARALLEL
!!
    call calc_alpha_beta2(ldim,alpha1,beta1,v1,v0)
    alpha(i)=localHam_max - alpha1 ! The i-1 th order
    beta(i) =beta1  ! 
    call eval_phys_val(i,alpha1,beta1,localHam_max,ldim,vg,v1,v0,vt)
!!
    if(lme==0) write(*,'(A,I5,2e25.6)') 'ITR : ',i,dble(2*i)/alpha1,localHam_max-alpha1
!!
!  if(beta(i).lt.0.5d-30)then
!            print *,' #(E07)# Tridiagonalization unsuccessful in lnc1'
!            print *,'         Beta(i) is too small at i=  ',i
!            stop
!  endif
  if(mod(i,100)==1) call dump_alpha_beta_coef(dsf_comp,ind,i+1,v0,v1,vg)
  itr_counter(1)=itr_counter(1)+1
enddo
itr=max(max_itr-1,ldim)
call dump_alpha_beta_coef(dsf_comp,ind,itr+1,v0,v1,vg)
itr_counter(1)=2
!!
end subroutine parallelized_TPSz
!!
!!
subroutine calc_q(resdat)
use parameters, only : n, ldim, max_ldim, max_itr, itr_counter, localHam_max
use working_area, only : alpha, beta, coef
implicit none
integer :: i,j,k
double precision :: resdat(:,:)
double precision :: x0, x1, fm, n_inv_temp
integer, allocatable :: tmp_list(:)
double precision, allocatable :: f0(:)
allocate(f0(max_itr),tmp_list(max_itr))
!! beta
do i=1,itr_counter(2)
  coef(i,2)=dble(2*i) / (localHam_max - coef(i,3))
enddo
!!
do j=1,itr_counter(2)
  n_inv_temp=coef(j,2)
  do i=1,itr_counter(2)
    f0(i)=fac(i-1,n_inv_temp)
  enddo
  call lsort(f0,tmp_list)
  fm=maxval(f0(:))
  f0=f0-fm
  do i=1,itr_counter(2)
    f0(i)=dexp(f0(i))
  enddo
  x0=0d0
  x1=0d0
  do i=itr_counter(2),1,-1
    k=tmp_list(i)
    x0=x0+f0(k)*(1d0      +n_inv_temp/dble(2*(k-1)+1)* (localHam_max-alpha(k))/beta(k) )
    x1=x1+f0(k)*(coef(k,3)+n_inv_temp/dble(2*(k-1)+1)* coef(k,4)                       )
  enddo
  resdat(j,1)=x1/x0
  write(*,*) j,x1/x0
enddo

contains
function fac(k,Nbeta)
implicit none
integer :: k
double precision :: Nbeta, fac
integer :: i, j, m
double precision :: z1, z2
!!if(Nbeta>1d0)then
!!  z1=dble( 2*(k-max_itr) )*( dlog(Nbeta) )
!!else
  z1=dble( 2*k )*( dlog(Nbeta) )
!!endif
z2=0d0
do i=1,2*k+1
  z2=z2+dlog(dble(i))
enddo
if (k==0.and.z1>z2) write(*,*) 'warning'
fac=z1-z2
end function fac
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
end subroutine calc_q
!
!
!!

subroutine dump_alpha_beta_coef(comp,indx1,indx2,v1,v2,v3)
use parallel_lib, only : get_me
use working_area, only : alpha, beta, coef
implicit none
integer, intent(in) :: indx1,indx2
integer :: num_data1, num_data2
character(LEN=2), intent(in) :: comp
character(LEN=128) :: FILENAME, IDENTIFY
complex(KIND(0d0)), intent(in) :: v1(:),v2(:),v3(:) 
integer :: j, pme
pme=get_me()
write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'SalphaSbeta_me',pme,'seed',indx1,trim(comp)
FILENAME=trim(IDENTIFY)
j=100+pme+1
num_data1=ubound(alpha,1)
num_data2=ubound(coef,2)
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
write(j) num_data1
write(j) num_data2
write(j) alpha
write(j) beta
write(j) coef
close(j)
!write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'vec_me',pme,'seed',indx1,trim(comp)
!FILENAME=trim(IDENTIFY)
!open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
!write(j) indx2
!write(j) v1
!write(j) v2
!write(j) v3
!close(j)
end subroutine dump_alpha_beta_coef

subroutine read_alpha_beta_coef(comp,indx1,indx2,v1,v2,v3)
use parallel_lib, only : get_me
use working_area, only : alpha, beta, coef
implicit none
integer, intent(inout) ::  indx1, indx2
complex(KIND(0d0)), intent(inout) :: v1(:), v2(:), v3(:)
character(LEN=2), intent(in) :: comp
character(LEN=128) :: FILENAME, IDENTIFY
integer :: j, k, m, pme, num_data1, num_data2

double precision, allocatable :: ta(:), tb(:), tc(:,:)

pme=get_me()
j=100+pme+1

write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'SalphaSbeta_me',pme,'seed',indx1,trim(comp)
FILENAME=trim(IDENTIFY)
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
read(j) num_data1
read(j) num_data2
allocate(ta(num_data1),tb(num_data1),tc(num_data1,num_data2))
read(j) ta
read(j) tb
read(j) tc
close(j)
alpha(:)=0d0; beta(:)=0d0; coef(:,:)=0d0
do k =1 ,num_data1+1
  alpha(k)=ta(k); beta(k)=tb(k)
enddo
do m= 1 ,num_data2
do k =1 ,num_data1
  coef(k,m)=tc(k,m)
enddo
enddo
deallocate(ta,tb,tc)
write(IDENTIFY,'(a,i5.5,a,i5.5,a)') 'vec_me',pme,'seed',indx1,trim(comp)
FILENAME=trim(IDENTIFY)
open(j,file=trim(FILENAME),form='unformatted',access='sequential',status='unknown')
read(j) indx2
read(j) v1
read(j) v2
read(j) v3
close(j)
end subroutine read_alpha_beta_coef


!!
end module parallelized_TPS_method
!**********************************************************************************************************************
