module parallelized_diag
implicit none
private
logical, private, save :: ch_mk_list3=.true.
integer, private, save :: ihf0, ihf, ihfbit, ilft, irght
public p_mltply_v3, calc_vec_inprdct, calc_beta, calc_alpha_beta, parallelized_set_initial_vec_p, &
       mk_list12, mk_list3, alloc_listt, get_listt, dealloc_listt, get_mask, set_ch_mk_list3, calc_alpha_beta2
contains
!*************** initial vector *****************
subroutine set_ch_mk_list3(l)
implicit none
logical, intent(in) :: l
ch_mk_list3=l
end subroutine set_ch_mk_list3

subroutine parallelized_set_initial_vec_p(v)
use parameters, only : n, ldim, p_bit, b_bit, q_spin, total_sz, pi
use parallel_lib, only : allreduce_rdata, lme, barrier
use data_mem, only : list1
use rand_num1
use openmp_lib
implicit none
complex(KIND(0d0)), intent(inout) :: v(:)
integer :: i, s, p, j, m
integer :: tot_sz_zero
double precision :: dnorm, dnorm2
double precision :: res, r, x, y(2)
v(1:ldim)=dcmplx(0d0,0d0)
dnorm=0d0
m=0
tot_sz_zero = (q_spin-1)*n/2
!!!!
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,p,s,j)
!$OMP do reduction(+:dnorm)
do i=1,ldim
  p=0
  s=list1(i)
  do j=0,b_bit-1
    p=p+mod(s,q_spin)
    s=s/q_spin
  enddo
  s=lme
  do j=0,p_bit-1
    p=p+mod(s,q_spin)
    s=s/q_spin
  enddo
!  v(i)=0d0
!  if(p==tot_sz_zero.and.m==0)then
!  v(i)=dcmplx(1d0,0d0)
  call random_number(y)
  y(2)=2d0*pi*y(1)
  v(i)=dcmplx(dcos(y(2)),dsin(y(2)))
!  if((p==total_sz).and.((ldim/3<i).and.(i<ldim/2)))then
!  if((p==tot_sz_zero).and.lme==0.and.m==0)then
!  if(i==13.and.lme==0.and.m==0)then
!!  if((iand(p,1)==0).and.(ldim/3<i<ldim/2))then
!  if((iand(p,1)==0).and.lme==0.and.m==0)then
!  if(p==total_sz.and.m==0.and.lme==0)then
!    v(i)=dcmplx(1d0,0.d0)
!    m=m+1
!  endif
  dnorm=dnorm+dconjg(v(i))*v(i)
enddo
!$OMP END DO
!$OMP END PARALLEL
!!
!!!!!!!!
call allreduce_rdata(1,dnorm,res)
dnorm2=0.d0
r=1d0/dsqrt(res)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(r) &
!$OMP& PRIVATE(i)
!$OMP do reduction(+:dnorm2)
do i=1,ldim
  v(i)=v(i)*r
  dnorm2=dnorm2+dconjg(v(i))*v(i)
enddo
!$OMP END DO
!$OMP END PARALLEL
call allreduce_rdata(1,dnorm2,res)
r=1d0/dsqrt(res)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(r) &
!$OMP& PRIVATE(i)
!$OMP do
do i=1,ldim
  v(i)=v(i)*r
enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine parallelized_set_initial_vec_p

subroutine calc_vec_inprdct(data_num,x_in1,x_in2,norm)
!use parameters, only : ldim
use parallel_lib, only : allreduce_data
use openmp_lib
implicit none
integer, intent(in) :: data_num
complex(KIND(0d0)), intent(in) :: x_in1(:), x_in2(:)
complex(KIND(0d0)), intent(inout) :: norm
integer :: j
complex(KIND(0d0)) :: res
norm=dcmplx(0d0,0d0)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j)
!$OMP do reduction(+:norm)
do j=1,data_num
  norm=norm+dconjg(x_in1(j))*x_in2(j)
enddo
!$OMP END DO
!$OMP END PARALLEL
call allreduce_data(1,norm,res)
norm=res
end subroutine calc_vec_inprdct
!
subroutine calc_alpha_beta2(ldim,t_alp,t_bet,v1,v0)
use parallel_lib, only : barrier
implicit none
integer, intent(in) :: ldim
double precision, intent(inout) :: t_alp, t_bet
complex(KIND(0d0)), intent(inout) :: v1(:)
complex(KIND(0d0)), intent(inout) :: v0(:)
complex(KIND(0d0)) :: t_alp1
complex(KIND(0d0)) :: t_bet1
call barrier
call p_mltply_v3(v1,v0,t_alp1,0d0)
call calc_vec_inprdct(ldim,v0,v0,t_bet1)
t_alp=DREAL(t_alp1) ! p_mltply_v3 : v0 = v0 + H * v1
t_bet=dsqrt(DREAL(t_bet1))
end subroutine calc_alpha_beta2
!
subroutine calc_alpha_beta(ldim,t_alp,t_bet,v1,v0)
use parallel_lib, only : barrier
implicit none
integer, intent(in) :: ldim
double precision, intent(inout) :: t_alp, t_bet
complex(KIND(0d0)), intent(inout) :: v1(:)
complex(KIND(0d0)), intent(inout) :: v0(:)
complex(KIND(0d0)) :: t_alp1
double precision :: t_bet1
integer :: i
call barrier
call p_mltply_v3(v1,v0,t_alp1,0d0)
call calc_beta(ldim,t_alp1,t_bet1,v1,v0)
t_alp=DREAL(t_alp1)
t_bet=dsqrt(t_bet1)
end subroutine calc_alpha_beta

subroutine calc_beta(data_num,alpha,beta,v1,v0)
use parallel_lib, only : allreduce_rdata
use openmp_lib
implicit none
integer, intent(in) :: data_num
complex(KIND(0d0)), intent(in) :: v1(:), v0(:)
complex(KIND(0d0)), intent(in) :: alpha
double precision, intent(inout) :: beta
double precision :: res
integer :: i
complex(KIND(0d0)) :: K_alpha
integer :: K_data_num
beta=0d0
K_alpha=alpha
K_data_num=data_num
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(K_alpha,K_data_num) &
!$OMP& PRIVATE(i)
!$OMP do reduction(+:beta)
do i=1,K_data_num
  beta=beta+dconjg(v0(i)-K_alpha*v1(i))*(v0(i)-K_alpha*v1(i))
enddo
!$OMP END DO
!$OMP END PARALLEL
call allreduce_rdata(1,beta,res)
beta=res
end subroutine calc_beta

subroutine p_mltply_v3(v1,v0,prdct,eperbd)
use parameters, only : ldim, total_sz, max_ldim
use working_area, only : rr1
use parallel_lib
use timer
use openmp_lib
implicit none
complex(KIND(0d0)), intent(inout) :: v1(:)
complex(KIND(0d0)), intent(inout) :: v0(:)
complex(KIND(0d0)), intent(inout) :: prdct
double precision, intent(in) :: eperbd
integer :: k, real_target_cpu, real_pre_target
integer :: target_cpu(2), target_sz, crr_ind, lst_ind, local_loop_indx
integer :: tdim_l(0:num_tot_ex_cpus)
complex(KIND(0d0)) :: res, lprdct, prdct2
logical :: send_logical, recv_logical
integer :: itag1,itag2,i
integer :: ii
!!
prdct=dcmplx(0d0,0d0)
target_sz=total_sz
real_pre_target=lme*lnp+lme
local_loop_indx=1
if(ch_mk_list3) call mk_list3(lme,target_sz,tdim_l)
call timer_start(5)
tdim_l(0)=ldim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
crr_ind=0
lst_ind=0
do k=ex_loop_start, num_int_bond
!  call barrier
  if(.not.exlog2(k))cycle
  real_target_cpu=real_cpu_list(k)
  target_cpu(1)=cpu_list_send(k)
  target_cpu(2)=cpu_list_recv(k)
!!
  if(real_target_cpu.ne.real_pre_target)then
    real_pre_target=real_target_cpu
    crr_ind=crr_ind+1
    lst_ind=lst_ind+1
!$OMP PARALLEL
!$OMP WORKSHARE
    rr1(1:max_ldim)=dcmplx(0d0,0d0)
!$OMP END WORKSHARE
!$OMP END PARALLEL
    itag1=lnp*k+target_cpu(1)
    itag2=lnp*k+lme
    call mpi_isendrecv(target_cpu,max_ldim,v1,rr1,itag1,itag2)
    do
      call mpi_test_isendrecv(send_logical,recv_logical)
      if(send_logical.and.recv_logical) exit
!!
      if(local_loop_indx<=max_local_loop) then
        call local_prdct(local_loop_indx,0,v1,v1,v0,lprdct,eperbd)
        prdct=prdct+lprdct
        local_loop_indx=local_loop_indx+1
      endif
    enddo
    call mpi_mywait
  endif
!!
  call local_prdct(k,label_list(k),v1,rr1,v0,lprdct,eperbd)
  prdct=prdct+lprdct
enddo
!!
do while(local_loop_indx<=max_local_loop)
    call local_prdct(local_loop_indx,0,v1,v1,v0,lprdct,eperbd)
    prdct=prdct+lprdct
    local_loop_indx=local_loop_indx+1
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call allreduce_data(1,prdct,res)
prdct=res
call timer_stop(5)

end subroutine p_mltply_v3

subroutine local_prdct(p,ind,v1,vt,v0,lprdct,eperbd)
use parameters, only : ldim, max_ldim, q_spin
use data_mem, only : list1, list3
use parallel_lib
use timer
implicit none
integer, intent(in) :: p, ind
complex(KIND(0d0)), intent(in) :: v1(:)
complex(KIND(0d0)), intent(in) :: vt(:)
complex(KIND(0d0)), intent(inout) :: v0(:)
complex(KIND(0d0)), intent(inout) :: lprdct
double precision, intent(in) :: eperbd
integer :: isite1, isite2, is1, is2, ibd
integer :: ia, ib, ie, l1, l2
integer :: interaction, elemnt
double precision :: wt
!!
call timer_start(9)
  interaction=int_type(p)
  ibd=ibond_type_p2(p)
  wt = bondwt_p2(p)
  elemnt=matrix_element2(p)
  ia = 1-iand(ishft(interaction-1,-1),1)
  ib = 1-iand(ishft(interaction  ,-1),1)
  isite1=site_pair(1,p)
  is1=q_spin**isite1
  isite2=site_pair(2,p)
  is2=q_spin**isite2
  call H_dot_v(lme,elemnt,ia,ib,ldim,isite1,isite2, &
               is1,is2,ibd,wt,list1,list3(:,:,ind),v0,v1,vt,lprdct,eperbd)
call timer_stop(9)
!!
end subroutine local_prdct
!
!
subroutine H_dot_v(Iam,ml,ia,ib,dim,isite1,isite2, &
                   is1,is2,ibd,wt,list_s,list_lib,v0,v1,vt,lprdct,eperbd)
use parameters, only : max_num_list2, ldim, q_spin, size_of_local_matrix
use data_mem, only : xlocal_Ham, llocal_Ham
use openmp_lib
use parallel_lib
implicit none
integer, intent(in) :: Iam, dim, ml, ia, ib
integer, intent(in) :: isite1, isite2, is1, is2, ibd
integer, intent(in) :: list_s(:),list_lib(2,0:max_num_list2)
complex(KIND(0d0)), intent(inout) :: v0(:)
complex(KIND(0d0)), intent(in) :: v1(:)
complex(KIND(0d0)), intent(in) :: vt(:)
complex(KIND(0d0)), intent(out) :: lprdct
double precision, intent(in) :: wt,eperbd
double precision :: base
complex(KIND(0d0)) :: factor2
integer :: iex, l0, it0, it1, it2, it, istat, ibit1, ibit2, Iamia, Iamib
integer :: iais1, ibis2, ic, l1, l2
!! For K
integer :: Kia,Kib,Kibd,Kisite1,Kisite2,Kis1,Kis2,Kml0,Kml1,Kml2
double precision :: Kwt,Keperbd
Kisite1=isite1; Kisite2=isite2;
Kis1=is1; Kis2=is2; Kibd=ibd;
Kia=ia;Kib=ib
Keperbd=-eperbd; 
Kml0 = ml
Kml1 = mod(ml,q_spin)
Kml2 = ml/q_spin
!! End For K
iais1=ia*is1; ibis2=ib*is2;
Iamia=Iam*(1-ia); Iamib=Iam*(1-ib);
!!
!!
lprdct=dcmplx(0d0,0d0)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(Kibd,Kwt,Kisite1,Kisite2,Kml0,Kml1,Kml2,&
!$OMP& Iamia,Iamib,Kia,Kib,Kis1,Kis2,Keperbd,irght,ilft,q_spin,iais1,ibis2)&
!$OMP& PRIVATE(l0,ibit1,ibit2,istat,it0,it1,it2,it,factor2)
!$OMP do reduction(+:lprdct)
do l0=1,dim
    it0=list_s(l0)
    ibit1= mod((it0*Kia + Iamia)/Kis1, q_spin)
    ibit2= mod((it0*Kib + Iamib)/Kis2, q_spin)
    istat= ibit2*q_spin+ibit1
    if(llocal_Ham(Kml0,istat,Kibd)) cycle
    it0=it0+iais1*(Kml1-ibit1)+ibis2*(Kml2-ibit2)
    it1=mod(it0,irght)
    it2=it0/ilft
    it=list_lib(1,it1)+list_lib(2,it2)
    factor2 = xlocal_Ham(Kml0,istat,Kibd)
    if(Kml0==istat) factor2 = factor2 + Keperbd
    v0(l0)=v0(l0)+factor2               *vt(it)
    lprdct=lprdct+factor2*dconjg(v1(l0))*vt(it)
enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine H_dot_v
!
subroutine get_mask(right,left,half)
implicit none
integer, intent(out) :: right, left, half
right=irght
left=ilft
half=ihf
end subroutine get_mask
!!
subroutine mk_list12(tot_sz,Iam)
use parameters, only : n, ldim, b_bit, p_bit, max_num_list2, num_list3, parity_type, max_ldim, dsf_comp, q_spin, ipm
use parallel_lib, only : mpi_get_maxint, num_tot_ex_cpus
use data_mem, only : list1,list3,sz_list_p,sz_list_b
use openmp_lib
implicit none
integer, intent(in) :: tot_sz
integer, intent(in) :: Iam
integer :: icnt, ja, jb, ibpatn, ic
integer :: i, j, isz, ia, ib, is, tmp, psz, pstat, bsz, bstat, i1, i2
integer :: max_a, max_b
!!!! initialization
      ihf0=b_bit/2
      ihf=(b_bit+1)/2
      ihfbit=q_spin**ihf
      irght=ihfbit
      ilft=q_spin**ihf
!!
      allocate(sz_list_p(0:q_spin**p_bit-1),sz_list_b(0:q_spin**ihf-1))
      do j=0,q_spin**p_bit-1
        pstat=j
        psz=0
        do i=0,p_bit-1
          psz=psz+mod(pstat,q_spin)
          pstat=pstat/q_spin
        enddo 
        sz_list_p(j)=psz
      enddo
      do j=0,q_spin**ihf-1
        bstat=j
        bsz=0
        do i=0,ihf-1
          bsz=bsz+mod(bstat,q_spin)
          bstat=bstat/q_spin
        enddo 
        sz_list_b(j)=bsz
      enddo
!!
      is=ipm
!!      if(dsf_comp=='xx') is=1
!!
      icnt=0; ja=0; jb=0; ibpatn=0
      ic = 1; max_a=0; max_b=0
!!
!!!!
!! main loop
psz=sz_list_p(Iam)
do i2=0,q_spin**ihf0 -1
do i1=0,q_spin**ihf  -1
  isz=sz_list_b(i1)+sz_list_b(i2)
  if(parity_type==0)then
    if((psz+isz).ne.(is+tot_sz)) cycle
  endif
  icnt=icnt+1
  ia=i1
  ib=i2
  max_a=max(max_a,ia)
  max_b=max(max_b,ib)
  if(ib==ibpatn)then
    ja=ja+1
  else
    ibpatn=ib
    ja=1
    jb=icnt-1
  endif
enddo
enddo
ldim=icnt
tmp=max(max_a,max_b)
call mpi_get_maxint(tmp,max_num_list2)
tmp=icnt
call mpi_get_maxint(tmp,max_ldim)
!!write(*,*) 'ihf0, ihf',ihf0,ihf
write(*,*) 'Total Sz',Iam,psz,tot_sz
write(*,*) 'Size of local vectors,Maximum Size',icnt,max_ldim
write(*,*) 'ME, Maximum size of library',Iam,max_b,max_a
write(*,*) 'NUM tot EX CPUS',num_tot_ex_cpus
allocate(list1(ldim))
!!allocate(list3(2,0:max_num_list2,0:num_list3))
allocate(list3(2,0:max_num_list2,0:num_tot_ex_cpus))
list3(:,:,:)=-1
!list1(:)=-1
list1(:)=0
!!
if(ldim.lt.3)then
  print *,' #(E02A)# Incorrect ldim or n given to sz'
  stop
endif
!!
      icnt=0; ja=0; jb=0; ibpatn=0
      ic = 1
!!!! main loop
do i2=0,q_spin**ihf0 -1
do i1=0,q_spin**ihf  -1
  isz=sz_list_b(i1)+sz_list_b(i2)
  if(parity_type==0)then
    if((psz+isz).ne.(tot_sz+is)) cycle
  endif
  icnt=icnt+1
  if(icnt>ldim)then
    print *,' #(E02B)# Incorrect ldim or n given to sz'
    stop
  endif
  ia=i1
  ib=i2
  if(ib==ibpatn)then
    ja=ja+1
  else
    ibpatn=ib
    ja=1
    jb=icnt-1
  endif
  list1(icnt)=ib*(q_spin**ihf)+ia
  list3(1,ia,0)=ja
  list3(2,ib,0)=jb
enddo
enddo
if(icnt==ldim)  return
print *,' #(E02C)# Incorrect ldim or n given to sz'
stop
end subroutine mk_list12
!!
!!
subroutine mk_list3(Iam,target_sz,tdim_l)
use parameters, only : n, max_num_list2, parity_type, dsf_comp, q_spin, ipm
use data_mem, only : list3, sz_list_p, sz_list_b
use parallel_lib, only : tot_ex_cpus, num_tot_ex_cpus
use openmp_lib
implicit none
integer, intent(in) :: Iam,target_sz
integer, intent(out) :: tdim_l(0:num_tot_ex_cpus)
!!!!
integer :: target_cpu
integer :: icnt, ja, jb, ibpatn
integer :: isz, ia, ib, psz, j0, i1, i2, is
integer :: c_max
!!!!
is=ipm
!!; if(dsf_comp=='xx') is=1
!!!!
c_max = num_tot_ex_cpus
!! Does not work 2015 Aug 6
!!!!$OMP PARALLEL &
!!!!$OMP DEFAULT(PRIVATE) &
!!!!$OMP SHARED(Iam,c_max,is,target_sz,list3,sz_list_p,sz_list_b, &
!!!!$OMP&       max_num_list2,parity_type,tdim_l,tot_ex_cpus,ihf,ihf0)
!!!!$OMP do SCHEDULE(STATIC,1)
!!
!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j0,icnt,ja,jb,ibpatn,target_cpu,psz,isz,i1,i2,ia,ib)
!$OMP do SCHEDULE(STATIC,1)
do j0=1, c_max
!!!! initialization
  icnt=0; ja=0 ;jb=0 ;ibpatn=0
  target_cpu = tot_ex_cpus(j0)
!!
  list3(1:2,0:max_num_list2,j0)=-1
!!
  psz=sz_list_p(target_cpu)
!!!! main loop
!!
do i2=0,q_spin**ihf0 -1
do i1=0,q_spin**ihf  -1
  isz=sz_list_b(i1)+sz_list_b(i2)
  if(parity_type==0)then
    if((psz+isz).ne.(target_sz+is)) cycle
  endif
  icnt=icnt+1
  ia=i1
  ib=i2
  if(ib==ibpatn)then
    ja=ja+1
  else
    ibpatn=ib
    ja=1
    jb=icnt-1
  endif
  list3(1,ia,j0)=ja
  list3(2,ib,j0)=jb
enddo
enddo
tdim_l(j0)=icnt
enddo
!$OMP END DO
!$OMP END PARALLEL
call set_ch_mk_list3(.false.)
end subroutine mk_list3
!!
!!
subroutine alloc_listt(target_sz,Iam,tdim)
use parameters, only : max_tdim, parity_type, q_spin
use data_mem, only : sz_list_p, sz_list_b
use parallel_lib, only : mpi_get_maxint
use openmp_lib
implicit none
integer, intent(in) :: target_sz
integer, intent(in) :: Iam
integer, intent(inout) :: tdim
integer :: icnt, isz, psz, i1, i2, tmp
!!
icnt=0
psz=sz_list_p(Iam)
do i2=0,q_spin**ihf0 -1
do i1=0,q_spin**ihf  -1
  isz=sz_list_b(i1)+sz_list_b(i2)
  if(parity_type==0)then
    if((psz+isz).ne.target_sz) cycle
  endif
  icnt=icnt+1
enddo
enddo
tmp=icnt
tdim=icnt
call mpi_get_maxint(tmp,max_tdim)
end subroutine alloc_listt
!!
subroutine get_listt(target_cpu,target_sz,Iam,pdim)
use parameters, only : parity_type, max_tdim, q_spin
use data_mem, only : listt, sz_list_p, sz_list_b
use openmp_lib
implicit none
integer, intent(in) :: target_sz, target_cpu
integer, intent(in) :: Iam
integer, intent(inout) :: pdim
integer :: icnt, isz, psz, i1, i2
if(.not.allocated(listt)) allocate(listt(max_tdim))
listt(:)=-1
!!!! main loop
psz=sz_list_p(target_cpu)
icnt=0
do i2=0,q_spin**ihf0-1
do i1=0,q_spin**ihf -1
  isz=sz_list_b(i1)+sz_list_b(i2)
  if(parity_type==0)then
    if((psz+isz).ne.target_sz) cycle
  endif
  icnt=icnt+1
  listt(icnt)=i2*q_spin**ihf+i1
enddo
enddo
pdim=icnt
end subroutine get_listt
!!
subroutine dealloc_listt
use data_mem, only : listt
implicit none
if(allocated(listt))deallocate(listt)
end subroutine dealloc_listt
!!
end module parallelized_diag

