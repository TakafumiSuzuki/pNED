module parallelized_diag
implicit none
private
logical, private, save :: ch_mk_list3=.true.
integer, private, save :: ihf0, ihf, ihfbit, ilft, irght
public p_mltply_v3, calc_vec_inprdct, calc_beta, calc_alpha_beta, parallelized_set_initial_vec_p, &
       mk_list12, mk_list3, alloc_listt, get_listt, dealloc_listt, get_mask, set_ch_mk_list3, &
       calc_alpha_beta2
contains
!*************** initial vector *****************
subroutine set_ch_mk_list3(l)
implicit none
logical, intent(in) :: l
ch_mk_list3=l
end subroutine set_ch_mk_list3

subroutine parallelized_set_initial_vec_p(v)
use parameters, only : n, ldim, p_bit, b_bit, pi
use parallel_lib, only : allreduce_rdata, lme, barrier
use rand_num1
use data_mem, only : list1
use openmp_lib
implicit none
complex(KIND(0d0)), intent(inout) :: v(:)
integer :: i, s, p, j, m
!integer, parameter :: seed=12348789
double precision :: dnorm, dnorm2
double precision :: res, r, x, y(2)
v(1:ldim)=dcmplx(0d0,0d0)
dnorm=0d0
m=0
!!!!
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,p,s,j,y)
!$OMP do reduction(+:dnorm)
do i=1,ldim
  p=0
  s=list1(i)
  do j=0,b_bit-1
    p=p+iand(s,1)
    s=ishft(s,-1)
  enddo
  s=lme
  do j=0,p_bit-1
    p=p+iand(s,1)
    s=ishft(s,-1)
  enddo
!!  v(i)=0d0
  call random_number(y)
  y(2)=2d0*pi*y(1)
  v(i)=dcmplx(dcos(y(2)),dsin(y(2)))
!!  if((p==n/2).and.((ldim/3<i).and.(i<ldim/2)))then
!!  if((p==n/2).and.lme==0.and.m==0)then
!!  if(p==n)then
!!  if(i==13.and.lme==0.and.m==0)then
!!  if((iand(p,1)==0).and.(ldim/3<i<ldim/2))then
!!  if((p==n/2).and.lme==0.and.m==0)then
!!  if(p==n/2)then
!!    v(i)=dcmplx(1d0,0d0)
!!    m=m+1
!!  endif
  dnorm=dnorm+dconjg(v(i))*v(i)
enddo
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!
!if(lme==0)then
!  call random_number(y)
!  i = min(int(ldim*y(1))+1,ldim)
!  call random_number(y)
!  v(i)=dcmplx(y(1),y(2))
!  dnorm=DREAL(dconjg(v(i))*v(i))
!endif
!!!!!!!!!
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


subroutine calc_alpha_beta(ldim,t_alp,t_bet,v1,v0)
use parallel_lib, only : barrier
implicit none
integer, intent(in) :: ldim
double precision, intent(inout) :: t_alp, t_bet
complex(KIND(0d0)), intent(inout) :: v1(:)
complex(KIND(0d0)), intent(inout) :: v0(:)
complex(KIND(0d0)) :: t_alp1
double precision :: t_bet1
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
use parameters, only : ldim, total_sz, max_ldim, num_list3
use working_area, only : rr1
use parallel_lib
use timer
use openmp_lib
implicit none
complex(KIND(0d0)), intent(inout) :: v1(:)
complex(KIND(0d0)), intent(inout) :: v0(:)
complex(KIND(0d0)), intent(inout) :: prdct
double precision, intent(in) :: eperbd
integer :: k
integer :: target_cpu, target_sz, pre_target, crr_ind, lst_ind, local_loop_indx
integer :: tdim_l(0:num_tot_ex_cpus)
complex(KIND(0d0)) :: res, lprdct, prdct2
logical :: send_logical, recv_logical
integer :: itag1,itag2,i
!!
prdct=dcmplx(0d0,0d0)
target_sz=total_sz
pre_target=lme
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
!!  call barrier
  if(.not.exlog2(k))cycle
  target_cpu=cpu_list(k)
!!
  if(target_cpu.ne.pre_target)then
    pre_target=target_cpu
    crr_ind=crr_ind+1
    lst_ind=lst_ind+1
!$OMP PARALLEL
!$OMP WORKSHARE
    rr1(1:max_ldim)=dcmplx(0d0,0d0)
!$OMP END WORKSHARE
!$OMP END PARALLEL
    itag1=lnp*k+target_cpu
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
  call local_prdct(k,lst_ind,v1,rr1,v0,lprdct,eperbd)
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
use parameters, only : ldim, max_ldim
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
integer :: isgn(0:3)=(/1,1,-1,-1/)
!!
call timer_start(9)
  interaction=int_type(p)
  elemnt=matrix_element2(p)
  ia = 1-iand(ishft(interaction-1,-1),1)
  ib = 1-iand(ishft(interaction  ,-1),1)
  ie = 2*ieor(ia,ib)*isgn(elemnt)
  l1=iand(elemnt,1)
  l2=iand(ishft(elemnt,-1),1)
  isite1=site_pair(1,p)
  is1=ishft(1,isite1)
  isite2=site_pair(2,p)
  is2=ishft(1,isite2)
  ibd=ibond_type_p2(p)
  wt = bondwt_p2(p)
!!
  call H_dot_v(lme,elemnt,ia,ib,ie,l1,l2,ldim,isite1,isite2, &
               is1,is2,ibd,wt,list1,list3(:,:,ind),v0,v1,vt,lprdct,eperbd)
call timer_stop(9)
!!
end subroutine local_prdct
!
!
subroutine H_dot_v(Iam,ml,ia,ib,ie,l1,l2,dim,isite1,isite2, &
                   is1,is2,ibd,wt,list_s,list_lib,v0,v1,vt,lprdct,eperbd)
use parameters, only : max_num_list2, ldim
use data_mem, only : xlocal_Ham, llocal_Ham
use openmp_lib
implicit none
integer, intent(in) :: Iam, dim, ml, ia, ib, ie
integer, intent(in) :: isite1, isite2, is1, is2, ibd, l1, l2
integer, intent(in) :: list_s(:),list_lib(2,0:max_num_list2)
complex(KIND(0d0)), intent(inout) :: v0(:)
complex(KIND(0d0)), intent(in) :: v1(:)
complex(KIND(0d0)), intent(in) :: vt(:)
complex(KIND(0d0)), intent(out) :: lprdct
double precision, intent(in) :: wt,eperbd
double precision :: base
complex(KIND(0d0)) :: factor2, hb(0:3,0:3)
integer :: iex, l0, it0, it1, it2, it, istat, ibit1, ibit2, Iamia, Iamib
integer :: iais1, ibis2, ic
!! For K
integer :: Kml,Kia,Kib,Kie,Kibd,Kl1,Kl2,Kisite1,Kisite2,Kis1,Kis2
double precision :: Kwt,Keperbd
Kml=ml; Kie=ie; Kl1=l1; Kl2=l2;
Kisite1=isite1; Kisite2=isite2;
Kis1=is1; Kis2=is2; Kibd=ibd;
Keperbd=-eperbd; 
!! End For K
iais1=ia*is1; ibis2=ib*is2;
Iamia=Iam*(1-ia); Iamib=Iam*(1-ib);
!!
lprdct=dcmplx(0d0,0d0)
hb(0:3,0:3)=xlocal_Ham(0:3,0:3,ibd)
ic=2*ib+ia
select case(ic)
case(0)
    ibit1=ishft(iand(Iamia, Kis1), -Kisite1)
    ibit2=ishft(iand(Iamib, Kis2), -Kisite2)
    istat=2*ibit2+ibit1
    iex = Kml+ Kie*ibit2
    if(llocal_Ham(iex,istat,Kibd)) return
    factor2 = hb(iex,istat)
    if(iex==istat) factor2 = factor2 + Keperbd
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(iais1,ibis2,Kl1,Kl2,ibit1,ibit2,factor2,irght,ilft,ihf,hb)&
!$OMP& PRIVATE(it0,it1,it2,it)
!$OMP do reduction(+:lprdct)
do l0=1,dim
    it0=list_s(l0)
    it0=it0+iais1*(Kl1-ibit1)
    it0=it0+ibis2*(Kl2-ibit2)
    it1=iand(it0,irght)
    it2=ishft(iand(it0,ilft),-ihf)
    it=list_lib(1,it1)+list_lib(2,it2)
    v0(l0)=v0(l0)+factor2               *vt(it)
    lprdct=lprdct+factor2*dconjg(v1(l0))*vt(it)
enddo
!$OMP END DO
!$OMP END PARALLEL
case(1)
    ibit2=ishft(iand(Iamib,      Kis2), -Kisite2)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(Kml,Kie,Kibd,Kwt,Kl1,Kl2,Kisite1,Kisite2,&
!$OMP& Kis1,Keperbd,irght,ilft,ihf,iais1,ibis2,Iamia,ibit2,hb)&
!$OMP& PRIVATE(l0,ibit1,istat,iex,it0,it1,it2,it,factor2)
!$OMP do reduction(+:lprdct)
do l0=1,dim
    ibit1=ishft(iand(list_s(l0), Kis1), -Kisite1)
    istat=2*ibit2+ibit1
    iex = Kml+ Kie*ibit2
    if(llocal_Ham(iex,istat,Kibd)) cycle
    it0=list_s(l0)
    it0=it0+iais1*(Kl1-ibit1)
    it0=it0+ibis2*(Kl2-ibit2)
    it1=iand(it0,irght)
    it2=ishft(iand(it0,ilft),-ihf)
    it=list_lib(1,it1)+list_lib(2,it2)
    factor2 = hb(iex,istat)
    if(iex==istat) factor2 = factor2 + Keperbd
    v0(l0)=v0(l0)+factor2               *vt(it)
    lprdct=lprdct+factor2*dconjg(v1(l0))*vt(it)
enddo
!$OMP END DO
!$OMP END PARALLEL
case(2)
    ibit1=ishft(iand(Iamia,      Kis1), -Kisite1)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(Kml,Kie,Kibd,Kwt,Kl1,Kl2,Kisite1,Kisite2,&
!$OMP& Kis2,Keperbd,irght,ilft,ihf,iais1,ibis2,Iamib,ibit1,hb)&
!$OMP& PRIVATE(l0,ibit2,istat,iex,it0,it1,it2,it,factor2)
!$OMP do reduction(+:lprdct)
do l0=1,dim
    ibit2=ishft(iand(list_s(l0), Kis2), -Kisite2)
    istat=2*ibit2+ibit1
    iex = Kml+ Kie*ibit2
    if(llocal_Ham(iex,istat,Kibd)) cycle
    it0=list_s(l0)
    it0=it0+iais1*(Kl1-ibit1)
    it0=it0+ibis2*(Kl2-ibit2)
    it1=iand(it0,irght)
    it2=ishft(iand(it0,ilft),-ihf)
    it=list_lib(1,it1)+list_lib(2,it2)
    factor2 = hb(iex,istat)
    if(iex==istat) factor2 = factor2 + Keperbd
    v0(l0)=v0(l0)+factor2               *vt(it)
    lprdct=lprdct+factor2*dconjg(v1(l0))*vt(it)
enddo
!$OMP END DO
!$OMP END PARALLEL
case(3)
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) &
!$OMP& FIRSTPRIVATE(Kml,Kie,Kibd,Kwt,Kl1,Kl2,Kisite1,Kisite2,&
!$OMP& Kis1,Kis2,Keperbd,irght,ilft,ihf,iais1,ibis2,Iamia,Iamib,hb)&
!$OMP& PRIVATE(l0,ibit1,ibit2,istat,iex,it0,it1,it2,it,factor2)
!$OMP do reduction(+:lprdct)
do l0=1,dim
    ibit1=ishft(iand(list_s(l0), Kis1), -Kisite1)
    ibit2=ishft(iand(list_s(l0), Kis2), -Kisite2)
    istat=2*ibit2+ibit1
    iex = Kml+ Kie*ibit2
    if(llocal_Ham(iex,istat,Kibd)) cycle
    it0=list_s(l0)
    it0=it0+iais1*(Kl1-ibit1)
    it0=it0+ibis2*(Kl2-ibit2)
    it1=iand(it0,irght)
    it2=ishft(iand(it0,ilft),-ihf)
    it=list_lib(1,it1)+list_lib(2,it2)
    factor2 = hb(iex,istat)
    if(iex==istat) factor2 = factor2 + Keperbd
    v0(l0)=v0(l0)+factor2               *vt(it)
    lprdct=lprdct+factor2*dconjg(v1(l0))*vt(it)
enddo
!$OMP END DO
!$OMP END PARALLEL
end select
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
use parameters, only : n, ldim, b_bit, p_bit, max_num_list2, num_list3, parity_type, max_ldim, dsf_comp
use parallel_lib, only : mpi_get_maxint, num_tot_ex_cpus, system_mem_size
use data_mem, only : list1,list3,sz_list_p,sz_list_b
use openmp_lib
implicit none
integer, intent(in) :: tot_sz
integer, intent(in) :: Iam
integer :: icnt, ja, jb, ibpatn, ic
integer :: i, j, isz, ia, ib, is, tmp, psz, pstat, bsz, bstat, i1, i2
integer :: max_a, max_b
integer :: temp_size
!!!! initialization
      ihf0=b_bit/2
      ihf=(b_bit+1)/2
      ihfbit=ishft(1,ihf)
      irght=ihfbit -1
      ilft=ieor(ishft(1,b_bit)-1,irght)
!!
      allocate(sz_list_p(0:ishft(1,p_bit)-1),sz_list_b(0:ishft(1,ihf)-1))
      do j=0,ishft(1,p_bit)-1
        pstat=j
        psz=0
        do i=0,p_bit-1
          psz=psz+iand(pstat,1)
          pstat=ishft(pstat,-1)
        enddo 
        sz_list_p(j)=psz
      enddo
      do j=0,ishft(1,ihf)-1
        bstat=j
        bsz=0
        do i=0,ihf-1
          bsz=bsz+iand(bstat,1)
          bstat=ishft(bstat,-1)
        enddo 
        sz_list_b(j)=bsz
      enddo
!!
      is=0
      if(dsf_comp=='xx') is=1
!!
      icnt=0; ja=0; jb=0; ibpatn=0
      ic = 1; max_a=0; max_b=0
!!
!!!!
!! main loop
psz=sz_list_p(Iam)
do i2=0,ishft(1,ihf0) -1
do i1=0,ishft(1,ihf ) -1
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
write(*,*) 'Size of local vectors,Maximum Size',icnt,max_ldim
write(*,*) 'ME, Maximum size of library',Iam,max_b,max_a
!temp_size=max_ldim/1024
!write(*,*) 'Memoery size estimation (KByte)',16*5*temp_size,4*max_num_list2*num_tot_ex_cpus/1024
!if(system_mem_size< (16*5*temp_size+4*(max_num_list2*num_tot_ex_cpus)/1024) )then
!  write(*,'(A)') 'Error : Memory size'
!  write(*,*) system_mem_size
!  stop
!endif
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
do i2=0,ishft(1,ihf0) -1
do i1=0,ishft(1,ihf ) -1
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
  list1(icnt)=ishft(ib,ihf)+ia
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
use parameters, only : n, max_num_list2, parity_type, dsf_comp
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
is=0; if(dsf_comp=='xx') is=1
!!!!
c_max = num_tot_ex_cpus
!$OMP PARALLEL &
!$OMP DEFAULT(PRIVATE) &
!$OMP SHARED(Iam,c_max,is,target_sz,list3,sz_list_p, &
!$OMP&       sz_list_b,max_num_list2,parity_type,tdim_l,tot_ex_cpus,ihf,ihf0)
!$OMP do SCHEDULE(STATIC,1)
do j0=1, c_max
!!!! initialization
  icnt=0; ja=0 ;jb=0 ;ibpatn=0
  target_cpu = tot_ex_cpus(j0)
!!
psz=sz_list_p(target_cpu)
!!!! main loop
list3(1:2,0:max_num_list2,j0)=-1
!!
do i2=0,ishft(1,ihf0)-1
do i1=0,ishft(1,ihf )-1
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
use parameters, only : max_tdim, parity_type
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
do i2=0,ishft(1,ihf0)-1
do i1=0,ishft(1,ihf )-1
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
use parameters, only : parity_type, max_tdim
use data_mem, only : listt, sz_list_p, sz_list_b
use openmp_lib
implicit none
integer, intent(in) :: target_sz, target_cpu
integer, intent(in) :: Iam
integer, intent(inout) :: pdim
integer :: icnt, isz, psz, i1, i2
allocate(listt(max_tdim))
listt(:)=-1
!!!! main loop
psz=sz_list_p(target_cpu)
icnt=0
do i2=0,ishft(1,ihf0)-1
do i1=0,ishft(1,ihf )-1
  isz=sz_list_b(i1)+sz_list_b(i2)
  if(parity_type==0)then
    if((psz+isz).ne.target_sz) cycle
  endif
  icnt=icnt+1
  listt(icnt)=ishft(i2,ihf)+i1
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
!!
!!
end module parallelized_diag

