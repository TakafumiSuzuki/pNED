!** Numerical exact diagonalization for the Heisenberg-Kitaev model **
! For 4 sub unit
! If you want to change the lattice geometry, modify "set_bond_info"
!!
!  execution
!  %s SIM_TIME [SIZE] [MODEL] [PARAM] [OLD_NP] [COMP] [Q_START Q_END] [RAND INT]'
!
!!
module parameters
implicit none
public
!!
!! See also subroutine set_bond_info and set_xlocal_Ham
!!
  integer :: simple=0
  integer :: c_count=3
  double precision :: xJ1   = 0d0,      xK1   =-1d0
  double precision :: xI11  = 0d0,      xI12  = 0d0
  double precision :: xJ1p  = 0d0,      xK1p  =-1d0
  double precision :: xJ1pp = 0d0,      xK1pp =-1d0
  double precision :: xI11p = 0d0,      xI12p = 0d0
  double precision :: xI11pp= 0d0,      xI12pp= 0d0
  double precision :: xJ2   = 0d0,      xK2   = 0d0
  double precision :: xI21  = 0d0,      xI22  = 0d0
  double precision :: xJ3   = 0d0,      xK3   = 0d0
  double precision :: ratio_l=1d0
  double precision :: alp   = 0d0
!!
!!
double precision, parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164d0
integer,          parameter :: nvec=8
!! Setting of DSF
integer                     :: q_num=20 ! number of Q points
double precision, parameter :: split=2.5d-3 ! Window size along the omega direction
integer         , parameter :: i_range=2000 ! Data points along the omega direction
double precision, parameter :: delt=1d-2 ! Width of Lorentian
!!
!! Setting of iteration number for Lanczos method
integer,parameter          :: max_itr=1024 ! Upper limit of Lanczos iterations
!integer          :: s_maxitr=512 ! Upper limit of DSF Lanczos iterations
!integer          :: s_maxitr=4096 ! Upper limit of DSF Lanczos iterations
integer,parameter          :: s_maxitr=4096 ! Upper limit of DSF Lanczos iterations
!! 
double precision            :: localHam_max=1d2
!!
!!
!! progress; 1:lanc ene 2:lanc vec; 3:cg vec; 4:dsf; 5:recalc DSF, 6:extened (beta), 7: TPS
!logical :: progress(7)=(/.true.,.false.,.true.,.true.,.true.,.false.,.false./)
! for eigen vector
!logical :: progress(7)=(/.true.,.false.,.true.,.false.,.false.,.false.,.false./)
! for DSF
!logical :: progres(7)=(/.false.,.false.,.false.,.true.,.true.,.false.,.false./)
! for calc extra DSF
!logical :: progress(7)=(/.false.,.false.,.false.,.true.,.false.,.true.,.false./)
! for calc TPS
logical :: progress(7)=(/.false.,.false.,.false.,.true.,.false.,.false.,.true./)
!!
!! Computation time (default 23.5 h)
real :: max_comp_time=84000.0
character(LEN=2) :: dsf_comp='zz'
integer :: total_sz
integer :: n, ibond, Lx, Ly
!!
!!
!!
!! trivial parallelization for wave numbers
logical, private, parameter :: wave_parallel=.true.
!! Setting of iteration env.
integer, private :: q_start=1
integer, private :: q_end=nvec
integer, private :: lnc1z_start=2 
integer, private :: lnc1z_end=max_itr
integer, private :: lnc1v_start=2 
integer, private :: lnc1v_end=max_itr
integer, private :: inv1z_start=1 
integer, private :: inv1z_end=5000
integer, private :: cg1_start=1 
integer, private :: cg1_end=500
integer, private :: dsflnc_start=2
integer, private :: dsflnc_end=s_maxitr
!! Etc.
integer :: ldim, tdim, pdim, max_ldim, max_tdim, p_bit, b_bit
character(LEN=128) :: fname
integer :: num_list3=1
logical :: use_conserv=.true.
integer :: max_num_list2
integer :: parity_type
integer :: itr_counter(12)
!!
public set_num_list3, set_itr_counter, set_comp_time, get_cline, set_qrange

contains 
subroutine get_cline(d1,d2,d3,d4,d5,d6,d7,d8,d9)
implicit none
!! If K, the following definision is required.
!integer, external :: iargc
real,    intent(inout), optional :: d1
integer, intent(inout), optional :: d2
integer, intent(inout), optional :: d3
integer, intent(inout), optional :: d4
integer, intent(inout), optional :: d5
integer, intent(inout), optional :: d6
real(KIND(0d0)), intent(inout), optional :: d7
character(LEN=2), intent(inout), optional :: d8
integer, intent(inout), optional :: d9
!!
integer :: numc
character(LEN=128) :: csim_time, clsize, clrange,  cold_np, cqstart, cqend, cpm, ctyp, crand
!!
numc=iargc()
!!
if(numc==9)then 
  call getarg(1,csim_time)
  if(present(d1)) read(csim_time,*) d1
  call getarg(2,clsize)
  if(present(d2)) read(clsize,*) d2
  call getarg(3,clrange)
  if(present(d3)) read(clrange,*) d3
  call getarg(4,cpm)
  if(present(d7)) read(cpm,*) d7
  call getarg(5,cold_np)
  if(present(d4)) read(cold_np,*) d4
  call getarg(6,ctyp)
  if(present(d8)) read(ctyp,*) d8
  call getarg(7,cqstart)
  if(present(d5)) read(cqstart,*) d5
  call getarg(8,cqend)
  if(present(d6)) read(cqend,*) d6
  call getarg(9,crand)
  if(present(d9)) read(crand,*) d9
else
  write(*,'(A)') 'ERROR: %s SIM_TIME [SIZE] [MODEL] [PARAM] [OLD_NP] [COMP] [Q_START Q_END] [RAND INT]'
  write(*,'(A)') 'SIM_TIME : upper limit of the computational time'
  write(*,'(A)') 'SIZE     : system size'
  write(*,'(A)') 'MODEL    : see set_lattice subroutne. The integer variables "ind" corresponds to the value of MODEL.'
  write(*,'(A)') 'PARAM    : tuning parameter for each model'
  write(*,'(A)') 'OLD_NP   : number of flat mpi threads in the previous computation'
  write(*,'(A)') 'COMP     : "ee" for ground state calculation, "zz(xx)" for DSF calculations, and "tp" is Sugiura-Simizu method.'
  stop
endif
write(*,'(F16.4,2I4,F8.2,I4,A,2I4)') d1,d2,d3,d7,d4,d8,d5,d6
end subroutine get_cline

  subroutine set_lattice(ind0,ind,pm)
  implicit none
  integer, intent(in) :: ind,ind0
  real(KIND(0d0)) :: pm
  n=ind0
  total_sz=n/2

!! Original Kitae
  if(ind==0)then
   c_count=3
   simple=0
   alp = pm
   xJ1   = 1d0-alp;  xK1   = 1d0-3d0*alp
   xI11  = 0d0;      xI12  = 0d0
   xJ1p  = 1d0-alp;  xK1p  = 1d0-3d0*alp
   xJ1pp = 1d0-alp;  xK1pp = 1d0-3d0*alp
   xI11p = 0d0;      xI12p = 0d0
   xI11pp= 0d0;      xI12pp= 0d0
   xJ2   = 0d0;      xK2   = 0d0
   xI21  = 0d0;      xI22  = 0d0
   xJ3   = 0d0;      xK3   = 0d0
   ratio_l=1d0
   write(*,*) 'Original'
  elseif(ind==20)then
   c_count=3
   simple=0
   alp = 4d0*atan(1d0)/180d0*pm
   xJ1   = dcos(alp);xK1   = dcos(alp)+2d0*dsin(alp)
   xJ1p  = xJ1;      xK1p  = xK1
   xJ1pp = xJ1;      xK1pp = xK1
   xI11  = 0d0;      xI12  = 0d0
   xI11p = 0d0;      xI12p = 0d0
   xI11pp= 0d0;      xI12pp= 0d0
   xJ2   = 0d0;      xK2   = 0d0
   xI21  = 0d0;      xI22  = 0d0
   xJ3   = 0d0;      xK3   = 0d0
   ratio_l=1d0
   write(*,*) 'Jackeli',alp
  else
    write(*,*) 'original coupling is used.',ind
  endif
  end subroutine set_lattice


subroutine set_comp_time(tm)
implicit none
real, intent(in) :: tm
max_comp_time = tm
end subroutine set_comp_time

subroutine set_qrange(q1,q2,typ)
implicit none
integer, intent(in) :: q1, q2
character(LEN=2), intent(inout) :: typ
q_start=q1
q_end=q2
dsf_comp=typ
if(dsf_comp=='ee')then
  write(*,*) 'Only GS energy and vector is evaluated.'
  progress(:)=(/.true.,.false.,.true.,.false.,.false.,.false.,.false./)
elseif(dsf_comp=='zz')then
  write(*,*) 'Z component of DSF is calculated.'
  progress(:)=(/.false.,.false.,.false.,.true.,.true.,.false.,.false./)
elseif(dsf_comp=='xx')then
  write(*,*) 'X component of DSF is calculated.'
  progress(:)=(/.false.,.false.,.false.,.true.,.true.,.false.,.false./)
elseif(dsf_comp=='zr')then
  write(*,*) 'Z component of DSF is reevaluated for matrix extention.'
  progress(:)=(/.false.,.false.,.false.,.true.,.false.,.true.,.false./)
  dsf_comp='zz'
elseif(dsf_comp=='xr')then
  write(*,*) 'X component of DSF is reevaluated for matrix extention.'
  progress(:)=(/.false.,.false.,.false.,.true.,.false.,.true.,.false./)
  dsf_comp='xx'
elseif(dsf_comp=='tp')then
  write(*,*) 'X component of DSF is reevaluated for matrix extention.'
  progress(:)=(/.false.,.false.,.false.,.false.,.false.,.false.,.true./)
else
  write(*,*) 'Default parallelization conf is used.'
endif


end subroutine set_qrange

subroutine set_num_list3(mm)
implicit none
integer, intent(in) :: mm
num_list3=mm
end subroutine set_num_list3

subroutine set_itr_counter(w_parallel,i_color,old_np,pnp,pme,l_stat)
implicit none
logical,intent(out) :: w_parallel, l_stat
integer,intent(out) :: i_color
integer,intent(in) :: pnp, pme, old_np
integer ::  q_start1, q_end1
integer ::  t1, t2, c_num
if(wave_parallel.and.progress(4).and.(.not.(progress(1))).and.(.not.(progress(2))).and.(.not.(progress(3))))then
    w_parallel=.true.
    t1=max(1,pnp/old_np)
    t2=1
    if(mod(pnp,old_np)==0) t2=0
    c_num=(q_end-q_start+1)/t1 +t2
    if(t1<1)then
      l_stat=.false.
      write(*,*) "np is smaller than the previous computation."
      return
    else
      l_stat=.true.
      i_color=pme/old_np
      q_start1=c_num*i_color+q_start; q_end1=min(q_num,q_start1+c_num-1)
      write(*,*) 'split s',i_color,q_start1,q_end1
    endif
else
    l_stat=.true.
    i_color=0
    w_parallel=.false.
    q_start1=q_start; q_end1=q_end;
endif
itr_counter( 1)=lnc1z_start   ; itr_counter( 2)=lnc1z_end
itr_counter( 3)=lnc1v_start   ; itr_counter( 4)=lnc1v_end
itr_counter( 5)=inv1z_start   ; itr_counter( 6)=inv1z_end
itr_counter( 7)=cg1_start     ; itr_counter( 8)=cg1_end
itr_counter( 9)=q_start1      ; itr_counter(10)=q_end1
itr_counter(11)=dsflnc_start  ; itr_counter(12)=dsflnc_end
end subroutine set_itr_counter
end module parameters

module working_area
use parameters, only : ldim, max_ldim,nvec, itr_counter
implicit none
complex(KIND(0d0)), allocatable :: wwk(:,:) 
complex(KIND(0d0)), allocatable :: rr1(:)
double precision, allocatable :: alpha(:), beta(:), coef(:,:)
double precision, allocatable :: salpha(:), sbeta(:)
complex(KIND(0d0)), allocatable :: xvec(:)
contains
subroutine clear_mem
implicit none
allocate(wwk(max_ldim,4))
allocate(rr1(1:max_ldim))
allocate(alpha(itr_counter(2)), beta(itr_counter(2)), coef(itr_counter(2),nvec+1))
allocate(salpha(1:itr_counter(12)),sbeta(1:itr_counter(12)))
wwk=0.d0
alpha=0.d0
beta=0.d0
coef=0.d0
end subroutine clear_mem

subroutine prep_gsvec(mm)
implicit none
integer, intent(in) :: mm
allocate(xvec(1:mm))
xvec=0d0
end subroutine prep_gsvec

subroutine dealloc_mem
implicit none
if(allocated(wwk))deallocate(wwk)
if(allocated(rr1))deallocate(rr1)
if(allocated(alpha))deallocate(alpha)
if(allocated(beta))deallocate(beta)
if(allocated(coef))deallocate(coef)
if(allocated(salpha))deallocate(salpha)
if(allocated(sbeta))deallocate(sbeta)
if(allocated(xvec))deallocate(xvec)
end subroutine dealloc_mem
end module working_area

module data_mem
use parameters
implicit none
integer, allocatable :: ipair(:), ibond_type(:)
double precision, allocatable :: bondwt(:)
double precision, allocatable :: site_vec(:,:)
double precision, allocatable :: q_list(:,:)
integer, allocatable :: q_list_symmetric(:)
integer, allocatable :: listt(:)
integer, allocatable :: list1(:)
integer, allocatable :: list3(:,:,:)
integer, allocatable :: list_ext3(:,:)
integer, allocatable :: sz_list_p(:),sz_list_b(:)
!double precision :: xlocal_Ham(0:3,0:3,9)
complex(KIND(0d0)) :: xlocal_Ham(0:3,0:3,9)
logical :: llocal_Ham(0:3,0:3,9)

contains
subroutine allocate_lattice_mem(nn,bond,qnum)
implicit none
integer, intent(in) :: nn,bond,qnum
      allocate(ipair(2*bond))
      allocate(ibond_type(bond))
      allocate(site_vec(2,nn))
      allocate(q_list(2,qnum))
      allocate(bondwt(bond))
      allocate(q_list_symmetric(qnum))
      q_list_symmetric(:)=-1
      bondwt(:)=1d0
end subroutine allocate_lattice_mem

subroutine set_bond_info
implicit none
integer :: i, j, k, s, t
double precision :: x, y
double precision :: a1(2), a2(2), a3(2), base(2)
double precision :: k1(2), k2(2), k3(2)
double precision :: k1p(2), k2p(2), k3p(2)
integer :: L1, L2, L3, ic, c_count, iq
!!
c_count=3 ! cordination number
if(simple==1) c_count=8
if(simple==2) c_count=6
if(simple==3) c_count=12
if(simple==4) c_count=9
k1p(:)=0d0
k2p(:)=0d0
k3p(:)=0d0
!! Setting of lattice structure
!!! definition of bonds
!!! definition of unit vectors and recipocal vectors
!!
      a1(1)=0.5d0                  ; a1(2)= 0.5d0/dsqrt(3d0)
      a2(1)=0.5d0                  ; a2(2)=-0.5d0/dsqrt(3d0)
      a3(1)=0d0                    ; a3(2)= 1d0/dsqrt(3d0)
      k1(1)= 2d0*pi                ; k1(2)= 0d0
      k2(1)= 0d0                   ; k2(2)= 2d0*pi/dsqrt(3d0)
      k3(1)= 0d0                   ; k3(2)= 0d0
      k1p(1)= 0d0                  ; k1p(2)= 0d0
      k2p(1)= 0d0                  ; k2p(2)= 0d0
!!!
k=0
iq=0
if(n==32)then
      ibond=n/2*c_count
      q_num=13
      call allocate_lattice_mem(n,ibond,q_num)
      L1=4;  L2=4;
      ic=1
! Nearest neighbor
      do j=0,L2-1
      do i=0,L1-1
! s0
        s = 2*(j*L1+i)+0 +1
! x
        t = 2*(j             *L1+i             )+1 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=1
        ic=ic+1
! y 
        t = 2*(j*L1             +mod(i-1+L1,L1))+1 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=2
        ic=ic+1
! z
        t = 2*(mod(j-1+L2,L2)*L1+i             )+1 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=3
        ic=ic+1
      enddo
      enddo
! 3rd nearest neighbor
      if((simple.ne.0).and.(simple.ne.4))then
      do j=0,L2-1
      do i=0,L1-1
! s0
        s = 2*(j*L1+i)+0 +1
! x
        t = 2*(mod(j-1+L2,L2)*L1+mod(i-1+L1,L1))+0 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=7
        ic=ic+1
! y 
        t = 2*(mod(j-1+L2,L2)*L1+mod(i+1,   L1))+0 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=8
        ic=ic+1
! z
        t = 2*(mod(j+1,L1)   *L1+mod(i-1+L1,L1))+0 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=9
        ic=ic+1
      enddo
      enddo
      endif
! Next nearest neighbor
      if(simple==1.or.simple==3.or.simple==4)then
      do j=0,L2-1
      do i=0,L1-1
! s0
        s = 2*(j*L1+i)+0 +1
!!
      if(simple==3.or.simple==4)then
!
        t = 2*(mod(j-1+L2,L2)*L1+mod(i+1,   L1))+0 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=4
        ic=ic+1
! y 
        t = 2*(mod(j+1   ,L2)*L1+i             )+0 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=5
        ic=ic+1
      endif
! z
        t = 2*(j*L1             +mod(i-1+L1,L1))+0 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=6
        ic=ic+1
! s1
        s = 2*(j*L1+i)+1 +1
      if(simple==3.or.simple==4)then
! x
        t = 2*(mod(j+1,   L2)*L1+mod(i-1+L1,L1))+1 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=4
        ic=ic+1
! y 
        t = 2*(mod(j-1+L2,L2)*L1+i             )+1 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=5
        ic=ic+1
      endif
! z
        t = 2*(j*L1             +mod(i-1+L1,L1))+1 +1
        ipair(2*ic-1)=s ; ipair(2*ic)=t
        ibond_type(ic)=6
        ic=ic+1
      enddo
      enddo
      endif
!!

      base(1)=0.5d0*a3(1)
      base(2)=0.5d0*a3(2)
      k=0
      do j=0,L2-1
      do i=0,L1-1
        k=k+1
        site_vec(1,k)= base(1)        + (a1(1)+a2(1))*dble(i) + (a1(1)+a3(1))*dble(j)
        site_vec(2,k)= base(2)        + (a1(2)+a2(2))*dble(i) + (a1(2)+a3(2))*dble(j)
        k=k+1
        site_vec(1,k)= base(1)+ a1(1) + (a1(1)+a2(1))*dble(i) + (a1(1)+a3(1))*dble(j)
        site_vec(2,k)= base(2)+ a1(2) + (a1(2)+a2(2))*dble(i) + (a1(2)+a3(2))*dble(j)
      enddo
      enddo
!!
      k=0
      do j= 0,2
      do i= 0,2
        k=k+1
        q_list(1,k)=k1(1)/dble(2)*dble(i)+k2(1)/dble(2)*dble(j);
        q_list(2,k)=k1(2)/dble(2)*dble(i)+k2(2)/dble(2)*dble(j);
        if (mod(i,2)==0.and.mod(j,2)==0)then
          iq=iq+1
          q_list_symmetric(iq)=k
        endif
      enddo
      enddo
      do j= 0,1
      do i= 0,1
        k=k+1
        q_list(1,k)=k1(1)/dble(2)*dble(i)+k2(1)/dble(2)*dble(j) + (k1(1) + k2(1))/dble(4);
        q_list(2,k)=k1(2)/dble(2)*dble(i)+k2(2)/dble(2)*dble(j) + (k1(2) + k2(2))/dble(4);
        if (i==1.and.j==1)then
          iq=iq+1
          q_list_symmetric(iq)=k
        endif
      enddo
      enddo

!!! Symmetric
k=0
elseif(n==24)then
  ibond=c_count*n/2
  q_num=12
  call allocate_lattice_mem(n,ibond,q_num)
  L1=6; L2=2; L3=2;
  if(simple==0)then
    ipair= (/ &
! 1st N. x-bond
          1, 2,  3, 4,  5, 6,  7, 8,  9,10, 11,12, &
         13,14, 15,16, 17,18, 19,20, 21,22, 23,24, &
! 1st N. z-bond
          1,18,  2,13,  3,20,  4,15,  5,22,  6,17, &
          7,24,  8,19,  9,14, 10,21, 11,16, 12,23, &
! 1st N. y-bond
          1,12,  3, 2,  5, 4,  7, 6,  9, 8, 11,10, &
         13,24, 15,14, 17,16, 19,18, 21,20, 23,22 &
       /)
! type _Z
    ibond_type=(/ &
                1,1,1,1,1,1,1,1,1,1,1,1,  &
                3,3,3,3,3,3,3,3,3,3,3,3,  &
                2,2,2,2,2,2,2,2,2,2,2,2  &
            /)
!!
  elseif(simple==1)then
    ipair= (/ &
! 1st N. x-bond
          1, 2,  3, 4,  5, 6,  7, 8,  9,10, 11,12, &
         13,14, 15,16, 17,18, 19,20, 21,22, 23,24, &
! 1st N. z-bond
          1,18,  2,13,  3,20,  4,15,  5,22,  6,17, &
          7,24,  8,19,  9,14, 10,21, 11,16, 12,23, &
! 1st N. y-bond
          1,12,  3, 2,  5, 4,  7, 6,  9, 8, 11,10, &
         13,24, 15,14, 17,16, 19,18, 21,20, 23,22, &
! 2nd N. X-bond  
!          1,23,  2,24,  3,13,  4,14,  5,15,  6,16, &
!          7,17,  8,18,  9,19, 10,20, 11,21, 12,22, &
!          1,19,  2,20,  3,21,  4,22,  5,23,  6,24, &
!          7,13,  8,14,  9,15, 10,16, 11,17, 12,18, &
! 2nd N. Z-bond  
          1, 3,  2, 4,  3, 5,  4, 6,  5, 7,  6, 8, &
          7, 9,  8,10,  9,11, 10,12, 11, 1, 12, 2, &
         13,15, 14,16, 15,17, 16,18, 17,19, 18,20, &
         19,21, 20,22, 21,23, 22,24, 23,13, 24,14, &
! 2nd N. Y-bond  
!          1,13,  2,14,  3,15,  4,16,  5,17,  6,18, &
!          7,19,  8,20,  9,21, 10,22, 11,23, 12,24, &
!          1,17,  2,18,  3,19,  4,20,  5,21,  6,22, &
!          7,23,  8,24,  9,13, 10,14, 11,15, 12,16, &
! 3rd N. x-bond
          1,16,  2,15,  3,18,  4,17,  5,20,  6,19, &
          7,22,  8,21,  9,24, 10,23, 11,14, 12,13, &
! 3rd N. z-bond
          1,24,  2,19,  3,14,  4,21,  5,16,  6,23, &
          7,18,  8,13,  9,20, 10,15, 11,22, 12,17, &
! 3rd N. y-bond
          1,20,  2,23,  3,22,  4,13,  5,24,  6,15, &
          7,14,  8,17,  9,16, 10,19, 11,18, 12,21  &
       /)
! type _Z
    ibond_type=(/ &
                1,1,1,1,1,1,1,1,1,1,1,1,  &
                3,3,3,3,3,3,3,3,3,3,3,3,  &
                2,2,2,2,2,2,2,2,2,2,2,2,  &
!                4,4,4,4,4,4,4,4,4,4,4,4,  &
!                4,4,4,4,4,4,4,4,4,4,4,4,  &
                6,6,6,6,6,6,6,6,6,6,6,6,  &
                6,6,6,6,6,6,6,6,6,6,6,6,  &
!                5,5,5,5,5,5,5,5,5,5,5,5,  &
!                5,5,5,5,5,5,5,5,5,5,5,5,  &
                7,7,7,7,7,7,7,7,7,7,7,7,  &
                9,9,9,9,9,9,9,9,9,9,9,9,  &
                8,8,8,8,8,8,8,8,8,8,8,8   &
            /)
  elseif(simple==2)then
    ipair= (/ &
! 1st N. x-bond
          1, 2,  3, 4,  5, 6,  7, 8,  9,10, 11,12, &
         13,14, 15,16, 17,18, 19,20, 21,22, 23,24, &
! 1st N. z-bond
          1,18,  2,13,  3,20,  4,15,  5,22,  6,17, &
          7,24,  8,19,  9,14, 10,21, 11,16, 12,23, &
! 1st N. y-bond
          1,12,  3, 2,  5, 4,  7, 6,  9, 8, 11,10, &
         13,24, 15,14, 17,16, 19,18, 21,20, 23,22, &
! 2nd N. X-bond  
!          1,23,  2,24,  3,13,  4,14,  5,15,  6,16, &
!          7,17,  8,18,  9,19, 10,20, 11,21, 12,22, &
!          1,19,  2,20,  3,21,  4,22,  5,23,  6,24, &
!          7,13,  8,14,  9,15, 10,16, 11,17, 12,18, &
! 2nd N. Z-bond  
!          1, 3,  2, 4,  3, 5,  4, 6,  5, 7,  6, 8, &
!          7, 9,  8,10,  9,11, 10,12, 11, 1, 12, 2, &
!         13,15, 14,16, 15,17, 16,18, 17,19, 18,20, &
!         19,21, 20,22, 21,23, 22,24, 23,13, 24,14, &
! 2nd N. Y-bond  
!          1,13,  2,14,  3,15,  4,16,  5,17,  6,18, &
!          7,19,  8,20,  9,21, 10,22, 11,23, 12,24, &
!          1,17,  2,18,  3,19,  4,20,  5,21,  6,22, &
!          7,23,  8,24,  9,13, 10,14, 11,15, 12,16, &
! 3rd N. x-bond
          1,16,  2,15,  3,18,  4,17,  5,20,  6,19, &
          7,22,  8,21,  9,24, 10,23, 11,14, 12,13, &
! 3rd N. z-bond
          1,24,  2,19,  3,14,  4,21,  5,16,  6,23, &
          7,18,  8,13,  9,20, 10,15, 11,22, 12,17, &
! 3rd N. y-bond
          1,20,  2,23,  3,22,  4,13,  5,24,  6,15, &
          7,14,  8,17,  9,16, 10,19, 11,18, 12,21  &
       /)
! type _Z
    ibond_type=(/ &
                1,1,1,1,1,1,1,1,1,1,1,1,  &
                3,3,3,3,3,3,3,3,3,3,3,3,  &
                2,2,2,2,2,2,2,2,2,2,2,2,  &
!                4,4,4,4,4,4,4,4,4,4,4,4,  &
!                4,4,4,4,4,4,4,4,4,4,4,4,  &
!                6,6,6,6,6,6,6,6,6,6,6,6,  &
!                6,6,6,6,6,6,6,6,6,6,6,6,  &
!                5,5,5,5,5,5,5,5,5,5,5,5,  &
!                5,5,5,5,5,5,5,5,5,5,5,5,  &
                7,7,7,7,7,7,7,7,7,7,7,7,  &
                9,9,9,9,9,9,9,9,9,9,9,9,  &
                8,8,8,8,8,8,8,8,8,8,8,8   &
            /)
  elseif(simple==3)then
    ipair= (/ &
! 1st N. x-bond
          1, 2,  3, 4,  5, 6,  7, 8,  9,10, 11,12, &
         13,14, 15,16, 17,18, 19,20, 21,22, 23,24, &
! 1st N. z-bond
          1,18,  2,13,  3,20,  4,15,  5,22,  6,17, &
          7,24,  8,19,  9,14, 10,21, 11,16, 12,23, &
! 1st N. y-bond
          1,12,  3, 2,  5, 4,  7, 6,  9, 8, 11,10, &
         13,24, 15,14, 17,16, 19,18, 21,20, 23,22, &
! 2nd N. X-bond  
          1,23,  2,24,  3,13,  4,14,  5,15,  6,16, &
          7,17,  8,18,  9,19, 10,20, 11,21, 12,22, &
          1,19,  2,20,  3,21,  4,22,  5,23,  6,24, &
          7,13,  8,14,  9,15, 10,16, 11,17, 12,18, &
! 2nd N. Z-bond  
          1, 3,  2, 4,  3, 5,  4, 6,  5, 7,  6, 8, &
          7, 9,  8,10,  9,11, 10,12, 11, 1, 12, 2, &
         13,15, 14,16, 15,17, 16,18, 17,19, 18,20, &
         19,21, 20,22, 21,23, 22,24, 23,13, 24,14, &
! 2nd N. Y-bond  
          1,13,  2,14,  3,15,  4,16,  5,17,  6,18, &
          7,19,  8,20,  9,21, 10,22, 11,23, 12,24, &
          1,17,  2,18,  3,19,  4,20,  5,21,  6,22, &
          7,23,  8,24,  9,13, 10,14, 11,15, 12,16, &
! 3rd N. x-bond
          1,16,  2,15,  3,18,  4,17,  5,20,  6,19, &
          7,22,  8,21,  9,24, 10,23, 11,14, 12,13, &
! 3rd N. z-bond
          1,24,  2,19,  3,14,  4,21,  5,16,  6,23, &
          7,18,  8,13,  9,20, 10,15, 11,22, 12,17, &
! 3rd N. y-bond
          1,20,  2,23,  3,22,  4,13,  5,24,  6,15, &
          7,14,  8,17,  9,16, 10,19, 11,18, 12,21  &
       /)
! type _Z
    ibond_type=(/ &
                1,1,1,1,1,1,1,1,1,1,1,1,  &
                3,3,3,3,3,3,3,3,3,3,3,3,  &
                2,2,2,2,2,2,2,2,2,2,2,2,  &
                4,4,4,4,4,4,4,4,4,4,4,4,  &
                4,4,4,4,4,4,4,4,4,4,4,4,  &
                6,6,6,6,6,6,6,6,6,6,6,6,  &
                6,6,6,6,6,6,6,6,6,6,6,6,  &
                5,5,5,5,5,5,5,5,5,5,5,5,  &
                5,5,5,5,5,5,5,5,5,5,5,5,  &
                7,7,7,7,7,7,7,7,7,7,7,7,  &
                9,9,9,9,9,9,9,9,9,9,9,9,  &
                8,8,8,8,8,8,8,8,8,8,8,8   &
            /)
  elseif(simple==4)then
    ipair= (/ &
! 1st N. x-bond
          1, 2,  3, 4,  5, 6,  7, 8,  9,10, 11,12, &
         13,14, 15,16, 17,18, 19,20, 21,22, 23,24, &
! 1st N. z-bond
          1,18,  2,13,  3,20,  4,15,  5,22,  6,17, &
          7,24,  8,19,  9,14, 10,21, 11,16, 12,23, &
! 1st N. y-bond
          1,12,  3, 2,  5, 4,  7, 6,  9, 8, 11,10, &
         13,24, 15,14, 17,16, 19,18, 21,20, 23,22, &
! 2nd N. X-bond  
          1,23,  2,24,  3,13,  4,14,  5,15,  6,16, &
          7,17,  8,18,  9,19, 10,20, 11,21, 12,22, &
          1,19,  2,20,  3,21,  4,22,  5,23,  6,24, &
          7,13,  8,14,  9,15, 10,16, 11,17, 12,18, &
! 2nd N. Z-bond  
          1, 3,  2, 4,  3, 5,  4, 6,  5, 7,  6, 8, &
          7, 9,  8,10,  9,11, 10,12, 11, 1, 12, 2, &
         13,15, 14,16, 15,17, 16,18, 17,19, 18,20, &
         19,21, 20,22, 21,23, 22,24, 23,13, 24,14, &
! 2nd N. Y-bond  
          1,13,  2,14,  3,15,  4,16,  5,17,  6,18, &
          7,19,  8,20,  9,21, 10,22, 11,23, 12,24, &
          1,17,  2,18,  3,19,  4,20,  5,21,  6,22, &
          7,23,  8,24,  9,13, 10,14, 11,15, 12,16 &
! 3rd N. x-bond
!          1,16,  2,15,  3,18,  4,17,  5,20,  6,19, &
!          7,22,  8,21,  9,24, 10,23, 11,14, 12,13, &
! 3rd N. z-bond
!          1,24,  2,19,  3,14,  4,21,  5,16,  6,23, &
!          7,18,  8,13,  9,20, 10,15, 11,22, 12,17, &
! 3rd N. y-bond
!          1,20,  2,23,  3,22,  4,13,  5,24,  6,15, &
!          7,14,  8,17,  9,16, 10,19, 11,18, 12,21  &
       /)
! type _Z
    ibond_type=(/ &
                1,1,1,1,1,1,1,1,1,1,1,1,  &
                3,3,3,3,3,3,3,3,3,3,3,3,  &
                2,2,2,2,2,2,2,2,2,2,2,2,  &
                4,4,4,4,4,4,4,4,4,4,4,4,  &
                4,4,4,4,4,4,4,4,4,4,4,4,  &
                6,6,6,6,6,6,6,6,6,6,6,6,  &
                6,6,6,6,6,6,6,6,6,6,6,6,  &
                5,5,5,5,5,5,5,5,5,5,5,5,  &
                5,5,5,5,5,5,5,5,5,5,5,5  &
!                7,7,7,7,7,7,7,7,7,7,7,7,  &
!                9,9,9,9,9,9,9,9,9,9,9,9,  &
!                8,8,8,8,8,8,8,8,8,8,8,8   &
            /)
  end if
!!
  base(1)=0.5d0*a3(1)
  base(2)=0.5d0*a3(2)
  do i=1,12
    if(iand(i,1)==1)then
      site_vec(1,i)= base(1) + (a1(1)+a2(1))*0.5d0*dble(i-1)
      site_vec(2,i)= base(2) + (a1(2)+a2(2))*0.5d0*dble(i-1)
    else
      site_vec(1,i)= site_vec(1,i-1)+a1(1)
      site_vec(2,i)= site_vec(2,i-1)+a1(2)
    endif
  enddo
  base(1)=a1(1)+a3(1)+site_vec(1,1)
  base(2)=a1(2)+a3(2)+site_vec(2,1)
  do i=13,24
    if(iand(i,1)==1)then
      site_vec(1,i)= base(1) + (a1(1)+a2(1))*0.5d0*dble(i-13)
      site_vec(2,i)= base(2) + (a1(2)+a2(2))*0.5d0*dble(i-13)
    else
      site_vec(1,i)= site_vec(1,i-1)+a1(1)
      site_vec(2,i)= site_vec(2,i-1)+a1(2)
    endif
  enddo
  do j = 0,2,2
  do i = 0,3
    k=k+1
    q_list(1,k)=k1(1)/dble(3)*dble(i)+k2(1)/dble(2)*dble(j);
    q_list(2,k)=k1(2)/dble(3)*dble(i)+k2(2)/dble(2)*dble(j);
    if ((i==0.or.i==3).and.(j==0.or.j==2))then
      iq=iq+1
      q_list_symmetric(iq)=k
    endif
  enddo
  enddo
  do j= 1,1,2
  do i= 1,5,2
    k=k+1
    q_list(1,k)=k1(1)/dble(2*3)*dble(i)+k2(1)/dble(2)*dble(j);
    q_list(2,k)=k1(2)/dble(2*3)*dble(i)+k2(2)/dble(2)*dble(j);
    if (i==3.and.j==1)then
      iq=iq+1
      q_list_symmetric(iq)=k
    endif
  enddo
  enddo
!!
elseif(n==8)then
    ibond=n/2*c_count
    q_num=9
    call allocate_lattice_mem(n,ibond,q_num)
    if(simple==0)then
    write(*,*) 'size 8 simple'
    ipair= (/ &
            1,2,  3,7,  5,4,  6,8, &
            7,1,  2,3,  6,5,  4,8, &
            1,6,  8,2,  3,4,  5,7  &
           /)
    ibond_type=(/ &
               1,1,1,1,  &
               2,2,2,2,  &
               3,3,3,3   &
             /)
    elseif(simple==1)then
    ipair= (/ &
            1,2,  3,7,  5,4,  6,8, &
            7,1,  2,3,  6,5,  4,8, &
            1,6,  8,2,  3,4,  5,7, &
!            1,5,  1,5,  2,4,  2,4,  3,8,  3,8,  6,7,  6,7, &
!            1,8,  1,8,  2,6,  2,6,  3,5,  3,5,  4,7,  4,7, &
            1,3,  1,3,  2,7,  2,7,  4,6,  4,6,  5,8,  5,8, &
            1,4,  2,5,  3,6,  7,8,  &
            1,4,  2,5,  3,6,  7,8,  &
            1,4,  2,5,  3,6,  7,8   &
           /)
    ibond_type=(/ &
               1,1,1,1,  &
               2,2,2,2,  &
               3,3,3,3,  &
!               4,4,4,4, 4,4,4,4, &
!               5,5,5,5, 5,5,5,5, &
               6,6,6,6, 6,6,6,6, &
               7,7,7,7,  &
               8,8,8,8,  &
               9,9,9,9   &
             /)
    elseif(simple==2)then
    ipair= (/ &
            1,2,  3,7,  5,4,  6,8, &
            7,1,  2,3,  6,5,  4,8, &
            1,6,  8,2,  3,4,  5,7, &
!            1,5,  1,5,  2,4,  2,4,  3,8,  3,8,  6,7,  6,7, &
!            1,8,  1,8,  2,6,  2,6,  3,5,  3,5,  4,7,  4,7, &
!            1,3,  1,3,  2,7,  2,7,  4,6,  4,6,  5,8,  5,8, &
            1,4,  2,5,  3,6,  7,8,  &
            1,4,  2,5,  3,6,  7,8,  &
            1,4,  2,5,  3,6,  7,8   &
           /)
    ibond_type=(/ &
               1,1,1,1,  &
               2,2,2,2,  &
               3,3,3,3,  &
!               4,4,4,4, 4,4,4,4, &
!               5,5,5,5, 5,5,5,5, &
!               6,6,6,6, 6,6,6,6, &
               7,7,7,7,  &
               8,8,8,8,  &
               9,9,9,9   &
             /)

    elseif(simple==3)then
    ipair= (/ &
            1,2,  3,7,  5,4,  6,8, &
            7,1,  2,3,  6,5,  4,8, &
            1,6,  8,2,  3,4,  5,7, &
            1,5,  1,5,  2,4,  2,4,  3,8,  3,8,  6,7,  6,7, &
            1,8,  1,8,  2,6,  2,6,  3,5,  3,5,  4,7,  4,7, &
            1,3,  1,3,  2,7,  2,7,  4,6,  4,6,  5,8,  5,8, &
            1,4,  2,5,  3,6,  7,8,  &
            1,4,  2,5,  3,6,  7,8,  &
            1,4,  2,5,  3,6,  7,8   &
           /)
    ibond_type=(/ &
               1,1,1,1,  &
               2,2,2,2,  &
               3,3,3,3,  &
               4,4,4,4, 4,4,4,4, &
               5,5,5,5, 5,5,5,5, &
               6,6,6,6, 6,6,6,6, &
               7,7,7,7,  &
               8,8,8,8,  &
               9,9,9,9   &
             /)

    elseif(simple==4)then
    ipair= (/ &
            1,2,  3,7,  5,4,  6,8, &
            7,1,  2,3,  6,5,  4,8, &
            1,6,  8,2,  3,4,  5,7, &
            1,5,  1,5,  2,4,  2,4,  3,8,  3,8,  6,7,  6,7, &
            1,8,  1,8,  2,6,  2,6,  3,5,  3,5,  4,7,  4,7, &
            1,3,  1,3,  2,7,  2,7,  4,6,  4,6,  5,8,  5,8 &
!            1,4,  2,5,  3,6,  7,8,  &
!            1,4,  2,5,  3,6,  7,8,  &
!            1,4,  2,5,  3,6,  7,8   &
           /)
    ibond_type=(/ &
               1,1,1,1,  &
               2,2,2,2,  &
               3,3,3,3,  &
               4,4,4,4, 4,4,4,4, &
               5,5,5,5, 5,5,5,5, &
               6,6,6,6, 6,6,6,6 &
!               7,7,7,7,  &
!               8,8,8,8,  &
!               9,9,9,9   &
             /)
    endif 
!!
    site_vec(1,1)= 1.5d0*a3(1)+a1(1)
    site_vec(2,1)= 1.5d0*a3(2)+a1(2)
    site_vec(1,2)= site_vec(1,1)+a1(1)
    site_vec(2,2)= site_vec(2,1)+a1(2)
    site_vec(1,3)= site_vec(1,2)+a2(1)
    site_vec(2,3)= site_vec(2,2)+a2(2)
    site_vec(1,4)= site_vec(1,3)-a3(1)
    site_vec(2,4)= site_vec(2,3)-a3(2)
    site_vec(1,5)= site_vec(1,4)-a1(1)
    site_vec(2,5)= site_vec(2,4)-a1(2)
    site_vec(1,6)= site_vec(1,5)-a2(1)
    site_vec(2,6)= site_vec(2,5)-a2(2)
    site_vec(1,7)= site_vec(1,5)-a3(1)
    site_vec(2,7)= site_vec(2,5)-a3(2)
    site_vec(1,8)= site_vec(1,2)+a3(1)
    site_vec(2,8)= site_vec(2,2)+a3(2)
    k=0
    do j = 0,2,2
    do i = 0,2,2
      k=k+1
      q_list(1,k)=k1(1)/dble(2)*dble(i)+k2(1)/dble(2)*dble(j);
      q_list(2,k)=k1(2)/dble(2)*dble(i)+k2(2)/dble(2)*dble(j);
      iq=iq+1
      q_list_symmetric(iq)=k
    enddo
    enddo
    do j= 1,1,2
    do i= 1,1,2
      k=k+1
      q_list(1,k)=k1(1)/dble(2)*dble(i)+k2(1)/dble(2)*dble(j);
      q_list(2,k)=k1(2)/dble(2)*dble(i)+k2(2)/dble(2)*dble(j);
      iq=iq+1
      q_list_symmetric(iq)=k
    enddo
    enddo

    elseif(n==12)then
    ibond=n/2*c_count
    q_num=8
    call allocate_lattice_mem(n,ibond,q_num)
    if(simple==0)then
    ipair= (/ &
            1,2, 5,6, 9,10, 3,4, 7,8, 11,12,&
            2,5, 6,9, 10,1, 4,7, 8,11, 12,3,&
            1,4, 5,8, 9,12, 3,2, 7,6, 11,10 &
           /)
  ibond_type=(/ &
               1,1,1,1,1,1,  &
               2,2,2,2,2,2,  &
               3,3,3,3,3,3   &
             /)
    elseif(simple==1)then
    ipair= (/ &
            1,2, 5,6, 9,10, 3,4, 7,8, 11,12,&
            2,5, 6,9, 10,1, 4,7, 8,11, 12,3,&
            1,4, 5,8, 9,12, 3,2, 7,6, 11,10,&
!            1,11, 1,7, 2,12, 2,8, 5,3, 5,11, 6,4, 6,12, 9,7, 9,3, 10,8, 10,4,&
!            1,3, 1,3, 2,4, 2,4, 5,7, 5,7, 6,8, 6,8, 9,11, 9,11, 10,12, 10,12,&
            1,9, 2,10, 5,1, 6,2, 9,5, 10,6, 3,11, 4,12, 7,3, 8,4, 11,7, 12,8,&
            1,12, 2,7, 5,4, 6,11, 9,8, 10,3,&
            1,8, 2,11, 5,12, 6,3, 9,4, 10,7,&
            1,12, 2,7, 5,4, 6,11, 9,8, 10,3 &
           /)
  ibond_type=(/ &
               1,1,1,1,1,1,  &
               2,2,2,2,2,2,  &
               3,3,3,3,3,3,  &
!               4,4,4,4,4,4, 4,4,4,4,4,4,  &
!               5,5,5,5,5,5, 5,5,5,5,5,5,  &
               6,6,6,6,6,6, 6,6,6,6,6,6,  &
               7,7,7,7,7,7,  &
               8,8,8,8,8,8,  &
               9,9,9,9,9,9   &
             /)
    elseif(simple==2)then
    ipair= (/ &
            1,2, 5,6, 9,10, 3,4, 7,8, 11,12,&
            2,5, 6,9, 10,1, 4,7, 8,11, 12,3,&
            1,4, 5,8, 9,12, 3,2, 7,6, 11,10,&
!            1,11, 1,7, 2,12, 2,8, 5,3, 5,11, 6,4, 6,12, 9,7, 9,3, 10,8, 10,4,&
!            1,3, 1,3, 2,4, 2,4, 5,7, 5,7, 6,8, 6,8, 9,11, 9,11, 10,12, 10,12,&
!            1,9, 2,10, 5,1, 6,2, 9,5, 10,6, 3,11, 4,12, 7,3, 8,4, 11,7, 12,8,&
            1,12, 2,7, 5,4, 6,11, 9,8, 10,3,&
            1,8, 2,11, 5,12, 6,3, 9,4, 10,7,&
            1,12, 2,7, 5,4, 6,11, 9,8, 10,3 &
           /)
  ibond_type=(/ &
               1,1,1,1,1,1,  &
               2,2,2,2,2,2,  &
               3,3,3,3,3,3,  &
!               4,4,4,4,4,4, 4,4,4,4,4,4,  &
!               5,5,5,5,5,5, 5,5,5,5,5,5,  &
!               6,6,6,6,6,6, 6,6,6,6,6,6,  &
               7,7,7,7,7,7,  &
               8,8,8,8,8,8,  &
               9,9,9,9,9,9   &
             /)
    elseif(simple==3)then
    ipair= (/ &
            1,2, 5,6, 9,10, 3,4, 7,8, 11,12,&
            2,5, 6,9, 10,1, 4,7, 8,11, 12,3,&
            1,4, 5,8, 9,12, 3,2, 7,6, 11,10,&
            1,11, 1,7, 2,12, 2,8, 5,3, 5,11, 6,4, 6,12, 9,7, 9,3, 10,8, 10,4,&
            1,3, 1,3, 2,4, 2,4, 5,7, 5,7, 6,8, 6,8, 9,11, 9,11, 10,12, 10,12,&
            1,9, 2,10, 5,1, 6,2, 9,5, 10,6, 3,11, 4,12, 7,3, 8,4, 11,7, 12,8,&
            1,12, 2,7, 5,4, 6,11, 9,8, 10,3,&
            1,8, 2,11, 5,12, 6,3, 9,4, 10,7,&
            1,12, 2,7, 5,4, 6,11, 9,8, 10,3 &
           /)
  ibond_type=(/ &
               1,1,1,1,1,1,  &
               2,2,2,2,2,2,  &
               3,3,3,3,3,3,  &
               4,4,4,4,4,4, 4,4,4,4,4,4,  &
               5,5,5,5,5,5, 5,5,5,5,5,5,  &
               6,6,6,6,6,6, 6,6,6,6,6,6,  &
               7,7,7,7,7,7,  &
               8,8,8,8,8,8,  &
               9,9,9,9,9,9   &
             /)
    elseif(simple==4)then
    ipair= (/ &
            1,2, 5,6, 9,10, 3,4, 7,8, 11,12,&
            2,5, 6,9, 10,1, 4,7, 8,11, 12,3,&
            1,4, 5,8, 9,12, 3,2, 7,6, 11,10,&
            1,11, 1,7, 2,12, 2,8, 5,3, 5,11, 6,4, 6,12, 9,7, 9,3, 10,8, 10,4,&
            1,3, 1,3, 2,4, 2,4, 5,7, 5,7, 6,8, 6,8, 9,11, 9,11, 10,12, 10,12,&
            1,9, 2,10, 5,1, 6,2, 9,5, 10,6, 3,11, 4,12, 7,3, 8,4, 11,7, 12,8&
!            1,12, 2,7, 5,4, 6,11, 9,8, 10,3,&
!            1,8, 2,11, 5,12, 6,3, 9,4, 10,7,&
!            1,12, 2,7, 5,4, 6,11, 9,8, 10,3 &
           /)
  ibond_type=(/ &
               1,1,1,1,1,1,  &
               2,2,2,2,2,2,  &
               3,3,3,3,3,3,  &
               4,4,4,4,4,4, 4,4,4,4,4,4,  &
               5,5,5,5,5,5, 5,5,5,5,5,5,  &
               6,6,6,6,6,6, 6,6,6,6,6,6  &
!               7,7,7,7,7,7,  &
!               8,8,8,8,8,8,  &
!               9,9,9,9,9,9   &
             /)
    endif
!!
    k=0
    do i=1,3
      base(1)=0.5d0*a3(1) + dble(i-1)*(a1(1)+a2(1))
      base(2)=0.5d0*a3(2) + dble(i-1)*(a1(2)+a2(2))
      k=k+1
        site_vec(1,k)= base(1) 
        site_vec(2,k)= base(2) 
      k=k+1
        site_vec(1,k)= site_vec(1,k-1) + a1(1)
        site_vec(2,k)= site_vec(2,k-1) + a1(2)
      k=k+1
        site_vec(1,k)= site_vec(1,k-1) + a3(1)
        site_vec(2,k)= site_vec(2,k-1) + a3(2)
      k=k+1
        site_vec(1,k)= site_vec(1,k-1) + a1(1)
        site_vec(2,k)= site_vec(2,k-1) + a1(2)
    enddo
    k=0
    do j = 0,2,2
    do i = 0,3
      k=k+1
      q_list(1,k)=k1(1)/dble(3)*dble(i)+k2(1)/dble(2)*dble(j);
      q_list(2,k)=k1(2)/dble(3)*dble(i)+k2(2)/dble(2)*dble(j);
      if((i==0.or.i==3).and.(j==0.or.j==2))then
        iq=iq+1
        q_list_symmetric(iq)=k
      endif
    enddo
    enddo
!    do j= 1,1,2
!    do i= 1,5,2
!      k=k+1
!      q_list(1,k)=k1(1)/dble(6)*dble(i)+k2(1)/dble(2)*dble(j);
!      q_list(2,k)=k1(2)/dble(6)*dble(i)+k2(2)/dble(2)*dble(j);
!    enddo
!    enddo
    elseif(n==18)then
    ibond=n/2*c_count
    q_num=8
    allocate(ipair(2*ibond), ibond_type(ibond))
    allocate(bondwt(ibond))
    allocate(site_vec(2,n))
    allocate(q_list(2,q_num))
    bondwt(:)=1d0
    if(simple==0)then
    ipair= (/ &
           1,2, 3,4, 5,6, 7,8, 9,10, 11,12, 13,14, 15,16, 17,18, &
           6,1, 2,3, 4,5, 12,7, 8,9, 10,11, 18,13, 14,15, 16,17, &
           1,14, 3,16, 5,18, 7,2, 9,4, 11,6, 13,8, 15,10, 17,12  &
           /)
    ibond_type=(/ &
               1,1,1,1,1,1,1,1,1,  &
               2,2,2,2,2,2,2,2,2,  &
               3,3,3,3,3,3,3,3,3   &
             /)
    else
    ipair= (/ &
           1,2, 3,4, 5,6, 7,8, 9,10, 11,12, 13,14, 15,16, 17,18, &
           6,1, 2,3, 4,5, 12,7, 8,9, 10,11, 18,13, 14,15, 16,17, &
           1,14, 3,16, 5,18, 7,2, 9,4, 11,6, 13,8, 15,10, 17,12, &
!           1,11,  2,12,  3,7,  4,8,  5,9,  6,10,  7,17,  8,18,  9,13,  &
!           10,14,  11,15,  12,16,  13,5,  14,6,  15,1,  16,2,  17,3,  18,4,  &
!           1,7,  2,8,  3,9,  4,10,  5,11,  6,12,  7,13,  8,14,  9,15,  &
!           10,16,  11,17,  12,18,  13,1,  14,2,  15,3,  16,4,  17,5,  18,6,  &
           1,3,  2,4,  3,5,  4,6,  5,1,  6,2,  7,9,  8,10,  9,11,  10,12,  &
           11,7,  12,8,  13,15,  14,16,  15,17,  16,18,  17,13,  18,14,  &
           1,18,  3,14,  5,17,  7,6,  9,2,  11,4,  13,12,  15,8,  17,10,  &
           1,15,  3,17,  5,14,  7,4,  9,6,  11,2,  13,10,  15,12,  17,7,  &
           1,12,  3,8,  5,10,  7,18,  9,14,  11,16,  13,6,  15,2,  17,4   &
           /)
! type _Z
    ibond_type=(/ &
               1,1,1,1,1,1,1,1,1,  &
               2,2,2,2,2,2,2,2,2,  &
               3,3,3,3,3,3,3,3,3,  &
!               4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,  &
!               5,5,5,5,5,5,5,5,5, 5,5,5,5,5,5,5,5,5,  &
               6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,  &
               7,7,7,7,7,7,7,7,7,  &
               8,8,8,8,8,8,8,8,8,  &
               9,9,9,9,9,9,9,9,9   &
             /)
    endif
!!
    base(1)=0.5d0*a3(1)
    base(2)=0.5d0*a3(2)
    do i=1,6
      k=k+1
      if(mod(i,2)==1)then
        site_vec(1,k)= base(1) + 0.5d0*dble(i-1)*(a1(1)+a2(1))
        site_vec(2,k)= base(2) + 0.5d0*dble(i-1)*(a1(2)+a2(2))
      else
        site_vec(1,k)= site_vec(1,k-1)+a1(1)
        site_vec(2,k)= site_vec(2,k-1)+a1(2)
      endif
    enddo
    base(1)=1.5d0*a3(1)+a1(1)
    base(2)=1.5d0*a3(2)+a1(2)
    do i=1,6
      k=k+1
      if(mod(i,2)==1)then
        site_vec(1,k)= base(1) + 0.5d0*dble(i-1)*(a1(1)+a2(1))
        site_vec(2,k)= base(2) + 0.5d0*dble(i-1)*(a1(2)+a2(2))
      else
        site_vec(1,k)= site_vec(1,k-1)+a1(1)
        site_vec(2,k)= site_vec(2,k-1)+a1(2)
      endif
    enddo
    base(1)=2.5d0*a3(1)+2d0*a1(1)
    base(2)=2.5d0*a3(2)+2d0*a1(2)
    do i=1,3
      k=k+1
      if(mod(i,2)==1)then
        site_vec(1,k)= base(1) + 0.5d0*dble(i-1)*(a1(1)+a2(1))
        site_vec(2,k)= base(2) + 0.5d0*dble(i-1)*(a1(2)+a2(2))
      else
        site_vec(1,k)= site_vec(1,k-1)+a1(1)
        site_vec(2,k)= site_vec(2,k-1)+a1(2)
      endif
    enddo

    do j = 0,2,2
    do i = 0,3
      k=k+1
      q_list(1,k)=k1(1)/dble(3)*dble(i)+k2(1)/dble(2)*dble(j);
      q_list(2,k)=k1(2)/dble(3)*dble(i)+k2(2)/dble(2)*dble(j);
    enddo
    enddo
!!
!    do j= 1,1,2
!    do i= 1,5,2
!      k=k+1
!      q_list(1,k)=k1(1)/dble(6)*dble(i)+k2(1)/dble(2)*dble(j);
!      q_list(2,k)=k1(2)/dble(6)*dble(i)+k2(2)/dble(2)*dble(j);
!    enddo
!    enddo
!!
!!
    elseif(n==16)then
    ibond=n/2*c_count
    L1=2; L2=2;
    q_num=(L1+1)*(L2+1)
    call allocate_lattice_mem(n,ibond,q_num)
!
    ic=1
    do j=0,L2-1
    do i=0,L1-1
! x  
      ipair(2*ic-1)= (4*L1)*j + 4*i + 0 + 1
      ipair(2*ic  )= (4*L1)*j + 4*i + 1 + 1
      ibond_type(ic)=1
      ic=ic+1
      ipair(2*ic-1)= (4*L1)*j + 4*i             + 2 + 1
      ipair(2*ic  )= (4*L1)*j + 4*(mod(i+1,L1)) + 3 + 1
      ibond_type(ic)=1
      ic=ic+1
! z
      ipair(2*ic-1)= (4*L1)*j + 4*i + 1 + 1
      ipair(2*ic  )= (4*L1)*j + 4*i + 2 + 1
      ibond_type(ic)=3
      ic=ic+1
      ipair(2*ic-1)= (4*L1)*j             + 4*i + 3 + 1
      ipair(2*ic  )= (4*L1)*(mod(j+1,L2)) + 4*i + 0 + 1
      ibond_type(ic)=3
      ic=ic+1
! y
      ipair(2*ic-1)= (4*L1)*j + 4*i + 1 + 1
      ipair(2*ic  )= (4*L1)*j + 4*mod(i+1,L1) + 0 + 1
      ibond_type(ic)=2
      ic=ic+1
      ipair(2*ic-1)= (4*L1)*j + 4*i + 3 + 1
      ipair(2*ic  )= (4*L1)*j + 4*i + 2 + 1
      ibond_type(ic)=2
      ic=ic+1
!!
!!
      base(1)=0.5d0*a3(1) + (a1(1)+a2(1))*dble(i) + a3(1)*dble(3*j)
      base(2)=0.5d0*a3(2) + (a1(2)+a2(2))*dble(i) + a3(2)*dble(3*j)
! 0
      k=4*(j*L1+i)+1
        site_vec(1,k+0)= base(1) 
        site_vec(2,k+0)= base(2) 
! 1
        site_vec(1,k+1)= base(1) + a1(1)
        site_vec(2,k+1)= base(2) + a1(2)
! 2
        site_vec(1,k+2)= base(1) + a1(1) + a3(1)
        site_vec(2,k+2)= base(2) + a1(2) + a3(2)
! 3
        site_vec(1,k+3)= base(1) + a1(1) + a3(1) - a2(1)
        site_vec(2,k+3)= base(2) + a1(2) + a3(2) - a2(2)
    enddo
    enddo
!!
    if(simple.ne.0)then
      if(simple.ne.4)then
!! 3rd N.N.
      do j=0,L2-1
      do i=0,L1-1
! x  
        ipair(2*ic-1)= (4*L1)*j                + 4*i                + 0 + 1
        ipair(2*ic  )= (4*L1)*(mod(j-1+L2,L2)) + 4*(mod(i-1+L1,L1)) + 3 + 1
        ibond_type(ic)=7
        ic=ic+1
        ipair(2*ic-1)= (4*L1)*j                + 4*i                + 2 + 1
        ipair(2*ic  )= (4*L1)*j                + 4*(mod(i-1+L1,L1)) + 1 + 1
        ibond_type(ic)=7
        ic=ic+1
! z
        ipair(2*ic-1)= (4*L1)*j + 4*i + 0 + 1
        ipair(2*ic  )= (4*L1)*j + 4*i + 3 + 1
        ibond_type(ic)=9
        ic=ic+1
        ipair(2*ic-1)= (4*L1)*j                + 4*i + 1 + 1
        ipair(2*ic  )= (4*L1)*(mod(j-1+L2,L2)) + 4*i + 2 + 1
        ibond_type(ic)=9
        ic=ic+1
! y
        ipair(2*ic-1)= (4*L1)*j                + 4*i                + 1 + 1
        ipair(2*ic  )= (4*L1)*(mod(j  +L2,L2)) + 4*(mod(i-1+L1,L1)) + 2 + 1
        ibond_type(ic)=8
        ic=ic+1
        ipair(2*ic-1)= (4*L1)*j                + 4*i                + 3 + 1
        ipair(2*ic  )= (4*L1)*(mod(j+1+L2,L2)) + 4*(mod(i-1+L1,L1)) + 0 + 1
        ibond_type(ic)=8
        ic=ic+1
      enddo
      enddo
      endif

      if(simple==1.or.simple==3.or.simple==4)then
!! 2nd N.N.
        do j=0,L2-1
        do i=0,L1-1
! z
          ipair(2*ic-1)= (4*L1)*j + 4*i             + 0 + 1
          ipair(2*ic  )= (4*L1)*j + 4*(mod(i+1,L1)) + 0 + 1
          ibond_type(ic)=6
          ic=ic+1
          ipair(2*ic-1)= (4*L1)*j + 4*i             + 1 + 1
          ipair(2*ic  )= (4*L1)*j + 4*(mod(i+1,L1)) + 1 + 1
          ibond_type(ic)=6
          ic=ic+1
          ipair(2*ic-1)= (4*L1)*j + 4*i             + 2 + 1
          ipair(2*ic  )= (4*L1)*j + 4*(mod(i+1,L1)) + 2 + 1
          ibond_type(ic)=6
          ic=ic+1
          ipair(2*ic-1)= (4*L1)*j + 4*i             + 3 + 1
          ipair(2*ic  )= (4*L1)*j + 4*(mod(i+1,L1)) + 3 + 1
          ibond_type(ic)=6
          ic=ic+1
        enddo
        enddo
      endif
!!
      if(simple==3.or.simple==4)then
!! 2nd N.N.
        do j=0,L2-1
        do i=0,L1-1
! x
          ipair(2*ic-1)= (4*L1)*j                + 4*i             + 0 + 1
          ipair(2*ic  )= (4*L1)*(mod(j-1+L2,L2)) + 4*i             + 2 + 1
          ibond_type(ic)=4
          ic=ic+1
          ipair(2*ic-1)= (4*L1)*j                + 4*i             + 1 + 1
          ipair(2*ic  )= (4*L1)*(mod(j-1+L2,L2)) + 4*(mod(i+1,L1)) + 3 + 1
          ibond_type(ic)=4
          ic=ic+1
          ipair(2*ic-1)= (4*L1)*j                + 4*i             + 2 + 1
          ipair(2*ic  )= (4*L1)*j                + 4*(mod(i+1,L1)) + 0 + 1
          ibond_type(ic)=4
          ic=ic+1
          ipair(2*ic-1)= (4*L1)*j                + 4*i             + 3 + 1
          ipair(2*ic  )= (4*L1)*j                + 4*i             + 1 + 1
          ibond_type(ic)=4
          ic=ic+1
        enddo
        enddo
!! 2nd N.N.
        do j=0,L2-1
        do i=0,L1-1
! y
          ipair(2*ic-1)= (4*L1)*j                + 4*i                + 0 + 1
          ipair(2*ic  )= (4*L1)*(mod(j-1+L2,L2)) + 4*(mod(i-1+L1,L1)) + 2 + 1
          ibond_type(ic)=5
          ic=ic+1
          ipair(2*ic-1)= (4*L1)*j                + 4*i                + 1 + 1
          ipair(2*ic  )= (4*L1)*(mod(j-1+L2,L2)) + 4*i                + 3 + 1
          ibond_type(ic)=5
          ic=ic+1
          ipair(2*ic-1)= (4*L1)*j                + 4*i                + 2 + 1
          ipair(2*ic  )= (4*L1)*j                + 4*i                + 0 + 1
          ibond_type(ic)=5
          ic=ic+1
          ipair(2*ic-1)= (4*L1)*j                + 4*i                + 3 + 1
          ipair(2*ic  )= (4*L1)*j                + 4*(mod(i-1+L1,L1)) + 1 + 1
          ibond_type(ic)=5
          ic=ic+1
        enddo
        enddo
      endif
    endif
!!
    k=0
    do j= 0,L2
    do i= 0,L1
      k=k+1
      q_list(1,k)=k1(1)/dble(L1)*dble(i)+k2(1)/dble(L2)*dble(j);
      q_list(2,k)=k1(2)/dble(L1)*dble(i)+k2(2)/dble(L2)*dble(j);
      if((i==0.or.i==L1).and.(j==0.or.j==L2))then
        iq=iq+1
        q_list_symmetric(iq)=k
      endif
    enddo
    enddo
!!
    else
      write(*,*) 'Lattice is not supported'
      stop
    endif
!!
end subroutine set_bond_info
!!
subroutine show_bond_list
implicit none
integer :: i,j
    do i=1,ibond
      write(*,'(A,3I8)'),'Bond: ',i,ipair(2*i-1),ipair(2*i)
    enddo

!    do i=1,n
!      write(*,*) i,site_vec(1,i),site_vec(2,i)
!    enddo    

end subroutine show_bond_list

!!
!!
!!
!! Setting of local Hamiltonian
subroutine set_xlocal_Ham
implicit none
double precision, parameter :: eps=.5d-28
complex(KIND(0d0)) :: sx(0:1,0:1),sy(0:1,0:1),sz(0:1,0:1)
complex(KIND(0d0)) :: sxx(0:3,0:3),syy(0:3,0:3),szz(0:3,0:3)
complex(KIND(0d0)) :: sxy(0:3,0:3),syz(0:3,0:3),szx(0:3,0:3)
complex(KIND(0d0)) :: syx(0:3,0:3),szy(0:3,0:3),sxz(0:3,0:3)

integer :: i, j, k, i1, j1, i2, j2
complex(KIND(0d0)) :: x1x,x1y,x1z,x2x,x2y,x2z,x3x,x3y,x3z
double precision :: xh11,xh12,xh13,xh2,xh3,rl
complex(KIND(0d0)) :: pJ11,pJ12,pJ13,pJ2,pJ3,pK11,pK12,pK13,pK2,pK3,pI111,pI112,pI121,pI122,pI131,pI132,pI211,pI212
complex(KIND(0d0)) :: J1x(3,3),J1y(3,3),J1z(3,3),J2x(3,3),J2y(3,3),J2z(3,3),J3z(3,3)

rl=ratio_l

pJ11=dcmplx(xJ1,0d0)
pJ12=dcmplx(xJ1p,0d0)
pJ13=dcmplx(xJ1pp,0d0)
pK11=dcmplx(xK1,0d0)
pK12=dcmplx(xK1p,0d0)
pK13=dcmplx(xK1pp,0d0)

pJ2=dcmplx(xJ2,0d0)
pJ3=dcmplx(xJ3,0d0)
pK2=dcmplx(xK2,0d0)
pK3=dcmplx(xK3,0d0)
pI111=dcmplx(xI11,0d0)
pI112=dcmplx(xI12,0d0)
pI121=dcmplx(xI11p,0d0)
pI122=dcmplx(xI12p,0d0)
pI131=dcmplx(xI11pp,0d0)
pI132=dcmplx(xI12pp,0d0)
pI211=dcmplx(xI21,0d0)
pI212=dcmplx(xI22,0d0)

J1x(:,:)=dcmplx(0d0,0d0)
J1y(:,:)=dcmplx(0d0,0d0)
J1z(:,:)=dcmplx(0d0,0d0)
J2x(:,:)=dcmplx(0d0,0d0)
J2y(:,:)=dcmplx(0d0,0d0)
J2z(:,:)=dcmplx(0d0,0d0)
J3z(:,:)=dcmplx(0d0,0d0)
sx(:,:)=dcmplx(0d0,0d0)
sy(:,:)=dcmplx(0d0,0d0)
sz(:,:)=dcmplx(0d0,0d0)
sx(0,1)=dcmplx(0.5d0,0d0)
sx(1,0)=dcmplx(0.5d0,0d0)
sy(0,1)=dcmplx(0d0,-0.5d0)
sy(1,0)=dcmplx(0d0, 0.5d0)
sz(0,0)=dcmplx( 0.5d0,0d0)
sz(1,1)=dcmplx(-0.5d0,0d0)

J1z(1,1)=pJ11 *rl ; J1z(1,2)=pI111*rl ; J1z(1,3)=pI112*rl;
J1z(2,1)=pI111*rl ; J1z(2,2)=pJ11 *rl ; J1z(2,3)=pI112*rl;
J1z(3,1)=pI112*rl ; J1z(3,2)=pI112*rl ; J1z(3,3)=pK11    ;

J1x(1,1)=pK12     ; J1x(1,2)=pI132*rl ; J1x(1,3)=pI122*rl;
J1x(2,1)=pI132*rl ; J1x(2,2)=pJ13 *rl ; J1x(2,3)=pI121*rl;
J1x(3,1)=pI122*rl ; J1x(3,2)=pI121*rl ; J1x(3,3)=pJ12 *rl;

J1y(1,1)=pJ13 *rl ; J1y(1,2)=pI132*rl ; J1y(1,3)=pI121*rl;
J1y(2,1)=pI132*rl ; J1y(2,2)=pK12     ; J1y(2,3)=pI122*rl;
J1y(3,1)=pI121*rl ; J1y(3,2)=pI122*rl ; J1y(3,3)=pJ12 *rl;

J2x(1,1)=pK2  *rl ; J2x(1,2)=pI211*rl ; J2x(1,3)=pI212*rl;
J2x(2,1)=pI211*rl ; J2x(2,2)=pJ2  *rl ; J2x(2,3)=pI212*rl;
J2x(3,1)=pI212*rl ; J2x(3,2)=pI212*rl ; J2x(3,3)=pJ2  *rl;

J2y(1,1)=pJ2  *rl ; J2y(1,2)=pI211*rl ; J2y(1,3)=pI212*rl;
J2y(2,1)=pI211*rl ; J2y(2,2)=pK2  *rl ; J2y(2,3)=pI212*rl;
J2y(3,1)=pI212*rl ; J2y(3,2)=pI212*rl ; J2y(3,3)=pJ2  *rl;

J2z(1,1)=pJ2  *rl ; J2z(1,2)=pI211*rl ; J2z(1,3)=pI212*rl;
J2z(2,1)=pI211*rl ; J2z(2,2)=pJ2  *rl ; J2z(2,3)=pI212*rl;
J2z(3,1)=pI212*rl ; J2z(3,2)=pI212*rl ; J2z(3,3)=pK2  *rl;

J3z(1,1)=pJ3*rl   ; J3z(2,2)=pJ3*rl   ; J3z(3,3)=pJ3*rl;

do i = 0, 3
do j = 0, 3
        sxx(j,i)=0d0
        syy(j,i)=0d0
        szz(j,i)=0d0
        sxy(j,i)=0d0
        syz(j,i)=0d0
        szx(j,i)=0d0
        syx(j,i)=0d0
        szy(j,i)=0d0
        sxz(j,i)=0d0
        xlocal_Ham(j,i,1)=dcmplx(0d0,0d0)
        xlocal_Ham(j,i,2)=dcmplx(0d0,0d0)
        xlocal_Ham(j,i,3)=dcmplx(0d0,0d0)
        xlocal_Ham(j,i,4)=dcmplx(0d0,0d0)
        xlocal_Ham(j,i,5)=dcmplx(0d0,0d0)
        xlocal_Ham(j,i,6)=dcmplx(0d0,0d0)
        xlocal_Ham(j,i,7)=dcmplx(0d0,0d0)
        xlocal_Ham(j,i,8)=dcmplx(0d0,0d0)
        xlocal_Ham(j,i,9)=dcmplx(0d0,0d0)
enddo
enddo
!!!!
do j2 = 0, 1
do i2 = 0, 1
do j1 = 0, 1
do i1 = 0, 1
  sxx(2*j1+i1,2*j2+i2) = sx(i1,i2)*sx(j1,j2)
  syy(2*j1+i1,2*j2+i2) = sy(i1,i2)*sy(j1,j2)
  szz(2*j1+i1,2*j2+i2) = sz(i1,i2)*sz(j1,j2)
  sxy(2*j1+i1,2*j2+i2) = sx(i1,i2)*sy(j1,j2)
  syz(2*j1+i1,2*j2+i2) = sy(i1,i2)*sz(j1,j2)
  szx(2*j1+i1,2*j2+i2) = sz(i1,i2)*sx(j1,j2)
  syx(2*j1+i1,2*j2+i2) = sy(i1,i2)*sx(j1,j2)
  szy(2*j1+i1,2*j2+i2) = sz(i1,i2)*sy(j1,j2)
  sxz(2*j1+i1,2*j2+i2) = sx(i1,i2)*sz(j1,j2)
enddo
enddo
enddo
enddo
!do j2 = 0, 1
!do i2 = 0, 1
!do j1 = 0, 1
!do i1 = 0, 1
!  write(*,'(2I4,2E16.6)')2*j1+i1,2*j2+i2, szz(2*j1+i1,2*j2+i2) 
!enddo
!enddo
!enddo
!enddo
!!!!
do i = 0, 3
do j = 0, 3
     x1x= J1x(1,1)*sxx(j,i)+J1x(1,2)*sxy(j,i)+J1x(1,3)*sxz(j,i) &
        + J1x(2,1)*syx(j,i)+J1x(2,2)*syy(j,i)+J1x(2,3)*syz(j,i) &
        + J1x(3,1)*szx(j,i)+J1x(3,2)*szy(j,i)+J1x(3,3)*szz(j,i)

     x1y= J1y(1,1)*sxx(j,i)+J1y(1,2)*sxy(j,i)+J1y(1,3)*sxz(j,i) &
        + J1y(2,1)*syx(j,i)+J1y(2,2)*syy(j,i)+J1y(2,3)*syz(j,i) &
        + J1y(3,1)*szx(j,i)+J1y(3,2)*szy(j,i)+J1y(3,3)*szz(j,i)

     x1z= J1z(1,1)*sxx(j,i)+J1z(1,2)*sxy(j,i)+J1z(1,3)*sxz(j,i) &
        + J1z(2,1)*syx(j,i)+J1z(2,2)*syy(j,i)+J1z(2,3)*syz(j,i) &
        + J1z(3,1)*szx(j,i)+J1z(3,2)*szy(j,i)+J1z(3,3)*szz(j,i)

     x2x= J2x(1,1)*sxx(j,i)+J2x(1,2)*sxy(j,i)+J2x(1,3)*sxz(j,i) &
        + J2x(2,1)*syx(j,i)+J2x(2,2)*syy(j,i)+J2x(2,3)*syz(j,i) &
        + J2x(3,1)*szx(j,i)+J2x(3,2)*szy(j,i)+J2x(3,3)*szz(j,i)

     x2y= J2y(1,1)*sxx(j,i)+J2y(1,2)*sxy(j,i)+J2y(1,3)*sxz(j,i) &
        + J2y(2,1)*syx(j,i)+J2y(2,2)*syy(j,i)+J2y(2,3)*syz(j,i) &
        + J2y(3,1)*szx(j,i)+J2y(3,2)*szy(j,i)+J2y(3,3)*szz(j,i)

     x2z= J2z(1,1)*sxx(j,i)+J2z(1,2)*sxy(j,i)+J2z(1,3)*sxz(j,i) &
        + J2z(2,1)*syx(j,i)+J2z(2,2)*syy(j,i)+J2z(2,3)*syz(j,i) &
        + J2z(3,1)*szx(j,i)+J2z(3,2)*szy(j,i)+J2z(3,3)*szz(j,i)

     x3z= J3z(1,1)*sxx(j,i)+J3z(1,2)*sxy(j,i)+J3z(1,3)*sxz(j,i) &
        + J3z(2,1)*syx(j,i)+J3z(2,2)*syy(j,i)+J3z(2,3)*syz(j,i) &
        + J3z(3,1)*szx(j,i)+J3z(3,2)*szy(j,i)+J3z(3,3)*szz(j,i)
!!
     x3x=x3z; x3y=x3z
     if(simple==3)then
       if(simple==3)then
         x2x=x2z; x2y=x2z
       endif
     endif

      xlocal_Ham(j,i,1)=xlocal_Ham(j,i,1)+x1x
      xlocal_Ham(j,i,2)=xlocal_Ham(j,i,2)+x1y
      xlocal_Ham(j,i,3)=xlocal_Ham(j,i,3)+x1z
      xlocal_Ham(j,i,4)=xlocal_Ham(j,i,4)+x2x
      xlocal_Ham(j,i,5)=xlocal_Ham(j,i,5)+x2y
      xlocal_Ham(j,i,6)=xlocal_Ham(j,i,6)+x2z
      xlocal_Ham(j,i,7)=xlocal_Ham(j,i,7)+x3x
      xlocal_Ham(j,i,8)=xlocal_Ham(j,i,8)+x3y
      xlocal_Ham(j,i,9)=xlocal_Ham(j,i,9)+x3z
enddo
enddo
llocal_Ham(:,:,:)=.true.
!! parity check
parity_type=0
do j=1,9
if( abs( dconjg(xlocal_Ham(1,0,j))*xlocal_Ham(1,0,j) )>eps ) parity_type=ibset(parity_type,1)
if( abs( dconjg(xlocal_Ham(2,0,j))*xlocal_Ham(2,0,j) )>eps ) parity_type=ibset(parity_type,1)
if( abs( dconjg(xlocal_Ham(3,0,j))*xlocal_Ham(3,0,j) )>eps ) parity_type=ibset(parity_type,0)
enddo
!!
write(*,*) 'parity :',parity_type
llocal_Ham(:,:,:)=.true.
if(parity_type==0)then
  do j=1,9
    llocal_Ham(0,0,j)=.false.
    llocal_Ham(1,1,j)=.false.
    llocal_Ham(1,2,j)=.false.
    llocal_Ham(2,1,j)=.false.
    llocal_Ham(2,2,j)=.false.
    llocal_Ham(3,3,j)=.false.
  enddo
else
  do j=1,9
  do i=0,3
  do k=0,3
    write(*,'(3I4,2e16.6E3)'),k,i,j,xlocal_Ham(k,i,j)
    if( abs(dconjg(xlocal_Ham(k,i,j))*xlocal_Ham(k,i,j))>eps ) llocal_Ham(k,i,j)=.false.
    if(k==i) llocal_Ham(k,i,j)=.false.
  enddo
  enddo
  enddo
endif
!
end subroutine set_xlocal_Ham
!!
!!
subroutine change_to_TPS_hamiltonian
implicit none
integer :: i,j,k
localHam_max=0d0
do j=1,9
  do i=0,3
  do k=0,3
    if(k==i) localHam_max = localHam_max +  dabs(DREAL(dconjg(xlocal_Ham(k,i,j))*xlocal_Ham(k,i,j)))
  enddo
  enddo
enddo
!!
localHam_max = localHam_max / c_count
xlocal_Ham(:,:,:)=-xlocal_Ham(:,:,:)
do j=1,9
  do i=0,3
  do k=0,3
    if(k==i) xlocal_Ham(k,i,j) =  localHam_max + xlocal_Ham(k,i,j)
  enddo
  enddo
enddo
localHam_max=localHam_max * dble(ibond)
write(*,*) localHam_max
end subroutine change_to_TPS_hamiltonian
!!
!!
subroutine dealloc_mem_2
implicit none
if(allocated(list1))deallocate(list1)
if(allocated(list3))deallocate(list3)
if(allocated(sz_list_p))deallocate(sz_list_p)
if(allocated(sz_list_b))deallocate(sz_list_b)
if(allocated(list_ext3))deallocate(list_ext3)
end subroutine dealloc_mem_2
!!
end module data_mem
