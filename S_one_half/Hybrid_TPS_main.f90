!!
program main
use timer
use parameters
use parallel_lib
use data_mem
use working_area
use parallelized_cg_method
use parallelized_lanczos_method
use check_results
use parallelized_diag, only : mk_list12, alloc_listt
use DynamicStractureFactor
use parallelized_TPS_method
use def_output
use openmp_lib
use rand_num1
implicit none

integer :: itr, list_max
double precision :: ene(max_itr), pm
integer :: num_threads
logical :: restart=.false.
logical :: restart2=.true.
logical :: fread=.true.
integer :: pme, pnp, i_color, o_np, q1, q2, lsize, lrange, r_seed
logical :: w_parallel, l_stat
real :: sim_time
character(LEN=2) :: typ
!!
!! initialization of environment
call set_openmp_env
call get_thread_num(num_threads)
!!
call get_cline(d1=sim_time,d2=lsize,d3=lrange,d4=o_np,d5=q1,d6=q2,d7=pm,d8=typ,d9=r_seed)
!!
call my_mpi_start()
!!
call init_timer
call timer_start(10)
call timer_start(1)
!!
call check_system_memsize
call bcast_sim_data(num_threads,o_np,sim_time,lsize,lrange,q1,q2,pm,typ,r_seed)
call set_comp_time(sim_time)
call set_lattice(lsize,lrange,pm)
call set_qrange(q1,q2,typ)
pme=get_me()
pnp=get_np()
!!
call set_itr_counter(w_parallel,i_color,o_np,pnp,pme,l_stat)
call my_mpi_split(i_color,w_parallel)
call myset_num_threads(num_threads)
!!
call set_dump_name
call int_ran(19*pme+r_seed)
!!
call set_parallel_bits
call set_xlocal_Ham
call set_bond_info
!!
call barrier
!if(pme==0) call show_bond_list
!call barrier
!!
call make_cpu_list2
!!
call check_interval_dump(progress)
call timer_stop(1)
write(*,*) typ
write(*,*) progress
!!
!***********************************************************
if((.not.progress(1)).and.(.not.progress(3)).and.(.not.progress(7))) fread=.false.
if((.not.progress(2))) restart2=.false.
!***********************************************************
!!
if(dsf_comp=='zz'.or.dsf_comp=='ee'.or.dsf_comp=='tp')then
  call mk_list12(total_sz,lme)
  call clear_mem
  call prep_gsvec(max_ldim)
!!
!!
!!
!!*** Eigenvalues 
  if(progress(1).or.progress(3))then
    call barrier
    call timer_start(2)
    if(progress(1))then
      call parallelized_lnc1(nvec,itr,ene,xvec,restart)
      if(lme==0)then
 100  format(/' [Eigenvalues]  '/2x,8f14.8/' [Iteration number]'/i8)
        write(*,100)ene(1:8),itr
      endif
    endif
    call timer_stop(2)

!*** Ground-state eigenvector
    call barrier
    if(progress(2))then
      call parallelized_lncv1(1,ene,xvec,itr,restart)
      call check1(xvec)
    endif
call timer_start(3)
    if(progress(3)) call parallelized_inv1(ene(1),xvec,restart2)
call timer_stop(3)
!!
!!
  elseif(progress(7))then
    call change_to_TPS_hamiltonian
    call init_res_dat(lme,total_sz)
    call barrier
    call timer_start(2)
    call parallelized_TPS(nvec,itr,ene,restart)
    call timer_stop(2)
  endif
!!
!!
!!*** Precision check and correlation functions
  call readwrite_vec_data(lme,fread,max_ldim,ene(1),xvec)
  call check1(xvec)
  write(*,*) 'GS is checked.'
elseif(dsf_comp=='xx')then
  fread=.false.
  call mk_list12(total_sz,lme)
  call clear_mem
  call alloc_listt(total_sz,lme,tdim)
  call prep_gsvec(max_tdim)
  call readwrite_vec_data(lme,fread,max_tdim,ene(1),xvec)
!!
!!
else
  write(*,*) 'DSF comp error.'
  stop
endif
!!
!!
!! DSF calculation
call barrier
call timer_start(4)
if(lme==0) write(*,*) progress(:)
if(progress(5)) call calc_DSF(lme,ene(1),xvec,i_color)
if(progress(6)) call calc_extended_DSF(lme,ene(1),xvec,i_color)
call timer_stop(4)
!!
!!
call dealloc_mem_TPS
call dealloc_mem
call dealloc_mem_2
call dealloc_mem_parallel
call timer_stop(10)
call make_timer_log(pme)
call mpi_end
end
