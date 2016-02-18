module openmp_lib
!$ use omp_lib
implicit none
private
integer :: threads
public set_openmp_env, myset_num_threads, get_thread_num

contains
subroutine set_openmp_env
!$ use omp_lib
implicit none
threads=1
!$ threads=omp_get_max_threads()
!$ call omp_set_num_threads(threads)
end subroutine set_openmp_env

subroutine get_thread_num(mm)
implicit none
integer, intent(inout) :: mm
mm = threads
end subroutine get_thread_num

subroutine myset_num_threads(mm)
!$ use omp_lib
implicit none
integer, intent(in) :: mm
!$ call omp_set_num_threads(mm)
threads = mm
write(*,*) 'setting num threads :',threads
end subroutine myset_num_threads
end module openmp_lib
