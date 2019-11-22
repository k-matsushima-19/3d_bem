module misc
  !$ use omp_lib
  implicit none
  private
  public :: start_clock, stop_clock, get_memory, assert

  real(8) :: time ! [s]

  interface
     function getmem() result(out) bind(c)
       use iso_c_binding
       integer(c_int) :: out ! [Byte]
     end function getmem
  end interface
  
contains
  subroutine start_clock
    call cpu_time(time)
    !$ time = omp_get_wtime()
  end subroutine start_clock

  subroutine stop_clock(t)
    real(8) :: t
    
    call cpu_time(t)
    !$ t = omp_get_wtime() 

    t = t - time
  end subroutine stop_clock

  real(8) function get_memory()
    get_memory = getmem() / 1024.d0 ! [MB]
  end function get_memory

  subroutine assert(flag, message)
    logical,intent(in) :: flag
    character(*),intent(in),optional :: message

    if(.not. flag) then
       write(*,*) "Assertion failed"
       if(present(message)) write(*,*) message
         
#ifdef IFORT
       call tracebackqq
#endif
       
       call abort
    
    end if
  end subroutine assert
  
end module misc
