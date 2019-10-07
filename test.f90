


module mod_param


  implicit none
  public
  

  integer rank

end module



module mod_memory
implicit none

public 
  real(8), dimension(:), allocatable :: array

contains

  subroutine initMem()
    use mod_param
    implicit none
    integer i
    allocate(array(0:rank))

    do i=0,rank 
      array(i) = 1.0*i
    enddo
  end

end module




program main

use mod_param
use mod_memory
implicit none

integer i


rank = 134

write(*,*) rank

call something()

call initMem()

do i=0,rank 
  write(*,*) array(i)
enddo


end program








subroutine something()
use mod_param
implicit none

write(*,*) rank

end 


