! rename to ld1x_utils?

module ld1x_utils
implicit none

private 
public assert

contains

    subroutine assert(condition)
    ! If condition == .false., it aborts the program.
    !
    ! Arguments
    ! ---------
    !
    logical, intent(in) :: condition
    !
    ! Example
    ! -------
    !
    ! call assert(a == 5)

    if (.not. condition) stop "Assert failed in ld1x_utils"
    end subroutine
end module ld1x_utils
