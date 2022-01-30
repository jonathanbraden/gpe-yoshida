module utils

contains

  integer function newunit(unit)
    integer, intent(out), optional :: unit

    integer, parameter :: umin = 50, umax = 1000
    logical :: o
    integer :: u
    newunit = -1
    do u=umin,umax
       inquire(unit=u,opened=o)
       if (.not.o) then; newunit=u; exit; endif
    enddo
    if (present(unit)) unit=newunit
  end function newunit

end module utils
