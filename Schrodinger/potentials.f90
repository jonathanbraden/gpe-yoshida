module potentials
  use constants, only : dl, twopi
  implicit none
  
contains

  real(dl) function potential_free(x,y,params) result(pot)
    real(dl), intent(in) :: x,y
    real(dl), dimension(:), intent(in) :: params

    pot = 0._dl
  end function potential_free

  real(dl) function potential_quadratic(x,y) result(pot)
    real(dl), intent(in) :: x,y

    pot = 0.5_dl*(x**2+y**2)
  end function potential_quadratic
  
end module potentials
