module Initial_Conditions_Two_Field
  use constants, only : dl, twopi
  use Simulation, only : Lattice
  implicit none

contains
  
  ! Debug all of these things
  ! In particular, check real vs. imaginary part of fields, etc
  subroutine initialise_homogeneous_fv(this,rho_bg,nu,dg,gc)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: rho_bg
    real(dl), intent(in) :: nu, dg, gc
    real(dl) :: sig

    sig = -1._dl  ! Chooses false vacuum
    theta = find_max_angle(nu,dg,gc)
    this%psi(XIND,1,1) = (2._dl*rho_bg)**0.5*cos(theta); this%psi(XIND,2,1) = 0._dl ! This is field 1
    this%psi(XIND,1,2) = -(2._dl*rho_bg)**0.5*sin(theta); this%psi(XIND,2,2) = 0._dl ! This is field 2
  end subroutine initialise_homogeneous_fv
  
  real(dl) function find_max_angle(nu,dg,gc) result(theta)
    real(dl), intent(in) :: nu, dg, gc
    real(dl) :: nu_bg, dg_bg, dtheta
    integer :: l

    integer, parameter :: maxit = 16
    real(dl), parameter :: eps = 1.e-16

    theta = 0.125_dl*twopi + 0.25_dl*dg/(1._dl-gc-nu)
    nu_bg = nu / (1._dl-gc); dg_bg = dg / (1._dl-gc)

    do l=1,maxit
       dtheta = -dvTheta(theta,nu_bg,dg_bg) / d2vTheta(theta,nu_bg,dg_bg)
       theta = theta + dtheta
       if (abs(dtheta) < eps) exit
    enddo
    if (l==maxit) print*,"Newton method failed to find a background"
  end function find_max_angle

  real(dl) function dvTheta(theta, nu_bg, dg_bg) result(dv)
    real(dl), intent(in) :: theta, nu_bg, dg_bg
    real(dl), parameter :: sig = -1._dl  ! Choose true or false vacuum

    dv = cos(2._dl*theta)*sin(2._dl*theta) + sig*nu_bg*cos(2._dl*theta) + 0.5_dl*dg_bg*sin(2._dl*theta)
  end function dvTheta

  real(dl) function d2vTheta(theta, nu_bg, dg_bg) result(d2v)
    real(dl), intent(in) :: theta, nu_bg, dg_bg
    real(dl), parameter :: sig = -1._dl  ! Choose true or false vacuum
    d2v = 2._dl*(cos(2._dl*theta)**2 - sin(2._dl*theta)**2) - 2._dl*sig*nu_bg*sin(2._dl*theta) + dg_bg*cos(2._dl*theta)
  end function d2vTheta

  
end module Initial_Conditions_Two_Field
