module Equation_Splits
  use constants, only : dl, twopi
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby
#endif
  use Simultion
  
contains

  subroutine split_equations_nlse_w_trap(this, dt, term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call evolve_gradient_trap_real(this, dt)
    case (2)
       call evolve_potential(this, dt)
    case (3)
       call evolve_gradient_trap_imag(this, dt)
    end select
  end subroutine split_equations_nlse_w_trap

  subroutine split_equations_nlse_2_part(this, dt, term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call evolve_gradient_full(this, dt)
    case (2)
       call evolve_potential(this, dt)
    end select
  end subroutine split_equations_nlse_2_part

  subroutine split_equations_nlse_3_part(this, dt, term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call evolve_gradient_real(this, dt)
    case (2)
       call evolve_potential(this, dt)
    case (3)
       call evolve_gradient_imag(this, dt)
    end select
  end subroutine split_equations_nlse_3_part

  subroutine split_equation_schrodinger(this, dt, term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call evolve_gradient_trap_real(this, dt)
    case (2)
       call evolve_gradient_trap_real(this, dt)
    end select
  end subroutine split_equation_schrodinger
  
end module Equation_Splits
