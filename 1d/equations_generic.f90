module Equations_Generic
  use constants, only : dl
  
contains

  subroutine split_equations_example(this,dt,term,num_terms)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term
    integer, intent(out), optional :: n_terms

    integer, parameter :: n_terms = 3
    
    if (present(num_terms)) then
       num_terms = n_terms
       return
    endif

    select case (term)
    case(1,-1)
       call evolve_gradient_trap_real(this,dt)
    case(2,-2)
       call evolve_self_scattering(this,dt)
    case(3,-3)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_example
  
end module Equations_Generic
