module Fluctuations
  use constants, only : dl, twopi
  use gaussianRandomField

  implicit none
  
  type SpecParams
     real(dl) :: dk, rho
     integer :: nCut
     integer :: num_atoms
     character(8) :: type='BOGO'
  end type SpecParams
  
contains

  function make_spec_params(rho, len, nu, lamEff, type, modes, nCut, fv_) result(params)
    real(dl), intent(in) :: rho, len, nu, lamEff
    character(*), intent(in) :: type
    logical, dimension(1:2), intent(in) :: modes
    integer, intent(in) :: nCut
    logical, intent(in), optional :: fv_  ! Don't make optional, encode somewhere else

    type(SpecParams) :: params

    ! Finish this
  end function make_spec_params

  subroutine spectrum_bogoliubov(spec, params, is_relative)
    real(dl), dimension(:), intent(out) :: spec
    type(SpecParams), intent(in) :: params
    logical, intent(in) :: is_relative

    

    
end module Fluctuations
