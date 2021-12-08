program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use gaussianRandomField
  use Yoshida
  use Equations ! remove this for yoshida
  implicit none

  type(Lattice) :: mySim
  integer :: i
  
  call create_lattice(mySim, 16, 4._dl, 2)

  call imprint_sine(mySim,2)
  
  call write_lattice_data(mySim,50)
  do i=1,1000
     call step_lattice(mySim,0.00225,80)
     call write_lattice_data(mySim,50)
  enddo
  
contains

  subroutine imprint_sine(this, wave_num)
    type(Lattice), intent(inout) :: this
    integer, intent(in) :: wave_num

    real(dl) :: dx
    integer :: i_

    dx = this%dx

    do i_=1,this%nlat
       this%psi(i_,1,1) = sin(wave_num*twopi*(i_-1)*dx/this%lSize)
       this%psi(i_,2,2) = cos(wave_num*twopi*(i_-1)*dx/this%lSize)
    enddo
  end subroutine imprint_sine
  
end program Evolve_GPE
