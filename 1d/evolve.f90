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
  call set_model_parameters(1.,0.,0.001,1.)
  
  !call imprint_sine(mySim,2,1._dl)
  call imprint_mean_relative_phase(mySim,0.1_dl*twopi)
  
  call write_lattice_data(mySim,50)
  do i=1,1000
     call step_lattice(mySim,0.00225,80)
     call write_lattice_data(mySim,50)
  enddo
  
contains

  subroutine imprint_sine(this, wave_num, amp)
    type(Lattice), intent(inout) :: this
    integer, intent(in) :: wave_num
    real(dl), intent(in) :: amp

    real(dl) :: dx
    integer :: i_

    dx = this%dx

    do i_=1,this%nlat
       this%psi(i_,1,1) = amp*sin(wave_num*twopi*(i_-1)*dx/this%lSize)
       this%psi(i_,2,2) = amp*cos(wave_num*twopi*(i_-1)*dx/this%lSize)
    enddo
  end subroutine imprint_sine

  subroutine imprint_mean_relative_phase(this,phase)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: phase

    real(dl) :: rho

    rho = 1._dl
    
    this%psi(:,1,1) = rho**0.5
    this%psi(:,2,1) = 0._dl
    
    this%psi(:,1,2) = this%psi(:,1,1)*cos(phase) - this%psi(:,2,1)*sin(phase)
    this%psi(:,2,2) = this%psi(:,2,1)*cos(phase) + this%psi(:,2,1)*sin(phase)
  end subroutine imprint_mean_relative_phase
  
end program Evolve_GPE
