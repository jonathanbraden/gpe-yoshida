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
  
  !  call create_lattice(mySim, 16, 4._dl, 2)
  call create_lattice(mySim,512,200._dl,2)
  call set_model_parameters(-1.,0.,0.,0.)
  
  call imprint_bright_soliton(mySim,0.9,1.)
  !call imprint_sine(mySim,2,1._dl)
  !call imprint_mean_relative_phase(mySim,0.1_dl*twopi)
  
  call write_lattice_data(mySim,50)
  do i=1,100
     call step_lattice(mySim,0.02,10)
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

    this%psi = 0._dl
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

  subroutine imprint_bright_soliton(this,eta,kappa)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: eta, kappa

    integer :: i_; real(dl) :: x
    
    this%psi = 0._dl

    do i_=1,this%nlat
       x = (i_-1)*this%dx - 0.5*this%lSize
       this%psi(i_,1,1) = eta/cosh(eta*x)*cos(kappa*x)
       this%psi(i_,2,1) = eta/cosh(eta*x)*sin(kappa*x)
    enddo
  end subroutine imprint_bright_soliton
  
end program Evolve_GPE
