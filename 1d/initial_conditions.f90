module Initial_Conditions
  use constants, only : dl, twopi
  use Simulation, only : Lattice
  implicit none

contains
  
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
  
end module Initial_Conditions
