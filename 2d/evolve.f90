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
  real(dl) :: dt
  
  call create_lattice(mySim, 128, 16._dl, 1)
  call set_model_parameters(-1.,0.,0.,0.)
  
!  call imprint_sine(mySim,2)
!  call imprint_gray_soliton_pair(mySim,0.1,0.)
  call imprint_bright_soliton(mySim,2.,3.)
  
  call write_lattice_data(mySim,50)

  dt = mySim%dx**2/8._dl
  do i=1,500
     call step_lattice(mySim,dt,10)
     call write_lattice_data(mySim,50)
  enddo
  
contains

  subroutine imprint_sine(this, wave_num)
    type(Lattice), intent(inout) :: this
    integer, intent(in) :: wave_num

    real(dl) :: dx
    integer :: i_, ny

    dx = this%dx
    this%psi = 0._dl
    ny = size(this%yGrid)
    
    do i_=1,ny
       this%psi(:,i_,1,1) = sin(wave_num*twopi*(i_-1)*dx/this%lSize)
       !this%psi(:,i_,2,2) = cos(wave_num*twopi*(i_-1)*dx/this%lSize)
    enddo
  end subroutine imprint_sine

  subroutine imprint_bright_soliton(this,eta,kappa)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: eta, kappa

    integer :: i_, ny
    real(dl) :: x

    this%psi = 0._dl
    ny = size(this%yGrid)
    
    do i_=1,ny
       this%psi(:,i_,1,1) = eta/cosh(eta*this%xGrid)*cos(kappa*this%xGrid)
       this%psi(:,i_,2,1) = eta/cosh(eta*this%xGrid)*sin(kappa*this%xGrid)
    enddo
    
  end subroutine imprint_bright_soliton
  
  subroutine imprint_gray_soliton_pair(this,amp,phi)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: amp, phi

    integer :: i_, ny
    real(dl) :: x0

    x0 = -4._dl
    this%psi = 0._dl

    ny = size(this%yGrid)
    do i_ = 1,ny
       this%psi(:,i_,1,1) = amp*cos(phi) * ( tanh(amp*cos(phi)*(this%xGrid-x0)) - tanh(amp*cos(phi)*(this%xGrid+x0)) )
       this%psi(:,i_,2,1) = amp*sin(phi)
    enddo
  end subroutine imprint_gray_soliton_pair
  
end program Evolve_GPE
