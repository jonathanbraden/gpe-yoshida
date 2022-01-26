program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
!  use gaussianRandomField
  use Yoshida
  use Equations ! remove this for yoshida
  use Equations_Imag
  implicit none

  type(Lattice) :: mySim
  integer :: i
  real(dl) :: dt
  integer :: nf

  nf = 1
  call create_lattice(mySim, 256, 64._dl, 1)
!  call create_lattice(mySim,256,4._dl,nf)
  call initialize_model_parameters(nf)
  call set_model_parameters(10.,0.,0.,0.,nf)
  call initialize_trap(mySim,1.,3)
  
!  call imprint_gray_soliton(mySim,1.,0.)
!  call imprint_black_soliton_pair(mySim,8.)
!  call add_white_noise(mySim,0.02)
!  call imprint_bright_soliton(mySim,2.,3.)
!  call imprint_vortex(mySim,0._dl,0._dl)
  call imprint_gaussian_2d(mySim,1._dl)
  
  call write_lattice_data(mySim,50)
  call solve_background_w_grad_flow(mySim,1.e-15,1.e-15)
  call write_lattice_data(mySim,50)
  
!  dt = minval(mySim%dx)**2/16._dl  
!  do i=1,500
!     call step_lattice(mySim,dt,20)
!     call write_lattice_data(mySim,50)
!  enddo
  
contains

  subroutine imprint_gaussian_2d(this,sig2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: sig2
    integer :: i_

    this%psi = 0._dl
    do i_ = 1,this%ny
       this%psi(:,i_,1,1) = exp(-0.5_dl*(this%xGrid**2+this%yGrid(i_)**2)/sig2)/(0.5_dl*twopi*sig2)**0.5
    enddo

    this%tPair%realSpace = this%psi(1:this%nx,1:this%ny,1,1)**2 + this%psi(1:this%nx,1:this%ny,2,1)**2
  end subroutine imprint_gaussian_2d
  
  subroutine imprint_sine(this, wave_num)
    type(Lattice), intent(inout) :: this
    integer, intent(in) :: wave_num

    real(dl) :: dx
    integer :: i_, ny

    dx = this%dx(1)
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

  subroutine imprint_gray_soliton(this, amp, phi)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: amp, phi

    integer :: i_
    real(dl) :: x0

    x0 = 0._dl
    this%psi = 0._dl

    do i_ = 1,this%ny
       this%psi(:,i_,1,1) = amp*cos(phi)*tanh(amp*cos(phi)*(this%xGrid-x0))
       this%psi(:,i_,2,1) = amp*sin(phi)
    enddo
  end subroutine imprint_gray_soliton
  
  subroutine imprint_black_soliton_pair(this,x0)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: x0

    integer :: i_, ny

    this%psi = 0._dl
    
    ny = size(this%yGrid)
    do i_ = 1,ny
       this%psi(:,i_,1,1) = ( tanh((this%xGrid-x0)) - tanh((this%xGrid+x0)) + 1._dl )
       this%psi(:,i_,2,1) = 0._dl
    enddo
  end subroutine imprint_black_soliton_pair

  !>@brief
  !> Imprint approximate vortex with winding number 1
  !
  ! \psi = r/sqrt(2+r**2)e^{i\phi} for wind = 1
  !      = (x+iy)/sqrt(2+r**2)
  !
  ! To Do : Need to add factors for different condensate normalization
  !         and different g normalizations
  subroutine imprint_vortex(this,x0,y0)
    type(Lattice), intent(inout) :: this
    real(dl) :: x0, y0
    
    integer :: i,j, l
    real(dl) :: rho2, y, x

    do l=1,this%nfld
       do j=1,this%ny
          y = this%yGrid(j) - y0
          do i=1,this%nx
             x = this%xGrid(i) - x0
             rho2 = x**2 + y**2
             this%psi(i,j,1,l) = x/sqrt(2._dl+rho2)
             this%psi(i,j,2,l) = y/sqrt(2._dl+rho2)
          enddo
       enddo
    enddo
    
  end subroutine imprint_vortex
  
  !>@brief
  !> Add white noise to the condensate fields
  !> The RMS fluctuation amplitude of each field is set by rms.  The RMS of each real part is thus rms/sqrt(2)
  subroutine add_white_noise(this,rms)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: rms

    integer :: i,l
    real(dl), dimension(1:this%nx) :: phase, amp

!    call initialize_rand(42,13)

    do l=1,this%nfld
       do i=1,this%ny
          call random_number(phase); call random_number(amp)
          this%psi(:,i,1,l) = this%psi(:,i,1,l) + rms*sqrt(-log(amp))*cos(twopi*phase)
          this%psi(:,i,2,l) = this%psi(:,i,2,l) + rms*sqrt(-log(amp))*sin(twopi*phase)
       enddo
    enddo
  end subroutine add_white_noise
  
end program Evolve_GPE
