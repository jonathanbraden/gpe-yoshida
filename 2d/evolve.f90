#include "macros.h"

program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
!  use gaussianRandomField
  use Yoshida
  use Equations ! remove this for yoshida
  use Equations_Imag
  implicit none

  real(dl) :: dt
  integer :: i, nf, u
  type(Lattice) :: mySim
  integer :: ierror
  real(dl) :: error
  real(dl) :: t2, t1
  
  nf = 1
  call create_lattice_rectangle(mySim, (/256,128/), (/64._dl,16._dl/), nf)
!  call create_lattice_rectangle(mySim, (/256,256/), (/16.*32._dl,16.*32._dl/), nf) ! preheating
  call initialize_model_parameters(nf)
  call set_model_parameters(0.,0.,0.,0.)
  call initialize_trap(mySim,(/16./),3)

!!!!
! Set up initial conditions
!!!!
! Solve for ICs in inhomogeneous bg
  call imprint_gaussian_2d(mySim,(/1._dl,0.25_dl/))
  call write_lattice_data(mySim,50)
  
!  call imprint_homogeneous_relative_phase(mySim,0.125*twopi,0.)

!  call set_chemical_potential(chemical_potential_full(mySim))
!  call set_chemical_potential(0.99)
!  call rotate_condensate(mySim,0.1_dl,2)
!  call add_white_noise(mySim,0.01)
  
  open(unit=newunit(u),file='chem_pot.dat')
  write(u,*) "# Chemical Potential for NLSE in asymmetric harmonic trap"
  write(u,*) "# V = 1/2*(x^2+16*y^2)"
  write(u,*) "# g/(hbar*omega)  mu   energy"
  write(u,*) "# Reference:  Gaussian initial state"
  write(u,*) "# mu = ", chemical_potential_full(mySim)," E = ", energy(mySim)
  do i=0,50
     call set_model_parameters(5.*i,0.,0.,0.)
     call solve_background_w_grad_flow(mySim,1.e-15,1.e-15,error,ierror)
     call write_lattice_data(mySim,50)
     write(u,*) 5.*i, chemical_potential_full(mySim), energy(mySim)
  enddo
  
  ! Add some calculation of timescales
  ! m2eff = 2\sqrt{nu} ! Back ground oscillations
  ! om_k ~ dx^2
  ! Use dispersion relationship to set dt
  
  dt = minval(mySim%dx)**2/32._dl
  call cpu_time(t1)
  do i=1,100
!     call step_lattice(mySim,dt,10)
!     call write_lattice_data(mySim,50)
  enddo
  call cpu_time(t2)

  print*,"runtime is ", (t2-t1)
  print*,"time per step i s", (t2-t1)/500./20.
  
contains

  subroutine rotate_condensate(this,dphi,ind)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dphi
    integer, intent(in) :: ind

    real(dl) :: amp, phi
    integer :: i,j

    do j=1,this%ny
       do i=1,this%nx
          amp = sqrt(this%psi(i,j,1,ind)**2 + this%psi(i,j,2,ind)**2)
          phi = atan2(this%psi(i,j,2,ind),this%psi(i,j,1,ind)) + dphi
          this%psi(i,j,1,ind) = amp*cos(phi) 
          this%psi(i,j,2,ind) = amp*sin(phi)
       enddo
    enddo
  end subroutine rotate_condensate

  subroutine global_rotation(this,dphi)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dphi

    real(dl) :: amp, phi
    integer :: i,j,l

    do l=1,this%nfld
       do j=1,this%ny
          do i=1,this%nx
             amp = sqrt(this%psi(i,j,1,l)**2 + this%psi(i,j,2,l)**2)
             phi = atan2(this%psi(i,j,2,l),this%psi(i,j,1,l)) + dphi
             this%psi(i,j,1,l) = amp*cos(phi)
             this%psi(i,j,2,l) = amp*sin(phi)
          enddo
       enddo
    enddo
  end subroutine global_rotation

  subroutine imprint_gaussian_2d(this,sig2)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:2), intent(in) :: sig2
    integer :: j,l

    this%psi = 0._dl
    do j = 1,this%ny
       do l=1,this%nfld
          this%psi(1:this%nx,j,1,l) = exp( -0.5_dl*(this%xGrid**2/sig2(1)+this%yGrid(j)**2/sig2(2)) ) / (0.5_dl**2*twopi**2*sig2(1)*sig2(2))**0.25
       enddo
    enddo

    this%tPair%realSpace = this%psi(1:this%nx,1:this%ny,1,1)**2 + this%psi(1:this%nx,1:this%ny,2,1)**2
  end subroutine imprint_gaussian_2d

  subroutine imprint_homogeneous_relative_phase(this,phi,phi_global)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: phi, phi_global
    real(dl) :: amp

    amp = 1._dl
    
    this%psi = 0._dl
    this%psi(XIND,1,1) = amp*cos(phi_global)
    this%psi(XIND,2,1) = amp*sin(phi_global)
    this%psi(XIND,1,2) = cos(phi+phi_global)*amp
    this%psi(XIND,2,2) = sin(phi+phi_global)*amp
  end subroutine imprint_homogeneous_relative_phase
  
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
