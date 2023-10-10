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

  real(dl) :: dt, dt_grad, dt_m
  integer :: i, nf, u
  type(Lattice) :: mySim
  integer :: ierror
  real(dl) :: error
  real(dl) :: t2, t1
  real(dl) :: nu_, mphi_

  real(dl) :: trap_ratio, grad_ratio, m2_, w_rat
!#ifdef TUNNEL
  nf = 1
  trap_ratio = 2.
  w_rat = 1. ! try various values to set hbar
  !call create_lattice_rectangle(mySim, (/256,256/), (/w_rat*4.*twopi,w_rat*4.*twopi/), nf)
  call create_lattice_rectangle(mySim, (/128,128/), (/32._dl,32._dl/), nf)
  call initialize_model_parameters(nf)
  call set_model_parameters(0.,0.,0.,0.)  ! check

  ! For harmonic trap aligned on axes
  trap_ratio = 2._dl
  call initialize_trap(mySim,(/trap_ratio/),3)
  !call imprint_gaussian_2d( mySim, (/1.,1./sqrt(trap_ratio)/) )
  call imprint_coherent_state( mySim, (/1.,1./sqrt(trap_ratio)/), 3. )
  mySim%v_trap = mySim%v_trap - 0.5_dl*(1.+sqrt(trap_ratio))
  
  ! For harmonic trap not aligned on axes
  !grad_ratio = 2.
  !call initialize_trap(mySim,(/grad_ratio/),2)
  !call imprint_gaussian_2d_diag( mySim, (/1., 1./sqrt(1.+grad_ratio)/) )
  !mySim%v_trap = mySim%v_trap - 0.5_dl*(1.+sqrt(1.+grad_ratio))

  ! For BEC potential
  !trap_ratio = 1.4
  !grad_ratio = 0.2_dl !0.5 ! / 10. for mid
  !m2_ = 1.-1._dl/trap_ratio**2
  !call initialize_trap(mySim,(/trap_ratio,grad_ratio,w_rat/),6)
  !call imprint_gaussian_2d_diag( mySim, (/1./sqrt(m2_), 1./sqrt(m2_+grad_ratio)/) )
  !mySim%v_trap = mySim%v_trap - 0.5_dl*(sqrt(m2_) + sqrt(m2_+grad_ratio))

  !trap_ratio = 100.
  !grad_ratio = 0.
  !m2_ = trap_ratio**2-1._dl
  !call initialize_trap(mySim,(/trap_ratio,grad_ratio/),5)
  !call imprint_gaussian_2d_diag( mySim, (/1./sqrt(m2_), 1./sqrt(m2_+grad_ratio)/) )
  !mySim%v_trap = mySim%v_trap - 0.5_dl*(sqrt(m2_) + sqrt(m2_+grad_ratio))
  
  ! call solve_background_w_grad_flow()  ! If I want automation

  call write_lattice_data(mySim, 50)
  ! call set_chemical_potential() ! Replace this with E_gs and subtract off potential

  dt = 1./512.
  do i=1,200
     call step_lattice(mySim,dt,64)
     call write_lattice_data(mySim,50)
     print*,"step ",i," time = ",mySim%time
  enddo
  
!#endif
  
#ifdef OLD  ! What is this one doing?
  nf = 2
  call create_lattice_rectangle(mySim, (/256,256/), (/32._dl,32._dl/), nf)
!  call create_lattice_rectangle(mySim, (/128,128/), (/ 8._dl, 2._dl /), nf)
  call initialize_model_parameters(nf)
  call set_model_parameters(100.,0.,0.0,0.)
  call initialize_trap(mySim,(/16./),3)
  
!!!!
! Set up initial conditions
!!!!
! Solve for ICs in inhomogeneous bg
  call imprint_gaussian_2d(mySim,(/1._dl,1._dl/))
  call solve_background_w_grad_flow(mySim, 1.e-15, 1.e-15, error, ierror)
  call write_lattice_data(mySim,50)
  call set_chemical_potential(chemical_potential_full(mySim))
   
  dt_grad = minval(mySim%dx)**2/16._dl  ! improve this
  dt_m = twopi/(2.*sqrt(0.01/100.))/dble(16.)
  dt = min(dt_grad,dt_m)
#endif

#ifdef PREHEATING_TRAPPED
!  call setup_preheating_sim_trapped(mySim, 0.01, 0.125*twopi, dt)
  call setup_preheating_sim_trapped(mySim, 0.05, 0.25*twopi, dt)
  ! Now I need to add thw noise on top
  call add_shaped_white_noise(mySim, 1.e-5)
  call write_lattice_data(mySim,50)
  
  call cpu_time(t1)
  do i=1,500
     call step_lattice(mySim,dt,512)
     call write_lattice_data(mySim,50)
     print*,"step ",i," time = ",mySim%time
  enddo
  call cpu_time(t2)

  print*,"runtime is ", (t2-t1)
  print*,"time per step is", (t2-t1)/100./10.
#endif

  ! Uncomment this once I'm done the Schrodinger mucking around
  !call run_preheating_sim_example(0.01, 0.125*twopi, 256)

#ifdef SINGLE_MODE
  ! Add a single-wave
  nu_ = 0.01; mphi_ = 2.*nu_**0.5
  nf = 2
  call create_lattice_rectangle( mySim, (/128,128/), (/50._dl/mphi_, 50._dl/mphi_/), nf )
  call initialize_model_parameters(nf)
  call set_model_parameters(1., 0., nu_, 0.)
  call initialize_trap(mySim, (/1./), 1)
  
  call set_chemical_potential(1._dl-nu_)
  call imprint_phase_wave_2d(mySim, 0.2*twopi, 0._dl, 0.0001, (/3,2/) ) ! sim, phi0, drho, amp, wavenum

  dt = twopi/128./mphi_/4.
  print*,"dt = ",dt
  print*,"dt_grad = ",mySim%dx**2/8.

  call write_lattice_data(mySim,50)
  do i=1,32*16
     call step_lattice(mySim,dt,8*4)
     call write_lattice_data(mySim,50)
  enddo
#endif
  
contains

  ! Need to work out all my unit conversions here
  ! Probably need to relate 2D g to 1D g (via an appropriate integral)
  subroutine setup_preheating_sim_trapped(this, nu_1d, phi0, dt)
    type(Lattice), intent(out) :: this
    real(dl), intent(in) :: nu_1d, phi0
    real(dl), intent(out) :: dt
    
    real(dl) :: error; integer :: ierror
    real(dl) :: geff
    real(dl) :: meff, dt_grad, dt_m
    real(dl) :: len_mass

    ! These should be temporary
    real(dl) :: nu_eff
    real(dl) :: psi_scl
    real(dl) :: len_x, len_y
    integer :: nx, ny
    real(dl) :: lperp, f_norm

    geff = 1.
    len_x = twopi*6/sqrt(nu_1d)/sqrt(1._dl-cos(phi0)) / (geff/sqrt(twopi)) ! add lperp part
    nx = 256
    len_y = 16.; ny = 64
    psi_scl = sqrt(len_x)
    
    ! Work out parameters here
    geff = geff*psi_scl**2
   
    nf = 2
    call create_lattice_rectangle(this, (/nx,ny/), (/len_x, len_y/), nf)
    call initialize_model_parameters(nf)
    call initialize_trap(this, (/32./), 4)

    ! Solve background without nu
    call set_model_parameters(geff, 0., 0., 0.)
    call imprint_transverse_gaussian_2d(this,1._dl) ! Fix scale
    call solve_background_w_grad_flow(this, 1.e-15, 1.e-15, error, ierror)
    ! Renormalize background to account for density.  Fix this later
    this%psi = this%psi * psi_scl 
    call set_chemical_potential(chemical_potential_full(mySim))

    ! Now turn on nu and rotate the condensate
    geff = geff / psi_scl**2
    
    ! Clean up this horrible inlined approach
    f_norm = sum((mySim%psi(1,:,1,1)**2+mySim%psi(1,:,2,1)**2))*this%dx(2)
    print*,"Field norm is ",f_norm
    lperp = sum((mySim%psi(1,:,1,1)**2+mySim%psi(1,:,2,1)**2)**2)*this%dx(2)
    lperp = f_norm / lperp
    print*,"lperp normalize = ",lperp/sqrt(twopi)
    nu_eff = geff*nu_1d / lperp
    print*,"nu_eff = ", nu_eff
    call set_model_parameters(geff, 0., nu_eff, 0.)
    ! Do I want to reevaluate the chemical potential here?
    call rotate_condensate(this,-0.5*phi0, 1)
    call rotate_condensate(this, 0.5*phi0, 2)

    ! Now work out dt (Fix this)
    dt_grad = minval(mySim%dx)**2/16._dl
    dt_m = twopi/(2.*sqrt(nu_1d))/dble(16.)
    print*,"dt_grad = ",dt_grad,", dt_m = ",dt_m
    print*,"geff = ",geff
    print*,"approx period = ", 0.5*twopi * lperp / sqrt(nu_1d)/geff
    dt = min(dt_grad,dt_m)  
  end subroutine setup_preheating_sim_trapped
  
  ! nu =0.01; phi0=0.125*twopi, 16*32 side lengths, n=256 give reasonable results
  ! what's the dt I need there?
  subroutine run_preheating_sim_example(nu, phi0, nLat)
    real(dl), intent(in) :: nu, phi0
    integer, intent(in) :: nLat
    ! Variables to include as input
    real(dl) :: noise_amp
    
    type(Lattice) :: mySim
    integer :: nf
    integer :: nsteps, stepsize
    real(dl) :: m_eff
    integer :: n_samp, n_unstable
    real(dl) :: dt_grad, dt_m, dt_floq, dt
    real(dl) :: kmax, kpeak, lSide
    integer :: fNum
    real(dl) :: alph
    real(dl) :: lyap
    
    ! Variables to turn to input
    noise_amp = 1.e-3 !0.0001 ! Fix this to account for kuv: sigms^2 = k_UV^2/twopi * A w/ A the PS amplitude
    
    fNum = 50
    alph = 8._dl
    
    n_samp = 16    ! Sampling of effective mass frequency
    n_unstable = 8 ! number of unstable modes to resolve
    
    nf = 2
    m_eff = 2._dl*sqrt(nu)
    kmax  = sqrt(0.5*(1.-cos(phi0)))
    kpeak = kmax/2.**0.5  ! approximation for the peak

    lSide = twopi*(n_unstable+0.5)/kpeak
    lSide = lSide / m_eff  ! Convert to program units

    lSide = 16.*32.
    
    call create_lattice_rectangle(mySim, (/nLat,nLat/), (/lSide,lSide/), nf) 
    call initialize_model_parameters(nf)
    call set_model_parameters(1.,0.,nu,0.)
    call initialize_trap(mySim, (/1./), 1)

    call imprint_homogeneous_relative_phase(mySim, phi0, 0.)
    call add_white_noise(mySim,noise_amp)
!    call set_chemical_potential(chemical_potential_full(mySim))
    call set_chemical_potential(1._dl-nu)
    
    dt_grad = minval(mySim%dx)**2/alph  ! improve this
    dt_m = twopi/m_eff/dble(n_samp)
    lyap = phi0**2/8._dl*m_eff  ! check units here.  lyap * T = 1, T = 1/lyap/n_samp
    dt_floq = 1._dl/lyap/dble(n_samp)
    
    dt = min(dt_grad,dt_m)
    nsteps = 2000; stepsize=5

    dt = ((16.*32.)/256.)**2/32.  ! Need to get Floquet exponent
    
    print*,"dt = ",dt
    print*,"dt_out = ", dt*stepsize
    print*,"dt_mass = ", dt_m
    print*,"dt_floq = ", dt_floq
    
    call cpu_time(t1)
    call write_lattice_data(mySim,fNum)
    do i=1,nsteps
       call step_lattice(mySim,dt,stepsize)
       call write_lattice_data(mySim,fNum)
    enddo
    call cpu_time(t2)
  end subroutine run_preheating_sim_example
  
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

    ! What was this for?  Computing norm or something?
    this%tPair%realSpace = this%psi(1:this%nx,1:this%ny,1,1)**2 + this%psi(1:this%nx,1:this%ny,2,1)**2 
  end subroutine imprint_gaussian_2d

  ! Add momentum here
  subroutine imprint_coherent_state(this,sig2,x)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:2), intent(in) :: sig2
    real(dl), intent(in) :: x
    integer :: j,l

    this%psi = 0._dl
    do j = 1,this%ny
       do l=1,this%nfld
          this%psi(1:this%nx,j,1,l) = exp( -0.5_dl*((this%xGrid-x)**2/sig2(1)+this%yGrid(j)**2/sig2(2)) ) / (0.5_dl**2*twopi**2*sig2(1)*sig2(2))**0.25
       enddo
    enddo

    ! What was this for?  Computing norm or something?
    this%tPair%realSpace = this%psi(1:this%nx,1:this%ny,1,1)**2 + this%psi(1:this%nx,1:this%ny,2,1)**2 
  end subroutine imprint_coherent_state
  
  ! Imprint Gaussian initial state in the diagonalized quadratic basis
  ! Need to check all coefficients, etc. in here
  subroutine imprint_gaussian_2d_diag(this,sig2)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:2), intent(in) :: sig2
    real(dl), dimension(1:this%nx,this%nfld,2) :: coord
    real(dl) :: norm
    integer :: j,l

    norm = 1._dl / (0.5_dl**2*twopi**2*sig2(1)*sig2(2))**0.25
    this%psi = 0._dl
    do j = 1,this%ny
       do l=1,this%nfld
          coord(:,l,1) = (this%xGrid(:) + this%yGrid(j))/sqrt(2.)
          coord(:,l,2) = (this%yGrid(j) - this%xGrid(:))/sqrt(2.)
          
          this%psi(1:this%nx,j,1,l) = norm*exp( -0.5_dl*(coord(:,l,1)**2/sig2(1)+coord(:,l,2)**2/sig2(2)) )
       enddo
    enddo
    
  end subroutine imprint_gaussian_2d_diag
  
  subroutine imprint_transverse_gaussian_2d(this,sig2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: sig2
    integer :: j,l

    this%psi = 0._dl
    do l=1,this%nfld
       do j=1,this%ny
          this%psi(1:this%nx,j,1,l) = exp(-0.5_dl*this%yGrid(j)**2/sig2) / (0.5_dl*twopi*sig2)**0.25
       enddo
    enddo
  end subroutine imprint_transverse_gaussian_2d
  
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

  subroutine imprint_homogeneous_relative_density(this, drho, phi_global)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: drho, phi_global

    this%psi = 0._dl
    this%psi(XIND,1,1) = sqrt(1.-drho)*cos(phi_global)
    this%psi(XIND,2,1) = sqrt(1.-drho)*sin(phi_global)
    this%psi(XIND,1,2) = sqrt(1.+drho)*cos(phi_global)
    this%psi(XIND,2,2) = sqrt(1.+drho)*sin(phi_global)
  end subroutine imprint_homogeneous_relative_density

  subroutine imprint_phase_wave_2d(this, phi, drho, amp, wave_num)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: phi, drho, amp
    integer, dimension(1:2), intent(in) :: wave_num

    real(dl), dimension(this%nx,this%ny) :: sine_wave
    integer :: j

    do j=1,this%ny
       sine_wave(:,j) = sin(wave_num(2)*twopi*this%yGrid(j)/this%lSize(2))*sin(wave_num(1)*twopi*this%xGrid(:)/this%lSize(1))
    enddo
    sine_wave = amp*sine_wave

    this%psi(XIND,1,1) = sqrt(1.-drho)*cos(-0.5*(phi+sine_wave))
    this%psi(XIND,2,1) = sqrt(1.-drho)*sin(-0.5*(phi+sine_wave))
    this%psi(XIND,1,2) = sqrt(1.+drho)*cos(0.5*(phi+sine_wave))
    this%psi(XIND,2,2) = sqrt(1.+drho)*sin(0.5*(phi+sine_wave))
  end subroutine imprint_phase_wave_2d
  
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

    call initialize_rand(42,13)

    do l=1,this%nfld
       do i=1,this%ny
          call random_number(phase); call random_number(amp)
          this%psi(:,i,1,l) = this%psi(:,i,1,l) + rms*sqrt(-log(amp))*cos(twopi*phase)
          this%psi(:,i,2,l) = this%psi(:,i,2,l) + rms*sqrt(-log(amp))*sin(twopi*phase)
       enddo
    enddo
  end subroutine add_white_noise

  ! Currently only has fluctuations in x-direction
  subroutine add_shaped_white_noise(this, rms)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: rms

    integer :: i,l
    real(dl), dimension(1:this%nx) :: phase, amp

    do l=1,this%nfld
       call random_number(phase); call random_number(amp)
       do i=1,this%ny
          this%psi(:,i,1,l) = this%psi(:,i,1,l)*( 1.+rms*sqrt(-log(amp))*cos(twopi*phase) )
          this%psi(:,i,2,l) = this%psi(:,i,2,l)*( 1.+rms*sqrt(-log(amp))*sin(twopi*phase) )
       enddo
    enddo
  end subroutine add_shaped_white_noise

  ! Write this to imprint the phase wave
  subroutine imprint_shape_sine_wave(this, wn, amp)
    type(Lattice), intent(inout) :: this
    integer, intent(in) :: wn
    real(dl), intent(in) :: amp
  end subroutine imprint_shape_sine_wave
  
  subroutine initialize_rand(seed, seedfac)
    integer, intent(in) :: seed, seedfac
    integer :: nseed, i
    integer, allocatable, dimension(:) :: seeds

    call random_seed(size=nseed)
    allocate(seeds(1:nseed))
    seeds = seed + seedfac*(/ (i-1, i=1,nseed) /)
    call random_seed(put=seeds)
    deallocate(seeds)
  end subroutine initialize_rand
  
end program Evolve_GPE
