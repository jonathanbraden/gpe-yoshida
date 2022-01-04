#define XIND 1:n

module Equations
  use constants, only : dl, twopi
  use fftw3
  use Simulation

  implicit none

!  integer, parameter :: n_terms = 3
!  integer, parameter :: n_terms = 5 
!  integer, parameter :: n_terms = 2
  integer, parameter :: n_terms = 3
  
  real(dl) :: nu, g_c, g
  real(dl) :: delta, omega
  real(dl) :: mu

  real(dl), dimension(1:n_terms) :: t_loc

  real(dl), dimension(:), allocatable :: v_trap
  
contains

  subroutine set_model_parameters(g_,gc_,nu_,mu_)
    real(dl), intent(in) :: g_, gc_, nu_, mu_

    g = g_; g_c = gc_
    nu = nu_
    mu = mu_

    t_loc = 0.
  end subroutine set_model_parameters

  ! A bit wasteful having to pass in xGrid
  ! Currently a temporary fix since the grid isn't stored in the Fourier
  ! transform pairs if it's periodic
  subroutine initialize_trap_potential(this, xGrid, amp)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:this%nLat) :: xGrid
    real(dl), intent(in) :: amp
    
    integer :: i
    
    allocate( v_trap(1:this%nlat) )
    v_trap = -amp*cos( twopi*xGrid/this%lSize )
  end subroutine initialize_trap_potential
  
  subroutine split_equations(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case(1)
       call evolve_gradient_trap_real(this,dt)
    case(2)
       call evolve_potential(this,dt)
    case(3)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations
  
  subroutine split_equations_nlse_w_trap(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case(1)
       call evolve_gradient_trap_real(this,dt)
    case(2)
       call evolve_potential(this,dt)
    case(3)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_nlse_w_trap
  
  ! This isn't tested yet
  subroutine split_equations_gpe(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term
   
    select case (term)
    case (1)
       call evolve_gradient_real(this,dt)
    case (2)
       call evolve_potential(this,dt)
    case (3)
       call evolve_nu_1(this,dt)
    case (4)
       call evolve_nu_2(this,dt)
    case (5)
       call evolve_gradient_imag(this,dt)
    end select
    t_loc(term) = t_loc(term) + dt
  end subroutine split_equations_gpe

  subroutine split_equations_nlse_2_term(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case(1)
       call evolve_gradient_full(this,dt)
    case(2)
       call evolve_potential(this,dt)
    end select
  end subroutine split_equations_nlse_2_term

  subroutine split_equations_nlse_3_term(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case(1)
       call evolve_gradient_real(this,dt)
    case(2)
       call evolve_potential(this,dt)
    case(3)
       call evolve_gradient_imag(this,dt)
    end select
  end subroutine split_equations_nlse_3_term
  
  !>@brief
  !> O(2) symplectic step with operator fusion
  subroutine symp_o2_step(this,dt,w1,w2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt, w1, w2

    integer :: i

    do i=2,n_terms-1
       call split_equations(this, 0.5_dl*w1*dt, i)
    enddo
    call split_equations(this, w1*dt, n_terms)
    do i=n_terms-1,2,-1
       call split_equations(this, 0.5_dl*w1*dt,i)
    enddo
    call split_equations(this,0.5_dl*(w1+w2)*dt,1)
  end subroutine symp_o2_step

  !>@brief
  !> Solve the equation:
  !> \f[
  !>     i\dot{\psi}_i = -\frac{1}{2}\nabla^2\psi_i
  !> \f]
  !> We can alternatively write this as
  !> \f[
  !>    \ddot{R} + \frac{k^4}{4}R = 0
  !>    \ddot{I} + \frac{k^4}{4}I = 0
  !>    \dot{R} = \frac{k^2}{2}I
  !>    \dot{I} = -\frac{k^2}{2}R
  !> \f]
  !> with solutions
  !> \f[
  !>   R(t) = R(0)\cos(\omega t) + I(0)\sin(\omega t)
  !>   I(t) = I(0)\cos(\omega t) - R(0)\sin(\omega t)
  !> \f]
  !> where
  !> \f[
  !>   \omega = \frac{k^2}{2}
  !> \f]
  subroutine evolve_gradient_full(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n, j, nn
    real(dl) :: dk, omega
    complex(dl), dimension(1:this%nlat/2+1) :: fk_real, fk_imag
    
    n = this%nlat; nn = this%nlat/2+1
    dk = this%dk
    
    do i_ = 1,this%nfld
       this%tPair%realSpace(XIND) = this%psi(XIND,1,i_)
       call fftw_execute_dft_r2c(this%tPair%planf, this%tPair%realSpace, this%tPair%specSpace)
       fk_real = this%tPair%specSpace
       this%tPair%realSpace(XIND) = this%psi(XIND,2,i_)
       call fftw_execute_dft_r2c(this%tPair%planf, this%tPair%realSpace, this%tPair%specSpace)
       fk_imag = this%tPair%specSpace
       
       do j=1,nn
          omega = 0.5_dl*(j-1)**2*this%dk**2
          this%tPair%specSpace(j) = cos(omega*dt)*fk_real(j) &
               + fk_imag(j)*sin(omega*dt) 
       enddo
       call fftw_execute_dft_c2r(this%tPair%planb, this%tPair%specSpace, this%tPair%realSpace)
       this%psi(XIND,1,i_) = this%tPair%realSpace(XIND) / dble(n)
       
       do j=1,nn
          omega = 0.5_dl*(j-1)**2*this%dk**2 
          this%tPair%specSpace(j) = cos(omega*dt)*fk_imag(j) & 
               - fk_real(j)*sin(omega*dt) 
       enddo
       call fftw_execute_dft_c2r(this%tPair%planb, this%tPair%specSpace, this%tPair%realSpace)
       this%psi(XIND,2,i_) = this%tPair%realSpace(XIND) / dble(n)
    enddo
  end subroutine evolve_gradient_full

  !>TODO: Write this to do complex rotations
  subroutine evolve_gradient_full_rotation(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n, j, nn
    real(dl) :: dk, omega
    complex(dl), dimension(1:this%nlat/2+1) :: fk_real, fk_imag
    real(dl), dimension(1:this%nlat/2+1) :: amp, phase
    
    n = this%nlat; nn = this%nlat/2+1
    dk = this%dk
    
    do i_ = 1,this%nfld
       this%tPair%realSpace(XIND) = this%psi(XIND,1,i_)
       call fftw_execute_dft_r2c(this%tPair%planf, this%tPair%realSpace, this%tPair%specSpace)
       fk_real = this%tPair%specSpace
       this%tPair%realSpace(XIND) = this%psi(XIND,2,i_)
       call fftw_execute_dft_r2c(this%tPair%planf, this%tPair%realSpace, this%tPair%specSpace)
       fk_imag = this%tPair%specSpace

       do j=1,nn
          omega = 0.5_dl*(j-1)**2*this%dk**2
          
          this%tPair%specSpace(j) = cos(omega*dt)*fk_real(j) &
               + fk_imag(j)*sin(omega*dt) 
       enddo
       call fftw_execute_dft_c2r(this%tPair%planb, this%tPair%specSpace, this%tPair%realSpace)
       this%psi(XIND,1,i_) = this%tPair%realSpace(XIND) / dble(n)
       
       do j=1,nn
          omega = 0.5_dl*(j-1)**2*this%dk**2 
          this%tPair%specSpace(j) = cos(omega*dt)*fk_imag(j) & 
               - fk_real(j)*sin(omega*dt) 
       enddo
       call fftw_execute_dft_c2r(this%tPair%planb, this%tPair%specSpace, this%tPair%realSpace)
       this%psi(XIND,2,i_) = this%tPair%realSpace(XIND) / dble(n)
    enddo
  end subroutine evolve_gradient_full_rotation

  
  !>@brief
  !> Evolve the real part of the fields using the Laplacian term
  !
  !> \f[
  !>    \dot{\psi}_i^{\rm R} = -\frac{1}{2}\nabla^2\psi_i^{\rm I}
  !> ]\f
  subroutine evolve_gradient_real(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: fld_ind, grad_ind
    integer :: i_, n

    fld_ind = 1; grad_ind = 2
    n = this%nlat
    do i_ = 1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,grad_ind,i_)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_) - 0.5_dl*this%tPair%realSpace(XIND)*dt
    enddo
  end subroutine evolve_gradient_real

  !>@brief
  !> Evolve the imaginary part of the fields using the Laplacian term
  !
  !> \f[
  !>    \dot{\psi}_i^{\rm I} = \frac{1}{2}\nabla^2\psi_i^{\rm R}
  !> \f]
  subroutine evolve_gradient_imag(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: fld_ind, grad_ind
    integer :: i_, n

    fld_ind = 2; grad_ind = 1
    n = this%nlat
    do i_=1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,grad_ind,i_)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_) + 0.5_dl*this%tPair%realSpace(XIND)*dt
    enddo
  end subroutine evolve_gradient_imag

  subroutine evolve_gradient_trap_real(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer :: fld_ind, grad_ind
    integer :: i_, n

    fld_ind = 1; grad_ind = 2
    n = this%nlat
    do i_ = 1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,grad_ind,i_)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_) + ( v_trap(XIND)*this%psi(XIND,grad_ind,i_) - 0.5_dl*this%tPair%realSpace(XIND) )*dt
    enddo
  end subroutine evolve_gradient_trap_real
  
  !>@brief
  !> Evolve the imaginary part of the fields using the Laplacian term
  !
  !> \f[
  !>    \dot{\psi}_i^{\rm I} = \frac{1}{2}\nabla^2\psi_i^{\rm R}
  !> \f]
  subroutine evolve_gradient_trap_imag(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: fld_ind, grad_ind
    integer :: i_, n

    fld_ind = 2; grad_ind = 1
    n = this%nlat
    do i_=1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,grad_ind,i_)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_) + (v_trap(XIND)*this%psi(XIND,grad_ind,i_) + 0.5_dl*this%tPair%realSpace(XIND) )*dt
    enddo
  end subroutine evolve_gradient_trap_imag
  
  !>@brief
  !> Combines evolve_gradient_real and evolve_gradient_imag into a single call
  !>
  !> type selects which fields to evolve
  !>  type = 1 evolves the real parts
  !>  type = 2 evolves the imaginary parts
  subroutine evolve_gradient_single(this,dt,type)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: type

    integer :: fld_ind, grad_ind
    integer :: i_, n
    real(dl) :: sig

    select case(type)
    case(1)
       fld_ind = 1; grad_ind = 2
       sig = -1._dl
    case(2)
       fld_ind = 2; grad_ind = 1
       sig = 1._dl
    end select
       
    n = this%nlat
    
    do i_=1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,grad_ind,i_)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_) + 0.5_dl*sig*this%tPair%realSpace(XIND)*dt
    enddo
  end subroutine evolve_gradient_single

  !>@brief
  !> Evolve the 2->2 self-scattering part of the equations
  !>
  !> \f[
  !>    i\dot{\psi}_i = g\left|\psi_i\right|^2\psi_i - \mu\psi_i
  !> \f]
  subroutine evolve_potential(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_
    integer :: n
    real(dl), dimension(1:this%nfld) :: g_loc
    real(dl) :: mu_loc
    real(dl) :: g_cur
    real(dl), dimension(1:this%nlat) :: phase_shift
    
    g_loc = g; mu_loc = mu
    n = this%nlat
    
    do i_ = 1,this%nfld
       g_cur = g_loc(i_)
       phase_shift = this%psi(XIND,1,i_)**2 + this%psi(XIND,2,i_)**2
       phase_shift = (g_cur*phase_shift - mu_loc)*dt
       
       this%tPair%realSpace = this%psi(XIND,2,i_)
       this%psi(XIND,2,i_) = cos(phase_shift)*this%psi(XIND,2,i_) - sin(phase_shift)*this%psi(XIND,1,i_)
       this%psi(XIND,1,i_) = cos(phase_shift)*this%psi(XIND,1,i_) + sin(phase_shift)*this%tPair%realSpace(XIND)
    enddo
  end subroutine evolve_potential

  !>@brief
  !> Evolve the 2->2 self-scattering part of the equations
  !>
  !> \f[
  !>    i\dot{\psi}_i = g\left|\psi_i\right|^2\psi_i - \mu\psi_i
  !> \f]
  !>
  !> The approach of this subroutine is to do a complex rotation
  subroutine evolve_potential_rotation(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_
    integer :: n
    real(dl), dimension(1:this%nfld) :: g_loc
    real(dl) :: mu_loc
    real(dl) :: g_cur
    real(dl), dimension(1:this%nlat) :: phase_shift, rho, theta
    
    g_loc = g; mu_loc = mu
    n = this%nlat
    
    do i_ = 1,this%nfld
       g_cur = g_loc(i_)
       phase_shift = this%psi(XIND,1,i_)**2 + this%psi(XIND,2,i_)**2
       rho = sqrt(phase_shift)
       phase_shift = (g_cur*phase_shift - mu_loc)*dt
       theta = atan2( this%psi(XIND,2,i_),this%psi(XIND,1,i_) ) - phase_shift

       this%psi(XIND,1,i_) = rho*cos(theta)
       this%psi(XIND,2,i_) = rho*sin(theta)
    enddo
  end subroutine evolve_potential_rotation

  ! Combine this into a single subroutine with the next call
  subroutine evolve_nu_1(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n
    real(dl) :: gc_loc, nu_loc
    real(dl), dimension(1:this%nlat) :: rho
    integer :: fld_ind, cross_ind

    n = this%nlat
    fld_ind = 1; cross_ind = 2
    nu_loc = nu
    gc_loc = g_c

    rho = this%psi(XIND,1,fld_ind)**2 + this%psi(XIND,2,fld_ind)**2

    this%psi(XIND,1,fld_ind) = this%psi(XIND,1,fld_ind) - nu_loc*this%psi(XIND,2,cross_ind)*dt
    this%psi(XIND,2,fld_ind) = this%psi(XIND,2,fld_ind) + nu_loc*this%psi(XIND,1,cross_ind)*dt
  end subroutine evolve_nu_1

  subroutine evolve_nu_2(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n
    real(dl) :: gc_loc, nu_loc
    real(dl), dimension(1:this%nlat) :: rho
    integer :: fld_ind, cross_ind

    n = this%nlat
    fld_ind = 2; cross_ind = 1
    nu_loc = nu
    gc_loc = g_c

    this%psi(XIND,1,fld_ind) = this%psi(XIND,1,fld_ind) - nu_loc*this%psi(XIND,2,cross_ind)*dt
    this%psi(XIND,2,fld_ind) = this%psi(XIND,2,fld_ind) + nu_loc*this%psi(XIND,1,cross_ind)*dt
  end subroutine evolve_nu_2

  subroutine evolve_cross_1(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n
    real(dl) :: gc_loc, nu_loc
    real(dl), dimension(1:this%nlat) :: phase_shift
    integer :: fld_ind, cross_ind

    n = this%nlat
    fld_ind = 1; cross_ind = 2

    nu_loc = nu; gc_loc = g_c

    phase_shift = (this%psi(XIND,1,cross_ind)**2 + this%psi(XIND,2,cross_ind)**2)*gc_loc*dt

    this%tPair%realSpace = this%psi(XIND,2,fld_ind)
    this%psi(XIND,2,fld_ind) = cos(phase_shift)*this%psi(XIND,2,fld_ind) &
                             - sin(phase_shift)*this%psi(XIND,1,fld_ind) &
                             + nu_loc*this%psi(XIND,1,cross_ind)*dt
    this%psi(XIND,1,fld_ind) = cos(phase_shift)*this%psi(XIND,1,fld_ind) &
                              + sin(phase_shift)*this%tPair%realSpace(XIND) &
                              - nu_loc*this%psi(XIND,2,cross_ind)*dt
  end subroutine evolve_cross_1

  subroutine evolve_cross_2(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n
    real(dl) :: gc_loc, nu_loc
    real(dl), dimension(1:this%nlat) :: phase_shift
    integer :: fld_ind, cross_ind

    n = this%nlat
    fld_ind = 2; cross_ind = 1

    nu_loc = nu; gc_loc = g_c

    phase_shift = (this%psi(XIND,1,cross_ind)**2 + this%psi(XIND,2,cross_ind)**2)*g_c*dt

    this%tPair%realSpace = this%psi(XIND,2,fld_ind)
    this%psi(XIND,2,fld_ind) = cos(phase_shift)*this%psi(XIND,2,fld_ind) &
                             - sin(phase_shift)*this%psi(XIND,1,fld_ind) &
                             + nu_loc*this%psi(XIND,1,cross_ind)*dt
    this%psi(XIND,1,fld_ind) = cos(phase_shift)*this%psi(XIND,1,fld_ind) &
                             + sin(phase_shift)*this%tPair%realSpace(XIND) &
                             - nu_loc*this%psi(XIND,2,cross_ind)*dt
  end subroutine evolve_cross_2
    
    
end module Equations
