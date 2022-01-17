#define XIND 1:n
#include "macros.h"

module Equations
  use constants, only : dl, twopi
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby
#endif
  use Model_Params
  use Simulation

  implicit none

  integer, parameter :: n_terms = 5

  real(dl), dimension(1:n_terms) :: t_loc
  
contains

  subroutine split_equations(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term
   
    select case (term)
    case (1,-1)
       call evolve_gradient_trap_real(this,dt)
    case (2)
       call evolve_self_scattering(this,dt) ! fix this
    case (3)
       !call evolve_nu_1(this,dt)  ! check this
       call evolve_interspecies_conversion_single(this,dt,1)
    case (4)
       !call evolve_nu_2(this,dt)  ! check this
       !call evolve_interspecies_conversion_single(this,dt,2)
    case (5)
       call evolve_gradient_trap_imag(this,dt)
    end select
    t_loc(term) = t_loc(term) + dt
  end subroutine split_equations
  
  subroutine split_equations_gpe(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term
   
    select case (term)
    case (1,-1)
       call evolve_gradient_trap_real(this,dt)
    case (2)
       call evolve_self_scattering(this,dt) ! fix this
    case (3)
       !call evolve_nu_1(this,dt)  ! check this
       call evolve_interspecies_conversion_single(this,dt,1)
    case (4)
       !call evolve_nu_2(this,dt)  ! check this
       call evolve_interspecies_conversion_single(this,dt,2)
    case (5)
       call evolve_gradient_trap_imag(this,dt)
    end select
    t_loc(term) = t_loc(term) + dt
  end subroutine split_equations_gpe
  
  subroutine split_equations_schrodinger(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call evolve_gradient_trap_real(this,dt)
    case(2)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_schrodinger
  
  subroutine split_equations_nlse_w_trap(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case(1)
       call evolve_gradient_trap_real(this,dt)
    case(2)
       call evolve_self_scattering(this,dt)
    case(3)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_nlse_w_trap
  
  subroutine split_equations_nlse_3_term(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case(1)
       call evolve_gradient_real(this,dt)
    case(2)
       !call evolve_potential(this,dt)
    case(3)
       call evolve_gradient_imag(this,dt)
    end select
  end subroutine split_equations_nlse_3_term

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Building Block Splits of Equations of Motion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
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
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
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
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
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
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_)  &
            + ( v_trap(XIND)*this%psi(XIND,grad_ind,i_) - 0.5_dl*this%tPair%realSpace(XIND) )*dt
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
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_) - (v_trap(XIND)*this%psi(XIND,grad_ind,i_) - 0.5_dl*this%tPair%realSpace(XIND) )*dt
    enddo
  end subroutine evolve_gradient_trap_imag

  subroutine evolve_self_scattering(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    real(dl), dimension(1:this%nlat) :: phase_shift
    real(dl) :: g_cur, mu_loc
    integer :: i_, n

    mu_loc = mu
    n = this%nlat

    phase_shift = 0._dl
    do i_ = 1,this%nfld
       g_cur = g_self(i_)
       phase_shift = this%psi(XIND,1,i_)**2 + this%psi(XIND,2,i_)**2
       phase_shift = (g_cur*phase_shift-mu_loc)*dt
    
       this%tPair%realSpace(:) = this%psi(XIND,2,i_)
       this%psi(XIND,2,i_) = cos(phase_shift)*this%psi(XIND,2,i_) - sin(phase_shift)*this%psi(XIND,1,i_)
       this%psi(XIND,1,i_) = cos(phase_shift)*this%psi(XIND,1,i_) + sin(phase_shift)*this%tPair%realSpace(XIND)
    enddo
  end subroutine evolve_self_scattering

  !>@brief
  !> Evolve intraspecies 2->2 scattering term using complex rotation
  subroutine evolve_self_scattering_rotation(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    real(dl), dimension(1:this%nfld) :: g_cur
    real(dl), dimension(1:this%nfld) :: nu_cur
    integer :: i_,n

    n = this%nlat
    ! write this and test conservation properties
  end subroutine evolve_self_scattering_rotation

  !>@brief
  !> Currently being written.  No interspecies 2->2 scattering yet.
  !>
  !> With the assumptions currently being made, this can't appear as either the first or last term in the splitting scheme
  subroutine evolve_interspecies_interactions_single(this,dt,fld_ind)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: fld_ind

    real(dl), dimension(1:this%nfld) :: g_cur
    real(dl), dimension(1:this%nfld) :: nu_cur
    integer :: l, n

    n = this%nlat
    g_cur = g_cross(:,fld_ind)
    nu_cur = nu(:,fld_ind)

    do l = 1,this%nfld
       this%psi(XIND,1,fld_ind) = this%psi(XIND,1,fld_ind) &
            - nu_cur(l)*dt * this%psi(XIND,1,l)
    enddo
    do l=1,this%nfld
       this%psi(XIND,2,fld_ind) = this%psi(XIND,2,fld_ind) &
            + nu_cur(l)*dt * this%psi(XIND,2,l)
    enddo
  end subroutine evolve_interspecies_interactions_single

  ! Basic checks on the frequency of a homogeneous oscillation done
  subroutine evolve_interspecies_conversion_single(this,dt,fld_ind)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: fld_ind

    real(dl), dimension(1:this%nfld) :: nu_cur
    real(dl), dimension(1:this%nlat,1:2) :: dpsi
    integer :: n, l

    n = this%nlat
    nu_cur = nu(:,fld_ind)

    dpsi = 0._dl
    do l=1,this%nfld
       dpsi(1:n,1) = dpsi(1:n,1) - nu_cur(l) * dt * this%psi(XIND,2,l)
       dpsi(1:n,2) = dpsi(1:n,2) + nu_cur(l) * dt * this%psi(XIND,1,l)
    enddo
    this%psi(1:n,1:2,fld_ind) = this%psi(1:n,1:2,fld_ind) + dpsi  
  end subroutine evolve_interspecies_conversion_single
  
  subroutine evolve_interspecies_scattering_single(this,dt,fld_ind)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: fld_ind

    real(dl), dimension(1:this%nlat) :: phase_rot, temp
    real(dl), dimension(1:this%nfld) :: g_cur
    integer :: l,n

    n = this%nlat
    g_cur = g_cross(:,fld_ind)

    phase_rot = 0._dl
    do l=1,this%nfld
       phase_rot(1:n) = phase_rot(1:n) + g_cur(l)*( this%psi(1:n,1,l)**2 + this%psi(1:n,2,l)**2 )
    enddo
    phase_rot = phase_rot*dt

    temp = this%psi(1:n,1,fld_ind)
    this%psi(1:n,1,fld_ind) = this%psi(1:n,1,fld_ind)*cos(phase_rot) + this%psi(1:n,2,fld_ind)*sin(phase_rot)
    this%psi(1:n,2,fld_ind) = this%psi(1:n,2,fld_ind)*cos(phase_rot) - temp(1:n)*sin(phase_rot)
  end subroutine evolve_interspecies_scattering_single
  
  !>@brief
  !> Evolve full scattering including inter and intra-species for a single species
  subroutine evolve_scattering_single(this,dt,fld_ind)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: fld_ind

    real(dl), dimension(1:this%nlat) :: phase_rot, temp
    real(dl), dimension(1:this%nfld) :: g_cur
    integer :: l,n

    n = this%nlat
    g_cur = g_cross(:,fld_ind)
    g_cur(fld_ind) = g_self(fld_ind)

    phase_rot = 0._dl
    do l=1,this%nfld
       phase_rot(1:n) = phase_rot(1:n) &
            + g_cur(l)*( this%psi(1:n,1,l)**2+this%psi(1:n,2,l)**2 )
    enddo
    phase_rot = phase_rot*dt

    temp = this%psi(1:n,1,fld_ind)
    this%psi(1:n,1,fld_ind) = this%psi(1:n,1,fld_ind)*cos(phase_rot) + this%psi(1:n,2,fld_ind)*sin(phase_rot)
    this%psi(1:n,2,fld_ind) = this%psi(1:n,2,fld_ind)*cos(phase_rot) - temp*cos(phase_rot)
  end subroutine evolve_scattering_single
    
  subroutine evolve_interspecies_scattering_forward(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    real(dl), dimension(1:this%nlat) :: phase_rot, temp
    real(dl), dimension(1:this%nfld) :: g_cur
    integer :: i_,l,n

    n = this%nlat

    do i_=1,this%nfld
       g_cur = g_cross(:,i_)
       phase_rot = 0._dl
       do l=1,this%nfld
          phase_rot(1:n) = phase_rot(1:n) &
               + g_cur(l)*( this%psi(1:n,1,l)**2 + this%psi(1:n,2,l)**2 )

          temp = this%psi(1:n,1,i_)
          this%psi(1:n,1,i_) = cos(phase_rot)*this%psi(1:n,1,i_) + sin(phase_rot)*this%psi(1:n,2,i_)
          this%psi(1:n,2,i_) = cos(phase_rot)*this%psi(1:n,2,i_) - sin(phase_rot)*temp
       enddo
    enddo
  end subroutine evolve_interspecies_scattering_forward

  subroutine evolve_interspecies_scattering_backward(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    real(dl), dimension(1:this%nlat) :: phase_rot, temp
    real(dl), dimension(1:this%nfld) :: g_cur
    integer :: i_,l,n
    integer :: start, end, step
    logical :: forward

    forward = .false.
    if (forward) then
       start = 1; end = this%nfld; step = 1
    else
       start = this%nlat; end = 1; step = -1
    endif
    
    n = this%nlat

    do i_ = this%nfld, 1, -1
       g_cur = g_cross(:,i_)
       phase_rot = 0._dl
       do l=1,this%nfld
          phase_rot(1:n) = phase_rot(1:n) &
               + g_cur(l)*( this%psi(1:n,1,l)**2 + this%psi(1:n,2,l)**2 )

          temp = this%psi(1:n,1,i_)
          this%psi(1:n,1,i_) = cos(phase_rot)*this%psi(1:n,1,i_) + sin(phase_rot)*this%psi(1:n,2,i_)
          this%psi(1:n,2,i_) = cos(phase_rot)*this%psi(1:n,2,i_) - sin(phase_rot)*temp
       enddo
    enddo
  end subroutine evolve_interspecies_scattering_backward
  
  !>@brief
  !> Solve for the nonlinear local dynamics (excluding the potential) using GL10 integrator
  subroutine evolve_local_dynamics(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

  end subroutine evolve_local_dynamics
  
#ifdef OLD
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
#endif
  
end module Equations
