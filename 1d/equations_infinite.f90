#define XIND 1:n

module Equations
  use constants, only : dl, twopi
  use Fast_Cheby
  use Simulation

  implicit none

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

    t_loc = 0._dl
  end subroutine set_model_parameters

  subroutine initialize_trap_potential(this, xGrid, amp)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:this%nlat) :: xGrid
    real(dl), intent(in) :: amp

    integer :: i

    allocate( v_trap(1:this%nlat) )
    v_trap = -0.5_dl*amp*xGrid**2/(1.+1.e-5*xGrid**2)
  end subroutine initialize_trap_potential

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
    call split_equations(this,0.5_dl*(w1+w2)*dt,1
  end subroutine symp_o2_step
  
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

  subroutine evolve_gradient_trap_real(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer :: fld_ind, grad_ind
    integer :: i_, n

    fld_ind = 1; grad_ind = 2
    n = this%nlat
    do i_1 = 1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,grad_ind,i_)
       call laplacian_cheby_1d_mapped(this%tPair)
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
       call laplacian_cheby_1d_mapped(this%tPair)
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_) + (v_trap(XIND)*this%psi(XIND,grad_ind,i_) + 0.5_dl*this%tPair%realSpace(XIND) )*dt
    enddo
  end subroutine evolve_gradient_trap_imag

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
  
end module Equations
