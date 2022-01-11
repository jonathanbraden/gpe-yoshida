#define XIND 1:n
#define INFINITE 1
!#define BOUNDARIES 1

! TO DO : Need to fix the boundary conditions to avoid exponential growth

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

  subroutine initialize_trap_potential(this, xGrid, amp, type)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:this%nlat) :: xGrid
    real(dl), intent(in) :: amp
    integer, intent(in), optional :: type

    integer :: i
    integer :: type_

    type_ = 1; if (present(type)) type_ = type
    allocate( v_trap(1:this%nlat) )

    select case (type_)
    case (1)
       v_trap = 0._dl
    case(2)
       v_trap = -0.5_dl*amp*(amp+1._dl) / cosh(this%xGrid)**2
    case(3)
       do i=1,this%nlat
          v_trap(i) = min(0.5_dl*this%xGrid(i)**2,32.)
       enddo
    end select
       
    open(unit=99,file='trap.dat')
    do i=1,this%nlat
       write(99,*) this%xGrid(i), v_trap(i)
    enddo
    close(99)
  end subroutine initialize_trap_potential
  
  subroutine split_equations(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case(1)
       call evolve_gradient_trap_real(this,dt)
    case(2)
       call evolve_scattering_term(this,dt)
    case(3)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations

  ! Debug this, then use it to study 1,2, and 3 particle Schrodinger
  ! See how decay rate changes with N
  subroutine split_equations_schrodinger(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case(1)
       call evolve_gradient_trap_real(this,dt)
    case(2)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_schrodinger
  
  subroutine evolve_gradient_trap_real(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer :: fld_ind, grad_ind
    integer :: i_, n

    fld_ind = 1; grad_ind = 2
    n = this%nlat
    do i_ = 1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,grad_ind,i_)
       call laplacian_cheby_1d_mapped(this%tPair)
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_)   &
            + ( v_trap(XIND)*this%psi(XIND,grad_ind,i_)        &
            - 0.5_dl*this%tPair%realSpace(XIND) )*dt
    enddo
    
#ifdef BOUNDARIES
    this%psi(1,:,:) = 0._dl; this%psi(this%nlat,:,:) = 0._dl
#endif
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
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_)    &
            - ( v_trap(XIND)*this%psi(XIND,grad_ind,i_)         &
            - 0.5_dl*this%tPair%realSpace(XIND) )*dt
    enddo
    
#ifdef BOUNDARIES
    this%psi(1,:,:) = 0._dl; this%psi(this%nlat,:,:) = 0._dl
#endif
  end subroutine evolve_gradient_trap_imag

  !>@brief
  !> Evolve the 2->2 self-scattering part of the equations
  !>
  !> \f[
  !>    i\dot{\psi}_i = g\left|\psi_i\right|^2\psi_i - \mu\psi_i
  !> \f]
  subroutine evolve_scattering_term(this,dt)
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

#ifdef BOUNDARIES
    this%psi(1,:,:) = 0._dl; this%psi(this%nlat,:,:) = 0._dl
#endif
  end subroutine evolve_scattering_term
  
end module Equations
