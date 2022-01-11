#define PERIODIC 1
#define XIND 1:this%nlat

module Equations_imag
  use constants, only : dl, twopi
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby
#endif
  use Simulation

  implicit none

  integer, parameter :: n_terms = 3
  real(dl) :: nu, g_c, g
  real(dl) :: mu

  real(dl), dimension(:), allocatable :: v_trap
  
contains

  subroutine set_model_parameters(g_,gc_,nu_,mu_)
    real(dl), intent(in) :: g_, gc_, nu_, mu_

    g = g_; g_c = gc_
    nu = nu_
    mu = mu_
  end subroutine set_model_parameters

  subroutine initialize_trap_potential(this, amp)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: amp
    
    integer :: i
    
    allocate( v_trap(1:this%nlat) )
    v_trap = 0._dl
    !v_trap = -amp*cos( twopi*this%xGrid/this%lSize )
    do i=1,this%nlat
       v_trap = min(amp*0.5*this%xGrid**2,32.)
    enddo

    ! Posch-Teller
!    v_trap = -0.5_dl*amp*(amp+1._dl)/cosh(this%xGrid)**2

    open(unit=99,file='trap.dat')
    do i=1,this%nlat
       write(99,*) this%xGrid(i), v_trap(i)
    enddo
    close(99)
  end subroutine initialize_trap_potential

  subroutine gradient_flow(this,dtau,nsteps)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dtau
    integer, intent(in) :: nsteps

    integer :: i
    do i=1,nsteps
       call gradient_step(this,dtau)
       call renorm_field(this)
    enddo
  end subroutine gradient_flow
  
  subroutine split_equations_bg(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call evolve_heat_equation(this,dt)
    case (2)
       call evolve_nonlinear_interaction(this,dt)
    case (3)
       call renorm_field(this)
    end select
  end subroutine split_equations_bg
    
  subroutine evolve_heat_equation(this, dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i, j

    do i = 1,this%nFld
       this%tPair%realSpace = this%psi(XIND,1,i)
       call laplacian_1d_wtype(this%tPair,this%dk)
       this%psi(XIND,1,i) = this%psi(XIND,1,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,1,i) )*dt

       this%tPair%realSpace = this%psi(XIND,2,i)
       call laplacian_1d_wtype(this%tPair,this%dk)
       this%psi(XIND,2,i) = this%psi(XIND,2,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,2,i) )*dt
    enddo
  end subroutine evolve_heat_equation

  !>@brief
  !> Solve u_t = -\beta |u|^2u
  subroutine evolve_nonlinear_interaction(this, dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    real(dl), dimension(1:this%nlat) :: rho2
    integer :: i
    real(dl) :: g_loc

    g_loc = g
    do i=1,this%nfld
       rho2 = this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2
       rho2 = rho2 / (1._dl + 2._dl*g_loc*rho2*dt)
       this%psi(XIND,1,i) = this%psi(XIND,1,i)*sqrt(rho2)
       this%psi(XIND,2,i) = this%psi(XIND,2,i)*sqrt(rho2)
    enddo
  end subroutine evolve_nonlinear_interaction

  subroutine renorm_field(this)
    type(Lattice), intent(inout) :: this

    real(dl) :: norm
    integer :: i
    
    do i=1,this%nfld
       norm = sum(this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2) * this%dx
       this%psi = this%psi / sqrt(norm)
    enddo
  end subroutine renorm_field

  ! Not debugged
  subroutine gradient_step(this,dtau)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dtau

    integer :: i
    real(dl), dimension(1:this%nlat) :: rho2
    real(dl) :: g_loc

    g_loc = g
    
    do i=1,this%nfld
       rho2 = this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2
       
       this%tPair%realSpace = this%psi(XIND,1,i)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(XIND,1,i) = this%psi(XIND,1,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,1,i) - g_loc*rho2*this%psi(XIND,2,i) ) *dtau

       
       this%tPair%realSpace = this%psi(XIND,2,i)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(XIND,2,i) = this%psi(XIND,2,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,2,i) - g_loc*rho2*this%psi(XIND,2,i) ) * dtau
    enddo
       
  end subroutine gradient_step
  
end module Equations_Imag
