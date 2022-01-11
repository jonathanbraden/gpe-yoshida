#define XIND 1:nx,1:ny
#define PERIODIC

module Equations
  use constants, only : dl, twopi
  use fftw3
  use Simulation

  implicit none

  real(dl) :: g, g_c, nu, mu
  integer, parameter :: n_terms = 3
  real(dl), dimension(:,:), allocatable :: v_trap
  
contains

  subroutine set_model_parameters(g_,gc_,nu_,mu_)
    real(dl), intent(in) :: g_, gc_, nu_, mu_

    g = g_; g_c = gc_
    nu = nu_
    mu = mu_
  end subroutine set_model_parameters
  
  subroutine initialize_trap(this, amp)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: amp

    integer :: i, nx, ny

    nx = size(this%xGrid); ny = size(this%yGrid)

    allocate(v_trap(1:nx,1:ny))

    v_trap = 0._dl
    
!    do i=1,ny
!       v_trap(:,i) = 0._dl
!       v_trap(:,i) = min(0.5_dl*(this%xGrid(:)**2+this%yGrid(i)**2),32.)
!    enddo

    open(unit=99,file='trap.bin',access='stream')
    write(99) v_trap
    close(99)
    open(unit=99,file='xgrid.dat')
    do i=1,nx
       write(99,*) this%xGrid(i)
    enddo
    close(99)
    open(unit=99,file='ygrid.dat')
    do i=1,ny
       write(99,*) this%yGrid(i)
    enddo
    close(99)
  end subroutine initialize_trap
  
  subroutine split_equations(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call evolve_gradient_real(this,dt)
    case(2)
       call evolve_potential(this,dt)
    case(3)
       call evolve_gradient_imag(this,dt)
    end select
  end subroutine split_equations

  subroutine evolve_gradient_real(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, nx, ny

    nx = this%nx; ny = this%ny
    do i_=1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,2,i_)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       this%psi(XIND,1,i_) = this%psi(XIND,1,i_) - 0.5_dl*this%tPair%realSpace(XIND)*dt
    enddo
  end subroutine evolve_gradient_real

  subroutine evolve_gradient_imag(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, nx, ny

    nx = this%nx; ny = this%ny
    do i_=1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,1,i_)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       this%psi(XIND,2,i_) = this%psi(XIND,2,i_) + 0.5_dl*this%tPair%realSpace(XIND)*dt
    enddo
  end subroutine evolve_gradient_imag

  subroutine evolve_gradient_trap_real(this, dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, nx, ny

    nx = this%nx; ny = this%ny
    do i_=1, this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,1,i_)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       this%psi(XIND,2,i_) = this%psi(XIND,2,i_)   &
            + ( 0.5_dl*this%tPair%realSpace(XIND) - v_trap(XIND) )*dt
    enddo
  end subroutine evolve_gradient_trap_real

  subroutine evolve_gradient_trap_imag(this, dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, nx, ny

    nx = this%nx; ny = this%ny
    do i_=1, this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,1,i_)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       this%tPair%realSpace(XIND) = this%psi(XIND,1,i_) &
            - ( 0.5_dl*this%tPair%realSpace(XIND) - v_trap(XIND) )*dt
    enddo
  end subroutine evolve_gradient_trap_imag
    
  subroutine evolve_potential(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, nx, ny
    real(dl) :: g_loc, mu_loc
    real(dl), dimension(1:this%nx,1:this%ny) :: rho

    g_loc = g; mu_loc = mu
    
    nx = this%nx; ny = this%ny
    
    do i_ = 1,this%nfld
       rho = this%psi(XIND,1,i_)**2 + this%psi(XIND,2,i_)**2
       this%tPair%realSpace = this%psi(XIND,2,i_)
       this%psi(1:nx,1:ny,2,i_) = cos(g_loc*rho*dt)*this%psi(XIND,2,i_) - sin(g_loc*rho*dt)*this%psi(XIND,1,i_)
       this%psi(XIND,1,i_) = cos(g_loc*rho*dt)*this%psi(XIND,1,i_) + sin(g_loc*rho*dt)*this%tPair%realSpace(XIND)
    enddo
  end subroutine evolve_potential

  subroutine evolve_cross_coupling(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

  end subroutine evolve_cross_coupling
  
end module Equations
