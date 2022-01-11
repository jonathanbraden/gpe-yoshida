#define XIND 1:n,1:n

module Equations
  use constants, only : dl, twopi
  use fftw3
  use Simulation

  implicit none
  
  integer, parameter :: n_terms = 3
  real(dl), dimension(:,:), allocatable :: v_trap
  
contains

  subroutine initialise_fields(this)
    type(Lattice), intent(inout) :: this
  end subroutine initialise_fields

  subroutine initialize_trap(this, amp)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: amp

    integer :: i, nx, ny

    nx = size(this%xGrid); ny = size(this%yGrid)

    allocate(v_trap(1:this%nlat,1:this%nlat))
    
    do i=1,this%nlat
!       v_trap(:,i) = 0._dl
       v_trap(:,i) = min(0.5_dl*(this%xGrid(:)**2+this%yGrid(i)**2),32.)
    enddo

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
       call evolve_gradient_imag(this,dt)
    case(3)
       call evolve_potential(this,dt)
    end select
  end subroutine split_equations

  subroutine evolve_gradient(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n

    n = this%nlat
    do i_ = 1,this%nfld
#ifdef FULL_LIN
       this%tPair%realSpace(1:n) = this%psi(:,1,i_)
       call fftw_execute_dft_r2c(this%tPair%planf, this%tPair%realSpace, this%tPair%specSpace)
       ! Do multiplication
       call fftw_execute_dft_c2r(this%tPair%planf, this%tPair%specSpace, this%tPair%realSpace)
       
       this%tPair%realSpace(1:n) = this%psi(:,2,i_)
       call fftw_execute_dft_r2c(this%tPair%planf, this%tPair%realSpace, this%tPair%specSpace)
       ! Do multiplication
       call fftw_execute_dft_c2r(this%tPair%planf, this%tPair%specSpace, this%tPair%realSpace)
#endif
    enddo
  end subroutine evolve_gradient

  subroutine evolve_gradient_real(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n

    n = this%nlat
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

    integer :: i_, n

    n = this%nlat
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

    integer :: i_, n

    n = this%nlat
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

    integer :: i_, n

    n = this%nlat
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

    integer :: i_
    integer :: n
    real(dl) :: g
    real(dl), dimension(1:this%nlat,1:this%nlat) :: rho
    
    g = 1._dl  ! Adjust model parameter
    n = this%nlat
    
    do i_ = 1,this%nfld
       rho = this%psi(XIND,1,i_)**2 + this%psi(XIND,2,i_)**2
       this%tPair%realSpace = this%psi(XIND,2,i_)
       this%psi(1:n,1:n,2,i_) = cos(g*rho*dt)*this%psi(XIND,2,i_) - sin(g*rho*dt)*this%psi(XIND,1,i_)
       this%psi(XIND,1,i_) = cos(g*rho*dt)*this%psi(XIND,1,i_) + sin(g*rho*dt)*this%tPair%realSpace(XIND)
    enddo
  end subroutine evolve_potential

  subroutine evolve_cross_coupling(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

  end subroutine evolve_cross_coupling
  
end module Equations
