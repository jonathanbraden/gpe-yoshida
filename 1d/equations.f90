#define XIND 1:n

module Equations
  use constants, only : dl, twopi
  use fftw3
  use Simulation

  implicit none
  
!  integer, parameter :: n_terms = 3
  integer, parameter :: n_terms = 5 

  real(dl) :: nu, g_c, g
  
contains

  subroutine initialise_fields(this)
    type(Lattice), intent(inout) :: this
  end subroutine initialise_fields
  
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
    case(4)
       call evolve_cross_1(this,dt)
    case(5)
       call evolve_cross_2(this,dt)
    end select
  end subroutine split_equations

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

  subroutine evolve_gradient(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n

    n = this%nlat
    do i_ = 1,this%nfld
       this%tPair%realSpace(1:n) = this%psi(:,1,i_)
       call fftw_execute_dft_r2c(this%tPair%planf, this%tPair%realSpace, this%tPair%specSpace)
       ! Do multiplication
       call fftw_execute_dft_c2r(this%tPair%planf, this%tPair%specSpace, this%tPair%realSpace)
       
       this%tPair%realSpace(1:n) = this%psi(:,2,i_)
       call fftw_execute_dft_r2c(this%tPair%planf, this%tPair%realSpace, this%tPair%specSpace)
       ! Do multiplication
       call fftw_execute_dft_c2r(this%tPair%planf, this%tPair%specSpace, this%tPair%realSpace)
    enddo
  end subroutine evolve_gradient

  subroutine evolve_gradient_real(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n

    n = this%nlat
    do i_=1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,2,i_)
       call laplacian_1d_wtype(this%tPair, this%dk)
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
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(XIND,2,i_) = this%psi(XIND,2,i_) + 0.5_dl*this%tPair%realSpace(XIND)*dt
    enddo
  end subroutine evolve_gradient_imag
  
  subroutine evolve_potential(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_
    integer :: n
    real(dl), dimension(1:this%nfld) :: g
    real(dl) :: g_cur
    real(dl), dimension(1:this%nlat) :: rho
    
    g = 1._dl  ! Adjust model parameter
    n = this%nlat
    
    do i_ = 1,this%nfld
       g_cur = g(i_)
       rho = this%psi(XIND,1,i_)**2 + this%psi(XIND,2,i_)**2
       this%tPair%realSpace = this%psi(XIND,2,i_)
       this%psi(XIND,2,i_) = cos(g_cur*rho*dt)*this%psi(XIND,2,i_) - sin(g_cur*rho*dt)*this%psi(XIND,1,i_)
       this%psi(XIND,1,i_) = cos(g_cur*rho*dt)*this%psi(XIND,1,i_) + sin(g_cur*rho*dt)*this%tPair%realSpace(XIND)
    enddo
  end subroutine evolve_potential

  subroutine evolve_cross_1(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n
    real(dl) :: g_c, nu
    real(dl), dimension(1:this%nlat) :: rho
    integer :: fld_ind, cross_ind

    n = this%nlat
    fld_ind = 1; cross_ind = 2
    g_c = 0._dl; nu = 0.01

    rho = this%psi(XIND,1,fld_ind)**2 + this%psi(XIND,2,fld_ind)**2

    this%psi(XIND,1,fld_ind) = this%psi(XIND,1,fld_ind) - nu*this%psi(XIND,2,cross_ind)*dt
    this%psi(XIND,2,fld_ind) = this%psi(XIND,2,fld_ind) + nu*this%psi(XIND,1,cross_ind)*dt
  end subroutine evolve_cross_1

  subroutine evolve_cross_2(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n
    real(dl) :: g_c, nu
    real(dl), dimension(1:this%nlat) :: rho
    integer :: fld_ind, cross_ind

    n = this%nlat
    fld_ind = 2; cross_ind = 1
    g_c = 0._dl; nu = 0.01

    rho = this%psi(XIND,1,fld_ind)**2 + this%psi(XIND,2,fld_ind)**2

    this%psi(XIND,1,fld_ind) = this%psi(XIND,1,fld_ind) - nu*this%psi(XIND,2,cross_ind)*dt
    this%psi(XIND,2,fld_ind) = this%psi(XIND,2,fld_ind) + nu*this%psi(XIND,1,cross_ind)*dt
  end subroutine evolve_cross_2
  
end module Equations
