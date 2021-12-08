module Equations
  use constants, only : dl, twopi
  use fftw3
  use Simulation

  implicit none
  
  integer, parameter :: n_terms = 3
  
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
       this%tPair%realSpace(1:n) = this%psi(1:n,2,i_)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(1:n,1,i_) = this%psi(1:n,1,i_) - 0.5_dl*this%tPair%realSpace(1:n)*dt
    enddo
  end subroutine evolve_gradient_real

  subroutine evolve_gradient_imag(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_, n

    n = this%nlat
    do i_=1,this%nFld
       this%tPair%realSpace(1:n) = this%psi(1:n,1,i_)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(1:n,2,i_) = this%psi(1:n,2,i_) + 0.5_dl*this%tPair%realSpace(1:n)*dt
    enddo
  end subroutine evolve_gradient_imag
  
  subroutine evolve_potential(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_
    integer :: n
    real(dl) :: g
    real(dl), dimension(1:this%nlat) :: rho
    
    g = 1._dl  ! Adjust model parameter
    n = this%nlat
    
    do i_ = 1,this%nfld
       rho = this%psi(1:n,1,i_)**2 + this%psi(1:n,2,i_)**2
       this%tPair%realSpace = this%psi(1:n,2,i_)
       this%psi(1:n,2,i_) = cos(g*rho*dt)*this%psi(1:n,2,i_) - sin(g*rho*dt)*this%psi(1:n,1,i_)
       this%psi(1:n,1,i_) = cos(g*rho*dt)*this%psi(1:n,1,i_) + sin(g*rho*dt)*this%tPair%realSpace(1:n)
    enddo
  end subroutine evolve_potential

  subroutine evolve_cross_coupling(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

  end subroutine evolve_cross_coupling
  
end module Equations
