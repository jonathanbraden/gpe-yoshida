#define XIND 1:n

module Equations
  use constants, only : dl, twopi
  use fftw3
  use Simulation

  implicit none
  
!  integer, parameter :: n_terms = 3
  integer, parameter :: n_terms = 5 

  real(dl) :: nu, g_c, g
  real(dl) :: mu
  
contains

  subroutine set_model_parameters(g_,gc_,nu_,mu_)
    real(dl), intent(in) :: g_, gc_, nu_, mu_

    g = g_; g_c = gc_; nu = nu_
    mu = mu_ - 4.*nu_  ! This should probably be tuned depending on the situation
  end subroutine set_model_parameters
   
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
       call evolve_nu_1(this,dt)
    case(5)
       call evolve_nu_2(this,dt)
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

  subroutine evolve_gradient_full(this,dt)
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
  end subroutine evolve_gradient_full

  ! Combine this with the next one, and just put a parameter in to specify the indices
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

  ! Combine this one with the above, and just put in the indices
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

  ! Combined evolution for the above
  subroutine evolve_gradient_single(this,dt,type)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: type

    integer :: fld_ind, grad_ind
    integer :: i_, n

    select case(type)
    case(1)
       fld_ind = 1; grad_ind = 2
    case(2)
       fld_ind = 2; grad_ind = 1
    end select
       
    n = this%nlat
    
    do i_=1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,grad_ind,i_)
       call laplacian_1d_wtype(this%tPair, this%dk)
       this%psi(XIND,fld_ind,i_) = this%psi(XIND,fld_ind,i_) + 0.5_dl*this%tPair%realSpace(XIND)*dt
    enddo
  end subroutine evolve_gradient_single
  
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

    phase_shift = (this%psi(XIND,1,cross_ind)**2 + this%psi(XIND,2,cross_ind)**2)*g_c*dt

    this%tPair%realSpace = this%psi(XIND,2,fld_ind)
    this%psi(XIND,2,fld_ind) = cos(phase_shift)*this%psi(XIND,2,fld_ind) - sin(phase_shift)*this%psi(XIND,1,fld_ind) + nu_loc*this%psi(XIND,1,cross_ind)*dt
    this%psi(XIND,1,fld_ind) = cos(phase_shift)*this%psi(XIND,1,fld_ind) + sin(phase_shift)*this%tPair%realSpace(XIND) - nu_loc*this%psi(XIND,2,cross_ind)*dt
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
    this%psi(XIND,2,fld_ind) = cos(phase_shift)*this%psi(XIND,2,fld_ind) - sin(phase_shift)*this%psi(XIND,1,fld_ind) + nu_loc*this%psi(XIND,1,cross_ind)*dt
    this%psi(XIND,1,fld_ind) = cos(phase_shift)*this%psi(XIND,1,fld_ind) + sin(phase_shift)*this%tPair%realSpace(XIND) - nu_loc*this%psi(XIND,2,cross_ind)*dt
  end subroutine evolve_cross_2
    
    
end module Equations
