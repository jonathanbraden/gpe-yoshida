#include "macros.h"

module Model_Params
  use constants, only : dl
  use Simulation
  
  type Model
     real(dl), dimension(:,:), allocatable :: g_vals, g_cross, nu
     real(dl), dimension(:), allocatable :: g_self
     real(dl) :: delta, omega
     integer :: nfld = -1
  end type Model

  real(dl), dimension(:,:), allocatable :: v_trap 
  real(dl), dimension(:,:), allocatable :: g_vals, g_cross, nu
  real(dl), dimension(:), allocatable :: g_self
  real(dl) :: mu
  real(dl) :: delta, omega
  integer :: nfld
  integer :: type

contains

  subroutine create_model_parameters(this,nf)
    type(Model), intent(out) :: this
    integer, intent(in) :: nf

    call destroy_model_parameters(this)
    this%nfld = nf
    allocate(this%g_vals(1:nf,1:nf), this%nu(1:nf,1:nf), this%g_cross(1:nf,1:nf))
    allocate(this%g_self(1:nf))
  end subroutine create_model_parameters

  subroutine destroy_model_parameters(this)
    type(Model), intent(inout) :: this

    if (allocated(g_vals)) deallocate(g_vals)
    if (allocated(g_cross)) deallocate(g_cross)
    if (allocated(nu)) deallocate(nu)
    if (allocated(g_self)) deallocate(g_self)
    this%nfld = -1
  end subroutine destroy_model_parameters
  
  subroutine set_chemical_potential(mu_)
    real(dl), intent(in) :: mu_

    mu = mu_
  end subroutine set_chemical_potential

  subroutine initialize_model_parameters(nf)
    integer, intent(in) :: nf

    if (allocated(g_self)) deallocate(g_self)
    if (allocated(g_vals)) deallocate(g_vals)
    if (allocated(g_cross)) deallocate(g_cross)
    if (allocated(nu)) deallocate(nu)

    allocate(g_self(1:nf))
    allocate(g_vals(1:nf,1:nf))
    allocate(g_cross(1:nf,1:nf))
    allocate(nu(1:nf,1:nf))
    nfld = nf
  end subroutine initialize_model_parameters
  
  subroutine set_model_parameters_full(g_new, nu_new)
    real(dl), intent(in), dimension(1:nfld,1:nfld) :: g_new, nu_new

    integer :: i

    g_vals = g_new
    g_cross = g_vals
    nu = nu_new
    
    do i=1,nfld
       g_self(i) = g_vals(i,i)
       g_cross(i,i) = 0._dl
       nu(i,i) = 0._dl
    enddo
  end subroutine set_model_parameters_full
  
  subroutine set_model_parameters(g_,gc_,nu_,mu_,nf)
    real(dl), intent(in) :: g_, gc_, nu_, mu_
    integer, intent(in) :: nf

    integer :: i_,j_

    call create_model_parameters(nf)
    
    g_self = g_
    g_cross = gc_
    nu = nu_
    mu = mu_
    
    do i_=1,nf
       g_cross(i_,i_) = 0._dl
       nu(i_,i_) = 0._dl
    enddo
  end subroutine set_model_parameters

  ! Fix this to compute integral properly for cheby calculation
  real(dl) function chemical_potential(this) result(mu)
    type(Lattice), intent(inout) :: this
    
    real(dl), dimension(1:this%nx,1:this%ny) :: rho2, mu_loc
    integer :: l
    real(dl) :: g_loc

    mu = 0._dl
    mu_loc = 0._dl
    
    do l=1,this%nfld
       g_loc = g_self(i)
       rho2 = this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2
       
       this%tPair%realSpace = this%psi(XIND,1,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,1,l)

       this%tPair%realSpace = this%psi(XIND,2,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,2,l)

       mu_loc = mu_loc + v_trap*rho2 + g_loc*rho2**2
    enddo
#if defined(PERIODIC)
    mu = this%dx*sum(mu_loc)
#elif defined(INFINITE)
    mu = sum(mu_loc*this%tPair%quad_weights)
#endif
  end function chemical_potential

  real(dl) function energy(this) result(en)
    type(Lattice), intent(inout) :: this

    real(dl), dimension(1:this%nlat) :: rho2, en_loc
    integer :: i, l
    real(dl) :: g_loc

    en = 0._dl; en_loc = 0._dl

    do i=1,this%nfld
       g_loc = g_self(i)
       rho2 = this%psi(XIND,1,1)**2 + this%psi(XIND,2,i)**2

       do l=1,2
          this%tPair%realSpace = this%psi(XIND,l,i)
#if defined(PERIODIC)
          call laplacian_1d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
          call laplacian_cheby_1d_mapped(this%tPair)
#endif
          en_loc = en_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,l,i)
       enddo
       en_loc = en_loc + v_trap*rho2 + 0.5_dl*g_loc*rho2**2
    enddo
#if defined(PERIODIC)
    en = this%dx*sum(en_loc)
#elif defined(INFINITE)
    en = sum(en_loc*this%tPair%quad_weights)
#endif
  end function energy
    
end module model_params
