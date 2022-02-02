#include "macros.h"

module Model_Params
  use constants, only : dl
  
  type Model
     real(dl), dimension(:,:), allocatable :: g_vals, g_cross, nu
     real(dl), dimension(:), allocatable :: g_self
     real(dl) :: delta, omega
     integer :: nfld = -1
  end type Model

  real(dl), dimension(:,:), allocatable :: v_trap 
  real(dl), dimension(:,:), allocatable :: g_vals, g_cross, nu
  real(dl), dimension(:), allocatable :: g_self
  real(dl) :: delta, omega
  real(dl) :: mu
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
    
    call free_model_parameters()
    
    allocate(g_self(1:nf))
    allocate(g_vals(1:nf,1:nf))
    allocate(g_cross(1:nf,1:nf))
    allocate(nu(1:nf,1:nf))
    nfld = nf
  end subroutine initialize_model_parameters

  subroutine free_model_parameters()
    if (allocated(g_self)) deallocate(g_self)
    if (allocated(g_vals)) deallocate(g_vals)
    if (allocated(g_cross)) deallocate(g_cross)
    if (allocated(nu)) deallocate(nu)
  end subroutine free_model_parameters
  
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
  
  subroutine set_model_parameters(g_,gc_,nu_,mu_)
    real(dl), intent(in) :: g_, gc_, nu_, mu_

    integer :: i_,j_

    g_self = g_
    g_cross = gc_
    nu = nu_
    mu = mu_
    
    do i_=1,nf
       g_cross(i_,i_) = 0._dl
       nu(i_,i_) = 0._dl
    enddo
  end subroutine set_model_parameters
    
end module model_params
