#include "macros.h"
#define XIND 1:this%nlat

module Model_Params
  use constants, only : dl
  use Simulation
  
  type Model
     real(dl), dimension(:,:), allocatable :: g_cross, nu
     real(dl), dimension(:), allocatable :: g_self
     real(dl) :: delta, omega
     integer :: nfld
     integer :: type
  end type Model

  real(dl), dimension(:), allocatable :: v_trap
  real(dl), dimension(:,:), allocatable :: g_cross, nu
  real(dl), dimension(:), allocatable :: g_self
  real(dl) :: mu
  real(dl) :: delta, omega
  integer :: nfld
  integer :: type

contains

  subroutine set_chemical_potential(mu_)
    real(dl), intent(in) :: mu_

    mu = mu_
  end subroutine set_chemical_potential
  
  subroutine set_model_parameters(g_,gc_,nu_,mu_,nf)
    real(dl), intent(in) :: g_, gc_, nu_, mu_
    integer, intent(in) :: nf

    integer :: i_
    
    if (allocated(g_self))  deallocate(g_self)
    if (allocated(g_cross)) deallocate(g_cross)
    if (allocated(nu))      deallocate(nu)
    
    allocate(g_self(1:nf))
    allocate( g_cross(1:nf,1:nf), nu(1:nf,1:nf) )

    g_self = g_
    g_cross = gc_
    nu = nu_
    mu = mu_
    
    do i_=1,nf
       g_cross(i_,i_) = 0._dl
       nu(i_,i_) = 0._dl
    enddo
    
    t_loc = 0.
  end subroutine set_model_parameters

  subroutine initialize_trap_potential(this, amp, type)
    type(Lattice), intent(inout) :: this
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

  ! Fix this to compute integral properly for cheby calculation
  ! This is missing the contributions from mu and cross-coupling g
  real(dl) function chemical_potential(this) result(mu)
    type(Lattice), intent(inout) :: this
    
    real(dl), dimension(1:this%nlat) :: rho2, mu_loc
    integer :: i
    real(dl) :: g_loc

    mu = 0._dl
    mu_loc = 0._dl
    
    do i=1,this%nfld
       g_loc = g_self(i)
       rho2 = this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2
       
       this%tPair%realSpace = this%psi(XIND,1,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,1,i)

       this%tPair%realSpace = this%psi(XIND,2,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,2,i)

       mu_loc = mu_loc + v_trap*rho2 + g_loc*rho2**2
    enddo
#if defined(PERIODIC)
    mu = this%dx*sum(mu_loc)
#elif defined(INFINITE)
    mu = sum(mu_loc*this%tPair%quad_weights)
#endif
  end function chemical_potential

  ! Fix this to compute integral properly for cheby calculation
  ! Debug this to make sure it's correct
  real(dl) function chemical_potential_full(this) result(mu)
    type(Lattice), intent(inout) :: this
    
    real(dl), dimension(1:this%nlat,1:this%nfld) :: rho
    real(dl), dimension(1:this%nlat) :: mu_loc
    real(dl), dimension(1:this%nfld) :: g_loc, nu_loc
    integer :: i,l

    mu = 0._dl
    mu_loc = 0._dl

    rho = this%psi(XIND,1,1:this%nfld)**2 + this%psi(XIND,2,1:this%nfld)**2
    do i=1,this%nfld
       g_loc = g_cross(:,i); g_loc(i) = g_self(i)
       nu_loc = nu(:,i)
       
       this%tPair%realSpace = this%psi(XIND,1,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,1,i)

       this%tPair%realSpace = this%psi(XIND,2,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,2,i)

       mu_loc = mu_loc + v_trap*rho(:,i)
       do l=1,this%nfld
          mu_loc = mu_loc + g_loc(l)*rho(:,i)*rho(:,l) &
               - nu_loc(l)*( this%psi(XIND,1,i)*this%psi(XIND,1,l) + this%psi(XIND,2,i)*this%psi(XIND,2,l) )
       enddo
    enddo
#if defined(PERIODIC)
    mu = this%dx*sum(mu_loc)
#elif defined(INFINITE)
    mu = sum(mu_loc*this%tPair%quad_weights)
#endif
  end function chemical_potential_full

  real(dl) function field_norm(this) result(norm)
    type(Lattice), intent(inout) :: this
    integer :: l
    real(dl) :: norm_loc

    norm = 0._dl
    do l=1,this%nfld
#if defined(PERIODIC)
       norm_loc = this%dx*sum(this%psi(1:this%nlat,1:2,1:this%nfld)**2)
#elif defined(INFINITE)
       norm_loc = sum( (this%psi(1:this%nlat,1,l)**2 + this%psi(1:this%nlat,2,l)**2)*this%tPair%quad_weights )
#endif
       norm = norm + norm_loc
    enddo
  end function field_norm
  
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
