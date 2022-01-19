#include "macros.h"

module Equations
  use constants, only : dl, twopi
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby_2D
#endif
  use Simulation

  implicit none

  real(dl), dimension(:,:), allocatable :: g_cross, nu
  real(dl), dimension(:), allocatable :: g_self
  
  real(dl) :: g, g_c, mu
  integer, parameter :: n_terms = 3
  real(dl), dimension(:,:), allocatable :: v_trap
  
contains

  subroutine set_chemical_potential(mu_)
    real(dl), intent(in) :: mu_

    mu = mu_
  end subroutine set_chemical_potential
  
  subroutine set_model_parameters(g_,gc_,nu_,mu_,nf)
    real(dl), intent(in) :: g_, gc_, nu_, mu_
    integer, intent(in) :: nf

    integer :: i_
    
    if (allocated(g_self)) deallocate(g_self)
    if (allocated(g_cross)) deallocate(g_cross)
    if (allocated(nu)) deallocate(nu)

    allocate( g_self(1:nf) )
    allocate( g_cross(1:nf,1:nf), nu(1:nf,1:nf) )

    g_self = g_
    g_cross = gc_
    nu = nu_
    mu = mu_

    do i_ = 1,nf
       g_cross(i_,i_) = 0._dl
       nu(i_,i_) = 0._dl
    enddo
  end subroutine set_model_parameters
  
  subroutine initialize_trap(this, amp)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: amp

    integer :: i, nx, ny

    nx = size(this%xGrid); ny = size(this%yGrid)

    allocate(v_trap(1:this%nx,1:this%ny))

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
       call evolve_self_scattering(this,dt)
    case(3)
       call evolve_gradient_imag(this,dt)
    end select
  end subroutine split_equations

  subroutine split_equations_schrodinger(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1,-1)
       call evolve_gradient_trap_real(this,dt)
    case(2,-2)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_schrodinger
  
  subroutine evolve_gradient_real(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_

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

    integer :: i_

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

    integer :: i_

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

    integer :: i_

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
    
  subroutine evolve_self_scattering(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i_
    real(dl) :: g_loc, mu_loc
    real(dl), dimension(XIND) :: phase_shift

    g_loc = g; mu_loc = mu
    
    do i_ = 1,this%nfld
       phase_shift = this%psi(XIND,1,i_)**2 + this%psi(XIND,2,i_)**2
       phase_shift = (g_loc*phase_shift - mu)*dt

       this%tPair%realSpace = this%psi(XIND,2,i_)
       this%psi(XIND,2,i_) = cos(phase_shift)*this%psi(XIND,2,i_) - sin(phase_shift)*this%psi(XIND,1,i_)
       this%psi(XIND,1,i_) = cos(phase_shift)*this%psi(XIND,1,i_) + sin(phase_shift)*this%tPair%realSpace(XIND)
    enddo
  end subroutine evolve_self_scattering

  subroutine evolve_interspecies_conversion_single(this,dt,fld_ind)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: fld_ind

    real(dl), dimension(1:this%nx,1:this%ny,1:2) :: dpsi
    real(dl), dimension(1:this%nfld) :: nu_cur
    integer :: l

    nu_cur = nu(:,fld_ind)

    dpsi = 0._dl
    do l=1,this%nfld
       dpsi(XIND,1) = dpsi(XIND,1) - nu_cur(l) * dt * this%psi(XIND,2,l)
       dpsi(XIND,2) = dpsi(XIND,2) + nu_cur(l) * dt * this%psi(XIND,1,l)
    enddo
    this%psi(XIND,1:2,fld_ind) = this%psi(XIND,1:2,fld_ind) + dpsi
  end subroutine evolve_interspecies_conversion_single

  subroutine evolve_interspecies_conversion(this,dt,fwd)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    logical, intent(in) :: fwd

    real(dl), dimension(1:this%nx,1:this%ny,1:2) :: dpsi
    real(dl), dimension(1:this%nfld) :: nu_cur
    integer :: m,l
    integer :: st, en, step

    if (fwd) then
       st = 1; en = this%nfld; step = 1
    else
       st = this%nfld; en = 1; step = -1
    endif
    
    do m = st,en,step
       nu_cur = nu(:,m)
       dpsi = 0._dl
       do l=1,this%nfld
          dpsi(XIND,1) = dpsi(XIND,1) - nu_cur(l) * dt * this%psi(XIND,2,l)
          dpsi(XIND,2) = dpsi(XIND,2) + nu_cur(l) * dt * this%psi(XIND,1,l)
       enddo
       this%psi(XIND,1:2,m) = this%psi(XIND,1:2,m) + dpsi
    enddo
  end subroutine evolve_interspecies_conversion

  subroutine evolve_interspecies_scattering_single(this,dt,fld_ind)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: fld_ind

    real(dl), dimension(1:this%nx,1:this%ny) :: phase_rot, temp
    real(dl), dimension(1:this%nfld) :: g_cur
    integer :: l
    
    g_cur = g_cross(:,fld_ind); g_cur(fld_ind) = 0._dl
    phase_rot = 0._dl
    do l = 1,this%nfld
       phase_rot(1:this%nx,1:this%ny) = phase_rot(1:this%nx,1:this%ny) + g_cur(l)*( this%psi(XIND,1,l)**2 + this%psi(XIND,2,l) )**2
    enddo
    phase_rot = phase_rot * dt

    temp = this%psi(XIND,1,fld_ind)
    this%psi(XIND,1,fld_ind) = this%psi(XIND,1,fld_ind)*cos(phase_rot) + this%psi(XIND,2,fld_ind)*sin(phase_rot)
    this%psi(XIND,2,fld_ind) = this%psi(XIND,2,fld_ind)*cos(phase_rot) - temp(XIND)*sin(phase_rot)
  end subroutine evolve_interspecies_scattering_single

  subroutine evolve_interspecies_scattering(this,dt,fwd)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    logical, intent(in) :: fwd

    real(dl), dimension(1:this%nx,1:this%ny) :: phase_rot, temp
    real(dl), dimension(1:this%nfld) :: g_cur
    integer :: st, en, step
    integer :: l, m
    
    if (fwd) then
       st = 1; en = this%nfld; step = 1
    else
       st = this%nfld; en = 1; step = -1
    endif

    do m = st, en, step
       g_cur = g_cross(:,m); g_cur(m) = 0._dl
       phase_rot = 0._dl
       do l = 1,this%nfld
          phase_rot(XIND) = phase_rot(XIND) + g_cur(l)*( this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2 )
       enddo
       phase_rot = phase_rot * dt

       temp = this%psi(XIND,1,m)
       this%psi(XIND,1,m) = this%psi(XIND,1,m)*cos(phase_rot) + this%psi(XIND,2,m)*sin(phase_rot)
       this%psi(XIND,2,m) = this%psi(XIND,2,m)*cos(phase_rot) - temp(XIND)*sin(phase_rot)
    enddo
  end subroutine evolve_interspecies_scattering
  
end module Equations
