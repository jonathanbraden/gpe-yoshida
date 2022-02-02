#include "macros.h"

module Equations
  use constants, only : dl, twopi
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby_2D
#endif
  use Model_Params
  use Simulation
  
  implicit none

  integer, parameter :: n_terms = 5
 
contains
    
  subroutine initialize_trap(this, params, type)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: params
    integer, intent(in) :: type

    integer :: i, j 

    if (allocated(this%v_trap)) deallocate(this%v_trap)
    allocate(this%v_trap(1:this%nx,1:this%ny))

    select case (type)
    case (1)
       this%v_trap = 0._dl
    case(2)
       this%v_trap = 0._dl
    case(3) ! Harmonic trap in both directions
       do j=1,this%ny
          this%v_trap(:,j) = min(0.5_dl*(this%xGrid(:)**2+params(1)*this%yGrid(j)**2),32.)
       enddo
    case(4)  ! Trap only the y-direction
       do j=1,this%ny
          this%v_trap(:,j) = min(0.5_dl*params(1)*this%yGrid(j)**2,32.)
       enddo
    case default
       this%v_trap = 0._dl
    end select

    open(unit=99,file='trap.bin',access='stream')
    write(99) this%v_trap
    close(99)
    open(unit=99,file='xgrid.dat')
    do i=1,this%nx
       write(99,*) this%xGrid(i)
    enddo
    close(99)
    open(unit=99,file='ygrid.dat')
    do i=1,this%ny
       write(99,*) this%yGrid(i)
    enddo
    close(99)
  end subroutine initialize_trap

  subroutine split_equations(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (abs(term))
    case (1)
       call evolve_gradient_trap_real(this,dt)
    case (2)
       call evolve_self_scattering(this,dt)
    case(3)
       call evolve_interspecies_conversion_single(this,dt,1)
    case(4)
       call evolve_interspecies_conversion_single(this,dt,2)
    case(5)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations

  subroutine split_equations_schrodinger(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (abs(term))
    case (1)
       call evolve_gradient_trap_real(this,dt)
    case (2)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_schrodinger

  subroutine split_equations_nlse(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (abs(term))
    case (1)
       call evolve_gradient_real(this,dt)
    case (2)
       call evolve_self_scattering(this,dt)
    case (3)
       call evolve_gradient_imag(this,dt)
    end select
  end subroutine split_equations_nlse

  subroutine split_equations_nlse_w_pot(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (abs(term))
    case (1)
       call evolve_gradient_trap_real(this,dt)
    case (2)
       call evolve_self_scattering(this,dt)
    case (3)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_nlse_w_pot

  subroutine split_equations_gpe_2fld(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (abs(term))
    case (1)
       call evolve_gradient_trap_real(this,dt)
    case (2)
       call evolve_self_scattering(this,dt)
    case(3)
       call evolve_interspecies_conversion_single(this,dt,1)
    case(4)
       call evolve_interspecies_conversion_single(this,dt,2)
    case(5)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_gpe_2fld

  subroutine split_equations_gpe_nfld(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1,-1)
       call evolve_gradient_trap_real(this,dt)
    case (2,-2)
       call evolve_self_scattering(this,dt)
    case (3,-3)
       call evolve_interspecies_conversion(this,dt,term>0)
    case (4,-4)
       call evolve_gradient_trap_imag(this,dt)
    end select
  end subroutine split_equations_gpe_nfld
  
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

    integer :: l

    do l=1, this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,2,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       this%psi(XIND,1,l) = this%psi(XIND,1,l)   &
            - ( 0.5_dl*this%tPair%realSpace(XIND) - this%v_trap(XIND)*this%psi(XIND,2,l) )*dt
    enddo
  end subroutine evolve_gradient_trap_real

  subroutine evolve_gradient_trap_imag(this, dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: l

    do l=1, this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,1,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       this%psi(XIND,2,l) = this%psi(XIND,2,l) &
            + ( 0.5_dl*this%tPair%realSpace(XIND) - this%v_trap(XIND)*this%psi(XIND,1,l) )*dt
    enddo
  end subroutine evolve_gradient_trap_imag

  ! Performance note : Precomputing the sine and cosine seems to run about twice as fast
  subroutine evolve_self_scattering(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: l
    real(dl) :: g_loc, mu_loc
    real(dl), dimension(XIND) :: phase_shift
    real(dl), dimension(XIND) :: co, sn
    
    mu_loc = mu
    
    do l = 1,this%nfld
       g_loc = g_self(l)
       phase_shift = this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2
       phase_shift = (g_loc*phase_shift - mu)*dt
       co = cos(phase_shift); sn = sin(phase_shift)
       
       this%tPair%realSpace = this%psi(XIND,2,l)
       this%psi(XIND,2,l) = co*this%psi(XIND,2,l) - sn*this%psi(XIND,1,l)
       this%psi(XIND,1,l) = co*this%psi(XIND,1,l) + sn*this%tPair%realSpace(XIND)
    enddo
  end subroutine evolve_self_scattering

  !>@brief
  !> Evolves the self-scattering coupling of the GPE, but first splits into
  !> amplitude and phase, then rotates phase, then projects back rather
  !> than directly applying a rotation matrix
  subroutine evolve_self_scattering_amp_phase(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    real(dl), dimension(XIND) :: amp, phase
    real(dl) :: g_loc, mu_loc
    integer :: l

    mu_loc = mu

    do l = 1,this%nfld
       g_loc = g_self(l)
       amp = this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2
       phase = atan2( this%psi(XIND,2,l), this%psi(XIND,1,l) )

       phase = phase + (g_loc*amp-mu)*dt
       amp = sqrt(amp)
       this%psi(XIND,1,l) = amp*cos(phase)
       this%psi(XIND,2,l) = amp*sin(phase)
    enddo
  end subroutine evolve_self_scattering_amp_phase
    
  subroutine evolve_interspecies_conversion_single(this,dt,fld_ind)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: fld_ind

    real(dl), dimension(1:this%nx,1:this%ny,1:2) :: dpsi
    real(dl), dimension(1:this%nfld) :: nu_cur
    integer :: l

    nu_cur = nu(:,fld_ind)
    nu_cur(fld_ind) = 0._dl
    
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

    ! The accumulation order needs to be checked here for performance.
    ! To reduce memory footprint, try storing only x-axis, and adding an inner loop over y-indices here
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

  subroutine evolve_interspecies_conversion_w_osc_single(this,dt,fld_ind)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: fld_ind

    real(dl), dimension(XIND,1:2) :: dpsi
    real(dl) :: dnu
    integer :: l

    nu_cur = nu(:,fld_ind)
    nu_cur(fld_ind) = 0._dl

    dpsi = 0._dl
    do l=1,this%nfld
       dnu = nu_cur(l)*dt + delta * ( sin(om*(tcur+dt)) - sim(om*tcur) )
       dpsi(XIND,1) = dpsi(XIND,1) - dnu * this%psi(XIND,2,l)
       dpsi(XIND,2) = dpsi(XIND,2) + dnu * this%psi(XIND,1,l)
    enddo
    this%psi(XIND,1:2,fld_ind) = this%psi(XIND,1:2,fld_ind) + dpsi
  end subroutine evolve_interspecies_conversion_w_osc_single
  
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

  ! Currently doing a lot of potentially unnecessary matrix creation
  subroutine evolve_local_evolution(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt


    integer, parameter :: nit = 8, order = 5
!    real(dl), parameter :: a(order,order) =  ! copy this
!    real(dl), parameter :: b(order) =        ! copy this
    
    real(dl), dimension(1:2*this%nfld,order) :: g
    real(dl), dimension(1:2*this%nfld) :: y, dy
    integer :: i,o
    ! Write Gauss-Legendre integrator here to do local nonlinear evolution
    ! Need to define the appropriate matrices

    g = 0._dl
    do i=1,nit
       !g = matmul(g,a)
       do o=1,order
          !call derivs_local(y+g(:,o)*dt , g(:,i), this%nfld)
       enddo
    enddo
    
  end subroutine evolve_local_evolution

  subroutine derivs_local(yc,yp,nf)
    real(dl), dimension(1:2*nf), intent(in) :: yc
    real(dl), dimension(1:2*nf), intent(out) :: yp
    integer, intent(in) :: nf
    
    integer :: i

    do i=1,nf
       yp(2*i-1) = 0._dl
       yp(2*i)   = 0._dl
    enddo
  end subroutine derivs_local
  
  
end module Equations
