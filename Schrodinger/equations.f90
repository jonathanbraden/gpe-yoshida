#include "macros.h"

! Optimization choices
!#define LOOP T

module Equations
  use constants, only : dl, twopi
  use utils, only : newunit
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby_2D
#endif
  use Simulation
  
  implicit none

  integer, parameter :: n_terms = 2
 
contains

  ! Improvements for this subroutine
  ! Write the gradient energy as one term
  ! Write the potential as a different term
  ! Combine them.  This makes it easier to write the code
  !
  ! Add option to 
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
       
    case(2) ! Harmonic trap with gradient
       do j=1,this%ny
          this%v_trap(:,j) = min(0.5_dl*(this%xGrid(:)**2+this%yGrid(j)**2)  &
               + 0.25_dl*params(1)*(this%yGrid(j)-this%xGrid(:))**2, &
               32._dl)
       enddo
       
    case(3) ! Harmonic trap in both directions
       do j=1,this%ny
          this%v_trap(:,j) = min(0.5_dl*(this%xGrid(:)**2+params(1)*this%yGrid(j)**2),32.)
       enddo
       
    case(4)  ! Harmonic trap only the y-direction
       do j=1,this%ny
          this%v_trap(:,j) = min(0.5_dl*this%yGrid(j)**2,params(1))
       enddo
       
    case(5)  ! Drummond potential with gradient
       do j=1,this%ny
          do i=1,this%nx
             this%v_trap(i,j) = cos(this%yGrid(j)) + 0.5_dl*params(1)**2*sin(this%yGrid(j))**2 - 1._dl  &
                  + cos(this%xGrid(i)) + 0.5_dl*params(1)**2*sin(this%xGrid(i))**2 - 1._dl  &
                  + 0.25_dl*params(2)*(this%xGrid(i)-this%yGrid(j))**2
          enddo
       enddo

    case(6)  ! Drummond with different norm convention (Further fix)
       do j=1, this%ny
          do i=1,this%nx
             this%v_trap(i,j) = 0.5_dl*sin(this%yGrid(j)/params(3))**2 + (cos(this%yGrid(j)/params(3))-1._dl)/params(1)**2 &
                  + 0.5_dl*sin(this%xGrid(i)/params(3))**2 + (cos(this%xGrid(i)/params(3))-1._dl)/params(1)**2

             this%v_trap(i,j) = params(3)**2*this%v_trap(i,j)
             this%v_trap(i,j) = this%v_trap(i,j) + 0.25_dl*params(2)*(this%xGrid(i)-this%yGrid(j))**2
          enddo
       enddo

    case(7)  ! Add potential in mean and diff modes here
       this%v_trap = 0._dl
       
    case default
       this%v_trap = 0._dl
    end select

    call output_trap(this)
  end subroutine initialize_trap

  ! Fix up normalisation here.
  !
  ! Add choice of convention for axis in phi_1, phi_2 or phi_+ and phi_-
  subroutine add_trap_gradient_energy(this, dx)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dx

    integer :: i,j
    
    do j=1,this%ny
       do i=1,this%nx
          this%v_trap(i,j) = this%v_trap(i,j) +  &
               (1._dl/dx**2)*(this%xGrid(i)-this%yGrid(j))**2
       enddo
    enddo

    ! Add code here for phi_+, phi_- coordinate system
    ! the (this%xGrid(i) - this%yGrid(j)) -> this%yGrid(j) above
    
  end subroutine add_trap_gradient_energy

  ! Clip trap if the values get too big
  ! Change this to use where
  !
  ! Also add option to use cut*tanh(v_trap / cut)
  ! for a smooth potential
  subroutine normalize_trap(this, cut)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: cut

    integer :: i,j

    do j=1,this%ny; do i=1,this%nx
       this%v_trap(i,j) = max(this%v_trap(i,j), cut)
    enddo; enddo
    
  end subroutine normalize_trap
  
  subroutine output_trap(this)
    type(Lattice), intent(in) :: this
    integer :: i,u

    open(unit=newunit(u),file='trap.bin',access='stream')
    write(u) this%v_trap
    close(u)
    
    open(unit=newunit(u),file='xgrid.dat')
    do i=1,this%nx
       write(u,*) this%xGrid(i)
    enddo
    close(u)
    open(unit=newunit(u),file='ygrid.dat')
    do i=1,this%ny
       write(u,*) this%yGrid(i)
    enddo
    close(u)
  end subroutine output_trap
  
  subroutine split_equations(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (abs(term))
    case (1)
       call evolve_gradient_trap_real(this,dt)
    case (2)
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

  subroutine evolve_gradient_real(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: l
#ifdef LOOP
    integer :: i,j
#endif
    
    do l = 1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,2,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif

#ifdef LOOP
!$OMP PARALLEL DO FIRSTPRIVATE(dt,l) PRIVATE(i,j)
       do j=1,this%ny
          do i=1,this%nx
             this%psi(i,j,1,l) = this%psi(i,j,1,l) - 0.5_dl*this%tPair%realSpace(i,j)*dt
          enddo
       enddo
!$OMP END PARALLEL DO
#else
       this%psi(XIND,1,l) = this%psi(XIND,1,l) - 0.5_dl*this%tPair%realSpace(XIND)*dt
#endif
    enddo
  end subroutine evolve_gradient_real

  subroutine evolve_gradient_imag(this,dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: l
#ifdef LOOP
    integer :: i,j
#endif
    
    do l = 1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,1,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif

#ifdef LOOP
!$OMP PARALLEL DO FIRSTPRIVATE(dt,l) PRIVATE(i,j)
       do j=1,this%ny
          do i=1,this%nx
             this%psi(i,j,2,l) = this%psi(i,j,2,l) + 0.5_dl*this%tPair%realSpace(i,j)*dt
          enddo
       enddo
!$OMP END PARALLEL DO
#else
       this%psi(XIND,2,l) = this%psi(XIND,2,l) + 0.5_dl*this%tPair%realSpace(XIND)*dt
#endif
    enddo
  end subroutine evolve_gradient_imag

  subroutine evolve_gradient_trap_real(this, dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: l
#ifdef LOOP
    integer :: i,j
#endif
    
    do l = 1,this%nFld
       this%tPair%realSpace(XIND) = this%psi(XIND,2,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif

#ifdef LOOP
!$OMP PARALLEL DO FIRSTPRIVATE(dt,l) PRIVATE(i,j)
       do j=1,this%ny
          do i=1,this%nx
             this%psi(i,j,1,l) = this%psi(i,j,1,l) &
                  - ( 0.5_dl*this%tPair%realSpace(i,j) - this%v_trap(i,j)*this%psi(i,j,2,l) )*dt
          enddo
       enddo
!$OMP END PARALLEL DO
#else
       this%psi(XIND,1,l) = this%psi(XIND,1,l)   &
            - ( 0.5_dl*this%tPair%realSpace(XIND) - this%v_trap(XIND)*this%psi(XIND,2,l) )*dt
#endif
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
    
end module Equations
