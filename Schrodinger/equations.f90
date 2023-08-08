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

  ! To Do : Pass in potential energy function
  ! To Do : Fix keff normalization
  ! To Do : Implement flag for doing coordinate vs diagonalized axes
  subroutine intialize_trap_(this, params, fld_norm, m2_norm, type, keff, cut, diag)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: params
    real(dl), intent(in) :: fld_norm, m2_norm
    integer, intent(in) :: type
    real(dl), intent(in) :: keff
    real(dl), intent(in) :: cut
    logical, intent(in) :: diag

    real(dl), dimension(1:this%nx) :: xc, yc
    real(dl) :: norm
    integer :: j

    if (allocated(this%v_trap)) deallocate(this%v_trap)
    allocate(this%v_trap(1:this%nx,1:this%ny))
    this%v_trap = 0._dl

    ! Can probably combine the final normalization into a single call
    select case (type)

    case(1)
       do j=1,this%ny
       
       enddo

    case(2)
       norm = 0.5_dl*m2_norm  ! Special case where fld_norm cancels out
       do j=1,this%ny
          this%v_trap(:,j) = norm*(this%xGrid(:)**2 + this%yGrid(j)**2)
       enddo

    case(3)
       norm = m2_norm*fld_norm**2 
       do j=1,this%ny
          this%v_trap(:,j) = 0.5_dl*sin(this%yGrid(j)/fld_norm)**2 &
               + (cos(this%yGrid(j)/fld_norm)-1._dl)/params(1)**2 &
               + 0.5_dl*sin(this%xGrid(:)/fld_norm)**2 &
               + (cos(this%xGrid(:)/fld_norm)-1._dl)/params(1)**2
       enddo
       this%v_trap = norm * this%v_trap

    case(4)
       norm = m2_norm*fld_norm**2
       do j=1,this%ny
          this%v_trap(:,j) = 1._dl - cos(this%yGrid(j)/fld_norm) &
               + 1._dl - cos(this%xGrid(:)/fld_norm)
       enddo
       this%v_trap = norm * this%v_trap

    case(5)
       norm = m2_norm*fld_norm**2
       do j=1,this%ny
          this%v_trap(:,j) = 0.25_dl*(this%yGrid(j)**2-fld_norm**2)**2 + 0.25_dl*(this%xGrid(:)**2-fld_norm**2)**2
       enddo
       this%v_trap = norm*this%v_trap
       
    case default
       this%v_trap = 0._dl
    end select
       
    call add_trap_gradient_energy(this, 1./keff**2, diag) ! fix keff norm
    call normalize_trap(this, cut)
  end subroutine intialize_trap_
  
  ! Improvements for this subroutine
  ! Write the gradient energy as one term
  ! Write the potential as a different term
  ! Combine them.  This makes it easier to write the code
  !
  ! Add option to 
  subroutine initialize_trap(this, params, type, fld_norm, k2)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: params
    integer, intent(in) :: type
    real(dl), intent(in) :: fld_norm, k2

    integer :: i, j
    real(dl), dimension(1:this%nx) :: xp, xm
    logical :: diag ! Put this as an input flag

    diag = .true.
    
    if (allocated(this%v_trap)) deallocate(this%v_trap)
    allocate(this%v_trap(1:this%nx,1:this%ny))
    this%v_trap = 0._dl
 
    select case (type)
    case (1)
       this%v_trap = 0._dl

    case(2) ! Harmonic trap with gradient
       do j=1,this%ny
          this%v_trap(:,j) = 0.5_dl*(this%xGrid(:)**2+this%yGrid(j)**2)  
       enddo
       call add_trap_gradient_energy(this, 2._dl/sqrt(k2), diag)
       call normalize_trap(this, 32._dl)

    case(3) ! Check normalization in here
       do j=1, this%ny
          if (diag) then
             xp = sqrt(0.5)*( this%yGrid(j) + this%xGrid(:) ) / fld_norm
             xm = sqrt(0.5)*( this%yGrid(j) - this%xGrid(:) ) / fld_norm
          else
             xp = this%xGrid(:) / fld_norm
             xm = this%yGrid(j) / fld_norm
          endif
             
          this%v_trap(:,j) = 1. - cos(xp) + 1. - cos(xm)
       enddo
       this%v_trap = fld_norm**2 * this%v_trap

       call add_trap_gradient_energy(this, 2./sqrt(k2), diag)
       call normalize_trap(this, 32._dl)
       
    case(5)  ! Drummond potential with gradient
       do j=1,this%ny
          this%v_trap(:,j) = cos(this%yGrid(j)/fld_norm) + 0.5_dl*params(1)**2*sin(this%yGrid(j)/fld_norm)**2 - 1._dl  &
               + cos(this%xGrid(:)/fld_norm) + 0.5_dl*params(1)**2*sin(this%xGrid(:)/fld_norm)**2 - 1._dl
       enddo
       this%v_trap = fld_norm**2 * this%v_trap
       call add_trap_gradient_energy( this, 2._dl/sqrt(k2), .false. ) 
       call normalize_trap(this, 32._dl)

    case(6)  ! Drummond with different norm convention (Further fix)
       do j=1, this%ny
          this%v_trap(:,j) = 0.5_dl*sin(this%yGrid(j)/fld_norm)**2 + (cos(this%yGrid(j)/fld_norm)-1._dl)/params(1)**2 &
                  + 0.5_dl*sin(this%xGrid(:)/fld_norm)**2 + (cos(this%xGrid(:)/fld_norm)-1._dl)/params(1)**2
       enddo
       this%v_trap(:,j) = fld_norm**2*this%v_trap(:,j)

       call add_trap_gradient_energy(this, 2./sqrt(k2), .false.)
       call normalize_trap(this, 32._dl)
       
    case(7)   ! Drummond in diagonal coordinates
       do j=1, this%ny
          xp = sqrt(0.5)*( this%yGrid(j) + this%xGrid(:) ) / fld_norm
          xm = sqrt(0.5)*( this%yGrid(j) - this%xGrid(:) ) / fld_norm
          
          this%v_trap(:,j) = 0.5_dl*sin(xp(:))**2 + (cos(xp(:))-1._dl)/params(1)**2  &
                  + 0.5_dl*sin(xm(:))**2 + (cos(xm(:))-1._dl)/params(1)**2
       enddo
       this%v_trap = fld_norm**2 * this%v_trap
       
       call add_trap_gradient_energy(this, 2./sqrt(k2), .true.)
       call normalize_trap(this, 32._dl)

    case(8) ! Double-well in diagonal coordinates (fix normalization)
       do j=1, this%ny
          xp = sqrt(0.5)*( this%yGrid(j) + this%xGrid(:) )
          xm = sqrt(0.5)*( this%yGrid(j) - this%xGrid(:) )

          this%v_trap(:,j) = 0.25*(xp**2/fld_norm**2-1._dl)**2 + 0.25*(xm**2/fld_norm**2-1._dl)**2
       enddo
       this%v_trap = fld_norm**2 * this%v_trap
       
       call add_trap_gradient_energy(this, 2./sqrt(k2), .true.)
       call normalize_trap(this, 32._dl)
       
    case default
       this%v_trap = 0._dl
    end select
    
    call output_trap(this)
  end subroutine initialize_trap
  
  ! Fix up normalisation here.
  !
  ! Add choice of convention for axis in phi_1, phi_2 or phi_+ and phi_-
  subroutine add_trap_gradient_energy(this, dx, diag)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dx
    logical, intent(in) :: diag

    real(dl) :: norm
    integer :: j

    norm = 1./dx**2
    
    if (.not.diag) then
       do j=1,this%ny
             this%v_trap(:,j) = this%v_trap(:,j) +  &
                  norm*(this%xGrid(:)-this%yGrid(j))**2
       enddo
    else ! Check this norm again
       do j=1,this%ny
          this%v_trap(:,j) = this%v_trap(:,j) + 2._dl*norm*this%yGrid(j)**2
       enddo
    endif
  end subroutine add_trap_gradient_energy

  ! Bug check element-wise max
  !
  ! Clip trap is values are above cut
  ! If smooth_ is true, then replace potential with cut*tanh(v_trap / cut)
  subroutine normalize_trap(this, cut, smooth_)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: cut
    logical, optional, intent(in) :: smooth_ 

    integer :: i,j
    logical :: smooth

    smooth = .false.; if (present(smooth_)) then; smooth=smooth_; endif 

    if (smooth) then
       this%v_trap(:,:) = cut*tanh(this%v_trap(:,:)/cut)
    else
        this%v_trap(:,:) = min(this%v_trap(:,:), cut)
    endif
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
