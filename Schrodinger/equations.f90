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
  subroutine initialize_trap(this, params, pot, fld_norm, k2, diag_, cut_)

    abstract interface
       function pot_func(x,par) result(pot)
         import :: dl
         real(dl), intent(in) :: x
         real(dl), dimension(:), intent(in) :: par
         real(dl) :: pot
       end function pot_func
    end interface
    
    type(Lattice), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: params
    procedure(pot_func) :: pot
    real(dl), intent(in) :: fld_norm, k2
    logical, intent(in), optional :: diag_ 
    real(dl), intent(in), optional :: cut_

    real(dl) :: cut
    logical :: diag
    integer :: i, j
    real(dl), dimension(1:this%nx) :: xp, xm
    real(dl), dimension(1) :: dummy

    diag=.true.; if (present(diag_)) diag=diag_
    cut = 32._dl; if (present(cut_)) cut=cut_

    if (allocated(this%v_trap)) deallocate(this%v_trap)
    allocate(this%v_trap(1:this%nx,1:this%ny))
    this%v_trap = 0._dl
    
    do j=1,this%ny
       if (diag) then
          xp = sqrt(0.5)*( this%yGrid(j) + this%xGrid(:) ) / fld_norm
          xm = sqrt(0.5)*( this%yGrid(j) - this%xGrid(:) ) / fld_norm
       else
          xp = this%xGrid(:) / fld_norm
          xm = this%yGrid(j) / fld_norm
       endif

       do i=1,this%nx
          this%v_trap(i,j) = pot(xm(i),dummy) + pot(xp(i),dummy)
       enddo
    enddo
    this%v_trap = fld_norm**2*this%v_trap

    call add_trap_gradient_energy(this, 2._dl/sqrt(k2), diag)
    call normalize_trap(this, cut)  
    
    call output_trap(this)
  end subroutine initialize_trap
  
  
  ! Improvements for this subroutine
  ! Write the gradient energy as one term
  ! Write the potential as a different term
  ! Combine them.  This makes it easier to write the code
  !
  ! Add option to 
  subroutine initialize_trap_old(this, params, type, fld_norm, k2)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(:), intent(in) :: params
    integer, intent(in) :: type
    real(dl), intent(in) :: fld_norm, k2

    integer :: i, j
    real(dl), dimension(1:this%nx) :: xp, xm
    logical :: diag ! Put this as an input flag
    
    diag = .true.  ! Turn this into a variable
    
    if (allocated(this%v_trap)) deallocate(this%v_trap)
    allocate(this%v_trap(1:this%nx,1:this%ny))
    this%v_trap = 0._dl

    ! Modularise this to take the potential as a function pointer
    
    select case (type)
    case (1)
       this%v_trap = 0._dl

    case(2) ! Harmonic trap with gradient
       do j=1,this%ny
          if (diag) then
             xp = sqrt(0.5)*( this%yGrid(j) + this%xGrid(:) ) / fld_norm
             xm = sqrt(0.5)*( this%yGrid(j) - this%xGrid(:) ) / fld_norm
          else
             xp = this%xGrid(:) / fld_norm
             xm = this%yGrid(j) / fld_norm
          endif
          
          this%v_trap(:,j) = 0.5_dl*( xp**2 + xm**2 )
       enddo
       this%v_trap = fld_norm**2 * this%v_trap
       
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
             
          this%v_trap(:,j) = 1. - cos(xp) + 1. - cos(xm) - 2. ! centering
       enddo
       this%v_trap = fld_norm**2 * this%v_trap

       call add_trap_gradient_energy(this, 2./sqrt(k2), diag)
       call normalize_trap(this, 128._dl)
       
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
       this%v_trap = fld_norm**2 * this%v_trap ! I think this is in the wrong place ... (and should go after I add the gradient energy, or maybe not since it's quadratic ...)
       
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
  end subroutine initialize_trap_old
  
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

  subroutine split_equations_schrodinger_pml(this, dt, term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (abs(term))
    case (1)
       call evolve_gradient_trap_real(this, dt)
    case (2)
       call evolve_gradient_trap_imag(this, dt)
    !case (3)
    !   call evolve_potential_damp(this, dt)
    !case (4)
    !   call evolve
       
    end select
  end subroutine split_equations_schrodinger_pml
  
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

! Debug this loop
#ifdef LOOP
!$OMP PARALLEL DO FIRSTPRIVATE(dt,l) PRIVATE(i,j)
       do j=1,this%ny
          do i=1,this%nx
             this%psi(i,j,1,l) = this%psi(i,j,1,l) &
                  + ( 0.5_dl*this%tPair%realSpace(i,j) - this%v_trap(i,j)*this%psi(i,j,2,l) ) *dt
!$OMP END PARALLEL DO       
#else
       this%psi(XIND,2,l) = this%psi(XIND,2,l) &
            + ( 0.5_dl*this%tPair%realSpace(XIND) - this%v_trap(XIND)*this%psi(XIND,1,l) )*dt
#endif
    enddo
  end subroutine evolve_gradient_trap_imag

  real(dl) function pot_quad(x,params) result(pot)
    real(dl), intent(in) :: x
    real(dl), dimension(:), intent(in) :: params

    real(dl) :: m2

    m2 = params(1)
    
    pot = 0.5_dl*m2*x**2
  end function pot_quad

  real(dl) function pot_bec(x,params) result(pot)
    real(dl), intent(in) :: x
    real(dl), dimension(:), intent(in) :: params

    real(dl) :: lVal

    pot = cos(x) + 1._dl + 0.5_dl*lVal**2*sin(x)**2
  end function pot_bec

  real(dl) function pot_bec_norm(x,params) result(pot)
    real(dl), intent(in) :: x
    real(dl), dimension(:), intent(in) :: params

    real(dl) :: lVal

    pot = 0.5_dl*sin(x)**2 + (cos(x)-1._dl)/lVal**2
  end function pot_bec_norm
  
  real(dl) function pot_hertzberg(x,params) result(pot)
    real(dl), intent(in) :: x
    real(dl), dimension(:), intent(in) :: params

    pot = 0.5_dl*x**2*(1._dl-0.5_dl*x**2)/(1._dl+0.5*x**4)
  end function pot_hertzberg

  real(dl) function pot_sine_gordon(x,params) result(pot)
    real(dl), intent(in) :: x
    real(dl), dimension(:), intent(in) :: params

    pot = 1._dl - cos(x) 
  end function pot_sine_gordon

  real(dl) function pot_double_well(x,params) result(pot)
    real(dl), intent(in) :: x
    real(dl), dimension(:), intent(in) :: params

    pot = 0.25*(x**2-1._dl)**2
  end function pot_double_well
  
end module Equations
