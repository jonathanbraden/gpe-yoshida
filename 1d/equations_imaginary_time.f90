#define XIND 1:this%nlat
#include "macros.h"

module Equations_imag
  use constants, only : dl, twopi
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby
#endif
  use Model_Params
  use Simulation

  implicit none
   
contains

  ! To do : implement input of tolerances
  ! Add checkpointing option to write errors (and possibly field) to file
  subroutine solve_background_w_grad_flow(this, tol_psi, tol_grad)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: tol_psi, tol_grad
    
    real(dl) :: dtau, err_psi, err_grad, err
    integer :: i
    integer, parameter :: maxit = 200
    real(dl), parameter :: tol = 1.e-15
    logical, dimension(1:2) :: check_pt

    check_pt = .false.
    
    dtau = this%dx**2 / 8._dl
    do i=1,maxit
       err = gradient_flow(this,dtau,50)
       
       if (check_pt(1)) print*,"Error logging not implemented"
       if (check_pt(2)) print*,"Field logging not implemented"
       
       if (err < 1.e-15) then
          print*,"converged in ",i," steps of 50"
          exit
       endif
    enddo
  end subroutine solve_background_w_grad_flow

  real(dl) function gradient_flow(this,dtau,nsteps) result(error)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dtau
    integer, intent(in) :: nsteps

    real(dl), dimension(1:this%nlat,1:2,1:this%nfld) :: psi_prev
    integer :: i
    
    psi_prev(:,:,:) = this%psi(XIND,:,:)
    do i=1,nsteps
       call gradient_step(this,dtau)
       call renorm_field(this)
    enddo
    error = maxval(abs(this%psi(XIND,:,:)-psi_prev))
  end function gradient_flow
  
  subroutine split_equations_bg(this,dt,term)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: term

    select case (term)
    case (1)
       call evolve_heat_equation(this,dt)
    case (2)
       call evolve_nonlinear_interaction(this,dt)
    case (3)
       call renorm_field(this)
    end select
  end subroutine split_equations_bg
    
  subroutine evolve_heat_equation(this, dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    integer :: i, j

    do i = 1,this%nFld
       this%tPair%realSpace = this%psi(XIND,1,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,1,i) = this%psi(XIND,1,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,1,i) )*dt

       this%tPair%realSpace = this%psi(XIND,2,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,2,i) = this%psi(XIND,2,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,2,i) )*dt
    enddo
  end subroutine evolve_heat_equation

  !>@brief
  !> Solve u_t = -\beta |u|^2u
  subroutine evolve_nonlinear_interaction(this, dt)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt

    real(dl), dimension(1:this%nlat) :: rho2
    integer :: i
    real(dl) :: g_loc

    do i=1,this%nfld
       g_loc = g_self(i)
       rho2 = this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2
       rho2 = rho2 / (1._dl + 2._dl*g_loc*rho2*dt)
       this%psi(XIND,1,i) = this%psi(XIND,1,i)*sqrt(rho2)
       this%psi(XIND,2,i) = this%psi(XIND,2,i)*sqrt(rho2)
    enddo
  end subroutine evolve_nonlinear_interaction

  ! Fix this for multiple fields
  subroutine renorm_field(this)
    type(Lattice), intent(inout) :: this

    real(dl) :: norm
    integer :: i
    
    do i=1,this%nfld
#if defined(PERIODIC)
       norm = sum(this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2) * this%dx
#elif defined(INFINITE)
       norm = sum( (this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2)*this%tPair%quad_weights )
#endif
       this%psi(XIND,:,i) = this%psi(XIND,:,i) / sqrt(norm)
    enddo
  end subroutine renorm_field

  ! Fix this for the interacting field case
  ! Currently just for nonlinear Schrodinger
  subroutine gradient_step(this,dtau)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dtau

    integer :: i
    real(dl), dimension(1:this%nlat) :: rho2
    real(dl) :: g_loc

    do i=1,this%nfld
       g_loc = g_self(i)
       rho2 = this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2
       
       this%tPair%realSpace = this%psi(XIND,1,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,1,i) = this%psi(XIND,1,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,1,i) - g_loc*rho2*this%psi(XIND,1,i) ) *dtau

       
       this%tPair%realSpace = this%psi(XIND,2,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,2,i) = this%psi(XIND,2,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,2,i) - g_loc*rho2*this%psi(XIND,2,i) ) * dtau
    enddo
       
  end subroutine gradient_step
  
end module Equations_Imag
