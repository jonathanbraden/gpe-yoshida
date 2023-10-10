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
    real(dl) :: mu_prev, en_prev, mu, en
    integer :: i
    integer, parameter :: maxit = 200
    real(dl), parameter :: tol = 5.e-15
    logical, dimension(1:2) :: check_pt
    integer, parameter :: it_size = 50
    
    check_pt = .false.
    
    dtau = this%dx**2 / 8._dl
    mu = chemical_potential(this); en = energy(this)
    
    do i=1,maxit
       mu_prev = mu; en_prev = en
       
       err = gradient_flow(this,dtau,it_size)
       
       mu = chemical_potential(this); en = energy(this)
       print*,"err is ",err,", mu is ",chemical_potential(this),", energy is ", energy(this)
       
       if (check_pt(1)) print*,"Error logging not implemented"
       if (check_pt(2)) print*,"Field logging not implemented"
       
       if (err < tol) then
          print*,"converged in ",i," steps of 50"
          exit
       endif
    enddo
    if (i == maxit) print*,"Failed to converge.  Error is ",err
  end subroutine solve_background_w_grad_flow

  real(dl) function gradient_flow(this,dtau,nsteps) result(error)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dtau
    integer, intent(in) :: nsteps

    real(dl), dimension(1:this%nlat,1:2,1:this%nfld) :: psi_prev
    integer :: i
    
    psi_prev(:,:,:) = this%psi(XIND,:,:)
    do i=1,nsteps
       call gradient_step_w_nu(this,dtau)
       call renorm_field(this)
    enddo
    error = maxval(abs(this%psi(XIND,:,:)-psi_prev))
  end function gradient_flow

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

  ! This one isn't converging to a solution yet
  ! This now converges, but I need to check the sign on the nu variable
  subroutine gradient_step_w_nu(this,dtau)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dtau

    integer :: i,l
    real(dl), dimension(XIND) :: rho2
    real(dl), dimension(XIND,1:2,1:this%nfld) :: psi_cur
    real(dl) :: g_loc, nu_loc(1:this%nfld)

    psi_cur(XIND,1:2,:) = this%psi(XIND,1:2,:)

    do i=1,this%nfld
       rho2 = this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2
       g_loc = g_self(i)

       
       this%tPair%realSpace = this%psi(XIND,1,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,1,i) = this%psi(XIND,1,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,1,i) - g_loc*rho2*this%psi(XIND,1,i) ) * dtau

       this%tPair%realSpace = this%psi(XIND,2,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,2,i) = this%psi(XIND,2,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,2,i) - g_loc*rho2*this%psi(XIND,2,i) ) * dtau
    enddo

    do i=1,this%nfld
       nu_loc = nu(:,i)
       ! Check the signs here
       do l=1,this%nfld
          this%psi(XIND,1,i) = this%psi(XIND,1,i) - nu_loc(l) * dtau * psi_cur(XIND,1,l)
          this%psi(XIND,2,i) = this%psi(XIND,2,i) - nu_loc(l) * dtau * psi_cur(XIND,2,l)
       enddo
    enddo
    
  end subroutine gradient_step_w_nu

  subroutine gradient_step_full(this, dtau)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dtau

    integer :: i,l
    real(dl), dimension(XIND) :: rho2, rho2_cross
    real(dl), dimension(XIND,1:2,1:this%nfld) :: psi_cur
    real(dl) :: g_loc, g_cross

    psi_cur(XIND,1:2,:) = this%psi(XIND,1:2,:)

    do i=1,this%nfld
       rho2 = this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2
       rho2_cross = this%psi(XIND,1,)**2 + this%psi(XIND,2,)**2
       g_loc = g_self(i)
       g_cross = 0._dl
       
       this%tPair%realSpace = this%psi(XIND,1,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,1,i) = this%psi(XIND,1,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,1,i)  &
            - g_loc*rho2*this%psi(XIND,1,i)                            &
            - 0.5_dl*g_cross*rho2_cross*this%psi(XIND,1,i) ) * dtau

       this%tPair%realSpace = this%psi(XIND,2,i)
#if defined(PERIODIC)
       call laplacian_1d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_1d_mapped(this%tPair)
#endif
       this%psi(XIND,2,i) = this%psi(XIND,2,i) + &
            ( 0.5_dl*this%tPair%realSpace - v_trap*this%psi(XIND,2,i)  &
            - g_loc*rho2*this%psi(XIND,2,i)                            &
            - 0.5_dl*g_cross*rho2_cross*this%psi(XIND,2,i) ) * dtau
    enddo
    
  end subroutine gradient_step_full
  
  ! Fix this for multiple fields
  ! Add required normalization constant as input
  ! Currently normalized to number of fields (for simplifying definitions of g, etc
  subroutine renorm_field(this)
    type(Lattice), intent(inout) :: this

    real(dl) :: norm, norm_loc, norm_tot
    integer :: i

    norm_tot = dble(this%nfld)
    
    norm = 0._dl
    do i=1,this%nfld
#if defined(PERIODIC)
       norm_loc = sum(this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2) * this%dx
#elif defined(INFINITE)
       norm_loc = sum( (this%psi(XIND,1,i)**2 + this%psi(XIND,2,i)**2)*this%tPair%quad_weights )
#endif
       norm = norm + norm_loc
    enddo
    this%psi(XIND,1:2,1:this%nfld) = sqrt(norm_tot)*this%psi(XIND,1:2,1:this%nfld) / sqrt(norm) 
  end subroutine renorm_field
  
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

    integer :: i

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
  
end module Equations_Imag
