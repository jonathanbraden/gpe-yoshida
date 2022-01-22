#include "macros.h"

module Equations_Imag
  use constants, only : dl, twopi
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby_2d
#endif
  use Model_Params
  use Simulation

  implicit none

contains

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
       print*,"err is ", err

       if (check_pt(1)) print*,"Error logging not implemented"
       if (check_pt(2)) print*,"Field logging not implemented"

       if (err < tol) then
          print*,"converged in ",i," steps of 50"
          exit
       endif
    enddo
  end subroutine solve_background_w_grad_flow

  real(dl) function gradient_flow(this,dtau,nsteps) result(error)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dtau
    integer, intent(in) :: nsteps

    real(dl), dimension(XIND,1:2,1:this%nfld) :: psi_prev
    integer :: i

    psi_prev(:,:,:,:) = this%psi(XIND,:,:)
    do i=1,nsteps
       call gradient_step(this,dtau)
       call renorm_field(this)
    enddo
    error = maxval(abs(this%psi(XIND,:,:)-psi_prev))
  end function gradient_flow

  ! Currently just for nonlinear schrodinger
  subroutine gradient_step(this,dtau)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dtau

    integer :: i,l
    real(dl), dimension(1:this%nx) :: rho2
    real(dl) :: g_loc

    do l=1,this%nfld
       g_loc = g_self(l)
       
       this%tPair%realSpace = this%psi(XIND,1,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_mapped(this%tPair)
#endif
       do i=1,this%ny
          rho2 = this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2

          this%psi(1:this%nx,i,1,l) = this%psi(1:this%nx,i,1,l) + &
               ( &
               0.5_dl*this%tPair%realSpace(1:this%nx,i) &
               - v_trap(1:this%nx,i)*this%psi(1:this%nx,i,1,l) &
               - g_loc*rho2(1:this%nx)*this%psi(1:this%nx,1,i) &
               ) * dtau 
       enddo

       this%tPair%realSpace = this%psi(XIND,2,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair, this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_mapped(this%tPair)
#endif
       do i=1,this%ny
          rho2 = this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2

          this%psi(1:this%nx,i,2,l) = this%psi(1:this%nx,i,2,l) + &
               ( &
               0.5_dl*this%tPair%realSpace(1:this%nx,i) &
               - v_trap(1:this%nx,i)*this%psi(1:this%nx,i,2,l) &
               - g_loc*rho2(1:this%nx)*this%psi(1:this%nx,2,i) &
               ) * dtau 
       enddo
    enddo
  end subroutine gradient_step

  subroutine renorm_field(this)
    type(Lattice), intent(inout) :: this

    real(dl) :: norm, norm_loc
    integer :: l

    norm = 0._dl
    do l=1,this%nfld
#if defined(PERIODIC)
       norm_loc = sum( this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2 ) * this%dx**2  ! fix to dx*dy
#elif defined(INFINITE)
       norm_loc = 0._dl  ! fix this one
#endif
       norm = norm + norm_loc
    enddo
    this%psi(XIND,1:2,1:this%nfld) = this%psi(XIND,1:2,1:this%nfld) / sqrt(norm)
  end subroutine renorm_field
    
end module Equations_Imag
