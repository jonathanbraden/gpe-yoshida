module yoshida
  use Equations
  implicit none

  ! Coefficients to be used to simplify the loops below
  real(dl), dimension(0:1), parameter :: w_o4 = (/ &
       -1.7024143839193115468333417084068059921_dl, &
       1.35120719195965577341667085420340299606_dl &
       /)
  
  real(dl), dimension(0:3), parameter :: w_o6 =  (/ &
       1.31518632068391121888424972823886251_dl, &
       -1.17767998417887100694641568096431573_dl, &
       0.235573213359358133684793182978534602_dl, &
       0.784513610477557263819497633866349876_dl &
       /)

  real(dl), dimension(0:7), parameter :: w_o8 = (/ &
       0._dl, &  ! Fix this one                            ! w0
       0.74167036435061295344822780_dl, &                  ! w1
       -0.40910082580003159399730010_dl, &                 ! w2
       0.19075471029623837995387626_dl, &
       -0.57386247111608226665638733_dl, &
       0.29906418130365592384446354_dl, &
       0.33462491824529818378495798_dl, &
       0.31529309239676659663205666_dl  &
       /) 

  type Yoshida_Integrator
     integer :: order
     integer :: n_terms
!   contains
!     procedure :: split => null()
  end type Yoshida_Integrator
  
contains
    
  !>@brief
  !> Take nstep steps of size dt on the lattice
  !
  !>@param[inout] this - The Lattice to advance
  !>@param[in]  dt - Step size dt
  !>@param[in]  nstep - Number of steps to take
  subroutine step_lattice(this,dt,nstep)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nstep

    call symp8(this,dt,nstep)
    this%time = this%time + dt*nstep
  end subroutine step_lattice

  subroutine yoshida_scheme(this,dt,nstep,w_coeffs)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nstep
    real(dl), dimension(:), intent(in) :: w_coeffs

    integer :: n_stage, i

    n_stage = size(w_coeffs)

    do i=-n_stage,n_stage
       ! Fill in appropriate stepping here
    enddo    
  end subroutine yoshida_scheme

  !>@brief
  !> Second order accurate symplectic step with fusion
  subroutine symp_o2_step(this,dt,w1,w2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt, w1, w2
    
    integer :: i

    do i=2,n_terms-1
       call split_equations(this, 0.5_dl*w1*dt, i)
    enddo
    call split_equations(this, w1*dt, n_terms)
    do i=n_terms-1,2,-1
       call split_equations(this, 0.5_dl*w1*dt,-i) ! sign flip in case we want to merge operators elsewhere
    enddo
    call split_equations(this,0.5_dl*(w1+w2)*dt,1) 
  end subroutine symp_o2_step

  !>@brief
  !> Second order accurate symplectic step with fusion
  !
  !> This version doesn't distinguish if we are using the first
  !> or second occurence of an operator
  subroutine symp_o2_step_(this,dt,w1,w2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt, w1, w2

    integer :: i

    do i=2,n_terms-1
       call split_equations(this, 0.5_dl*w1*dt, i)
    enddo
    call split_equations(this, w1*dt, n_terms)
    do i=n_terms-1,2,-1
       call split_equations(this, 0.5_dl*w1*dt,i)
    enddo
    call split_equations(this,0.5_dl*(w1+w2)*dt,1)
  end subroutine symp_o2_step_
  
  subroutine yoshida2(this,dt,w1,w2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    real(dl), intent(in) :: w1, w2

    integer :: i

    do i=2,n_terms-1
       call split_equations(this, 0.5_dl*w1*dt, i)
    enddo
    call split_equations(this, w1*dt, n_terms)
    do i=n_terms-1,2,-1
       call split_equations(this, 0.5_dl*w1*dt, i)
    enddo
    call split_equations(this, 0.5_dl*(w1+w2)*dt, 1)
  end subroutine yoshida2
  
  subroutine symp2(this,dt,nsteps)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nsteps

    integer :: i

    call split_equations(this,0.5_dl*dt,1)
    do i=1,nsteps-1
       call symp_o2_step(this,dt,1._dl,1._dl)
    enddo
    call symp_o2_step(this,dt,1._dl,0._dl)
  end subroutine symp2

  subroutine symp4(this,dt,nsteps)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nsteps

    real(dl), parameter :: w1 = 1._dl/(2._dl-2._dl**(1._dl/3._dl))
    real(dl), parameter :: w0 = 1._dl - 2._dl*w1
    integer :: i
    
    call split_equations(this,0.5_dl*w1*dt,1)
    do i=1,nsteps
       call symp_o2_step(this,dt,w1,w0)
       call symp_o2_step(this,dt,w0,w1)
       if (i == nsteps) then
          call symp_o2_step(this,dt,w1,0._dl)
       else
          call symp_o2_step(this,dt,w1,w1)
       endif
    enddo
  end subroutine symp4

  subroutine symp6(this,dt,nsteps)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nsteps

    !real(dl), parameter :: w1 = -1.17767998417887_dl
    !real(dl), parameter :: w2 = 0.235573213359357_dl
    !real(dl), parameter :: w3 = 0.784513610477560_dl
    !real(dl), parameter :: w0 = 1._dl-2._dl*(w1+w2+w3)

    real(dl), parameter :: w0 = 1.31518632068391121888424972823886251_dl
    real(dl), parameter :: w1 = -1.17767998417887100694641568096431573_dl
    real(dl), parameter :: w2 =  0.235573213359358133684793182978534602_dl
    real(dl), parameter :: w3 = 0.784513610477557263819497633866349876_dl

    integer :: i
    
    call split_equations(this,0.5_dl*w3*dt,1)
    do i=1,nsteps
       call symp_o2_step(this,dt,w3,w2)
       call symp_o2_step(this,dt,w2,w1)
       call symp_o2_step(this,dt,w1,w0)
       call symp_o2_step(this,dt,w0,w1)
       call symp_o2_step(this,dt,w1,w2)
       call symp_o2_step(this,dt,w2,w3)
       if (i==nsteps) then
          call symp_o2_step(this,dt,w3,0._dl)
       else
          call symp_o2_step(this,dt,w3,w3)
       endif
    enddo
  end subroutine symp6

  subroutine symp8(this,dt,nsteps)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nsteps

    real(dl), parameter :: w1 = 0.74167036435061295344822780_dl
    real(dl), parameter :: w2 = -0.40910082580003159399730010_dl
    real(dl), parameter :: w3 = 0.19075471029623837995387626_dl
    real(dl), parameter :: w4 = -0.57386247111608226665638733_dl
    real(dl), parameter :: w5 = 0.29906418130365592384446354_dl
    real(dl), parameter :: w6 = 0.33462491824529818378495798_dl
    real(dl), parameter :: w7 = 0.31529309239676659663205666_dl
    real(dl), parameter :: w0 = 1._dl - 2._dl*(w1+w2+w3+w4+w5+w6+w7)

    integer :: i
    
    call split_equations(this,0.5_dl*w7*dt,1)   
    do i=1,nsteps
       call symp_o2_step(this,dt, w7, w6)
       call symp_o2_step(this,dt, w6, w5)
       call symp_o2_step(this,dt, w5, w4)
       call symp_o2_step(this,dt, w4, w3)
       call symp_o2_step(this,dt, w3, w2)
       call symp_o2_step(this,dt, w2, w1)
       call symp_o2_step(this,dt, w1, w0)
       call symp_o2_step(this,dt, w0, w1)
       call symp_o2_step(this,dt, w1, w2)
       call symp_o2_step(this,dt, w2, w3)
       call symp_o2_step(this,dt, w3, w4)
       call symp_o2_step(this,dt, w4, w5)
       call symp_o2_step(this,dt, w5, w6)
       call symp_o2_step(this,dt, w6, w7)
       if (i==nsteps) then
          call symp_o2_step(this,dt, w7, 0._dl)
       else
          call symp_o2_step(this,dt, w7, w7)
       endif
    enddo
  end subroutine symp8
  
end module yoshida
