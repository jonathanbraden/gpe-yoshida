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
       0._dl, &  ! Fix this one
       0.74167036435061295344822780_dl, &
       -0.40910082580003159399730010_dl, &
       0.19075471029623837995387626_dl, &
       -0.57386247111608226665638733_dl, &
       0.29906418130365592384446354_dl, &
       0.33462491824529818378495798_dl, &
       0.31529309239676659663205666_dl  &
       /) 

  interface
     subroutine split_equations(this,dt,term,n_terms)
       type(Lattice), intent(inout) :: this
       real(dl), intent(in) :: dt
       integer, intent(in) :: term
       integer, intent(out), optional :: n_terms
     end subroutine split_equations
  end interface
  
contains
    
  !>@brief
  !> Take nstep steps of size dt on the lattice
  !
  !>@param[inout] this - The Lattice to advance
  !>@param[in]  dt - Step size dt
  !>@param[in]  nstep - Number of steps to take
  subroutine step_lattice_w_func(split,this,dt,nstep)
    procedure(split_equations) :: split
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nstep

    call symp8_w_func(this,dt,nstep,split)
    this%time = this%time + dt*nstep
  end subroutine step_lattice_w_func

  subroutine yoshida_scheme(this,dt,nsteps,w_coeffs)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nsteps
    real(dl), dimension(:), intent(in) :: w_coeffs

    integer :: n_stage, i, s

    n_stage = size(w_coeffs)

    do i=1,nsteps
       ! take 0.5*w(n_stage)*dt step on 1
       do s=n_stage,-n_stage,-1
       ! Fill in appropriate stepping here
       ! First, pick abs(i) as dt scaling and pick from w_coeffs
       ! will need to pick next w as well, this is a bit annoying with the way I set up the o2 step
          !
          ! w_max, w_max-1
          ! w_max-1, w_max-2
          ! ...
          ! w_1, w_0
          ! w_0, w_1
          ! ...
          ! w_max-1, w_max
          !
          ! If last step
          ! w_max, 0
          ! else
          ! w_max,w_max
       enddo
    enddo
  end subroutine yoshida_scheme

  !>@brief
  !> Second order accurate symplectic step with fusion
  !> w1 is time step of current step
  !> w2 is time step of next time step
  subroutine symp_o2_step_w_func(split,this,dt,w1,w2)
    procedure(split_equations) :: split
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt, w1, w2

    integer :: i

    do i=2,n_terms-1
       call split(this, 0.5_dl*w1*dt, i)
    enddo
    call split(this, w1*dt, n_terms)
    do i=n_terms-1,2,-1
       call split(this, 0.5_dl*w1*dt,i)
    enddo
    call split(this,0.5_dl*(w1+w2)*dt,1)
  end subroutine symp_o2_step_w_func
    
  subroutine symp2_w_func(split,this,dt,nsteps)
    procedure(split_equations) :: split
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nsteps

    integer :: i

    call split_equations(this,0.5_dl*dt,1)
    do i=1,nsteps-1
       call symp_o2_step_w_func(this,dt,1._dl,1._dl)
    enddo
    call symp_o2_step_w_func(this,dt,1._dl,0._dl)
  end subroutine symp2_w_func

  subroutine symp4_w_func(split,this,dt,nsteps)
    procedure(split_equations) :: split
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt
    integer, intent(in) :: nsteps

    real(dl), parameter :: w1 = 1._dl/(2._dl-2._dl**(1._dl/3._dl))
    real(dl), parameter :: w0 = 1._dl - 2._dl*w1
    integer :: i
    
    call split_equations(this,0.5_dl*w1*dt,1)
    do i=1,nsteps
       call symp_o2_step_w_func(this,dt,w1,w0)
       call symp_o2_step_w_func(this,dt,w0,w1)
       if (i == nsteps) then
          call symp_o2_step_w_func(this,dt,w1,0._dl)
       else
          call symp_o2_step_w_func(this,dt,w1,w1)
       endif
    enddo
  end subroutine symp4_w_func

  subroutine symp6_w_func(split,this,dt,nsteps)
    procedure(split_equations) :: split
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
       call symp_o2_step_w_func(this,dt,w3,w2)
       call symp_o2_step_w_func(this,dt,w2,w1)
       call symp_o2_step_w_func(this,dt,w1,w0)
       call symp_o2_step_w_func(this,dt,w0,w1)
       call symp_o2_step_w_func(this,dt,w1,w2)
       call symp_o2_step_w_func(this,dt,w2,w3)
       if (i==nsteps) then
          call symp_o2_step_w_func(this,dt,w3,0._dl)
       else
          call symp_o2_step_w_func(this,dt,w3,w3)
       endif
    enddo
  end subroutine symp6_w_func

  subroutine symp8_w_func(split,this,dt,nsteps)
    procedure(split_equations) :: split
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
       call symp_o2_step_w_func(split,this,dt, w7, w6)
       call symp_o2_step_w_func(split,this,dt, w6, w5)
       call symp_o2_step_w_func(split,this,dt, w5, w4)
       call symp_o2_step_w_func(split,this,dt, w4, w3)
       call symp_o2_step_w_func(split,this,dt, w3, w2)
       call symp_o2_step_w_func(split,this,dt, w2, w1)
       call symp_o2_step_w_func(split,this,dt, w1, w0)
       call symp_o2_step_w_func(split,this,dt, w0, w1)
       call symp_o2_step_w_func(split,this,dt, w1, w2)
       call symp_o2_step_w_func(split,this,dt, w2, w3)
       call symp_o2_step_w_func(split,this,dt, w3, w4)
       call symp_o2_step_w_func(split,this,dt, w4, w5)
       call symp_o2_step_w_func(split,this,dt, w5, w6)
       call symp_o2_step_w_func(split,this,dt, w6, w7)
       if (i==nsteps) then
          call symp_o2_step_w_func(split,this,dt, w7, 0._dl)
       else
          call symp_o2_step_w_func(split,this,dt, w7, w7)
       endif
    enddo
  end subroutine symp8_w_func
  
end module yoshida
