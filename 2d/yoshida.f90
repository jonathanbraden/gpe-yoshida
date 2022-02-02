module yoshida
  use Equations
  implicit none

  real(dl), dimension(1:2) :: w_o4 = (/0._dl, 0._dl /)
!  real(dl), dimension(1:4), parameter :: w_o6 = (/ -1.17767998417887_dl, 0.235573213359357_dl, 0.784513610477560_dl, 1._dl-2._dl*(w1+w2+w3) /)
  real(dl), dimension(1:2) :: w_o8 = (/0._dl, 0._dl /)
 
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

    call symp4(this,dt,nstep)
  end subroutine step_lattice

  subroutine symp_o2_step(this,dt,w1,w2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: dt, w1, w2

    integer :: i

    do i=2,n_terms-1
       call split_equations(this, 0.5_dl*w1*dt, i)
    enddo
    this%time = this%time + 0.5_dl*w1*dt
    call split_equations(this, w1*dt, n_terms)
    !this%time_mid = this%time_mid + w1*dt
    do i=n_terms-1,2,-1
       call split_equations(this, 0.5_dl*w1*dt,-i)
    enddo
    this%time = this%time + 0.5_dl*w1*dt
    call split_equations(this,0.5_dl*(w1+w2)*dt,1)
    !this%time_first = this%time_first + 0.5_dl*(w1+w2)*dt
  end subroutine symp_o2_step
  
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

    real(dl), parameter :: w1 = -1.17767998417887_dl
    real(dl), parameter :: w2 = 0.235573213359357_dl
    real(dl), parameter :: w3 = 0.784513610477560_dl
    real(dl), parameter :: w0 = 1._dl-2._dl*(w1+w2+w3)

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
  end subroutine symp8
  
end module yoshida
