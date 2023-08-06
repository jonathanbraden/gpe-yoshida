#include "macros.h"

program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Yoshida
  use Equations 
  implicit none

  integer :: nf, nLat
  real(dl) :: lSize
  type(Lattice) :: mySim

  real(dl) :: dt
  integer :: i
  real(dl) :: trap_ratio, w_rat, grad_ratio, m2_
  
  nf = 1
  trap_ratio = 2.
  w_rat = 1. ! try various values to set hbar
  !call create_lattice_rectangle(mySim, (/256,256/), (/w_rat*4.*twopi,w_rat*4.*twopi/), nf)
  call create_lattice_rectangle(mySim, (/128,128/), (/32._dl,32._dl/), nf)

  ! For harmonic trap aligned on axes
  trap_ratio = 2._dl
  call initialize_trap(mySim,(/trap_ratio/),3)
  !call imprint_gaussian_2d( mySim, (/1.,1./sqrt(trap_ratio)/) )
  call imprint_coherent_state( mySim, (/1.,1./sqrt(trap_ratio)/), 3. )
  mySim%v_trap = mySim%v_trap - 0.5_dl*(1.+sqrt(trap_ratio))
  
  ! For harmonic trap not aligned on axes
  !grad_ratio = 2.
  !call initialize_trap(mySim,(/grad_ratio/),2)
  !call imprint_gaussian_2d_diag( mySim, (/1., 1./sqrt(1.+grad_ratio)/) )
  !mySim%v_trap = mySim%v_trap - 0.5_dl*(1.+sqrt(1.+grad_ratio))

  ! For BEC potential
  !trap_ratio = 1.4
  !grad_ratio = 0.2_dl !0.5 ! / 10. for mid
  !m2_ = 1.-1._dl/trap_ratio**2
  !call initialize_trap(mySim,(/trap_ratio,grad_ratio,w_rat/),6)
  !call imprint_gaussian_2d_diag( mySim, (/1./sqrt(m2_), 1./sqrt(m2_+grad_ratio)/) )
  !mySim%v_trap = mySim%v_trap - 0.5_dl*(sqrt(m2_) + sqrt(m2_+grad_ratio))

  !trap_ratio = 100.
  !grad_ratio = 0.
  !m2_ = trap_ratio**2-1._dl
  !call initialize_trap(mySim,(/trap_ratio,grad_ratio/),5)
  !call imprint_gaussian_2d_diag( mySim, (/1./sqrt(m2_), 1./sqrt(m2_+grad_ratio)/) )
  !mySim%v_trap = mySim%v_trap - 0.5_dl*(sqrt(m2_) + sqrt(m2_+grad_ratio))
  
  ! call solve_background_w_grad_flow()  ! If I want automation

  call write_lattice_data(mySim, 50)
  ! call set_chemical_potential() ! Replace this with E_gs and subtract off potential

  dt = 1./512.
  do i=1,200
     call step_lattice(mySim,dt,64)
     call write_lattice_data(mySim,50)
     print*,"step ",i," time = ",mySim%time
  enddo
  
contains

  subroutine imprint_gaussian_2d(this,sig2)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:2), intent(in) :: sig2
    integer :: j,l

    this%psi = 0._dl
    do j = 1,this%ny
       do l=1,this%nfld
          this%psi(1:this%nx,j,1,l) = exp( -0.5_dl*(this%xGrid**2/sig2(1)+this%yGrid(j)**2/sig2(2)) ) / (0.5_dl**2*twopi**2*sig2(1)*sig2(2))**0.25
       enddo
    enddo

    this%tPair%realSpace = this%psi(1:this%nx,1:this%ny,1,1)**2 + this%psi(1:this%nx,1:this%ny,2,1)**2
  end subroutine imprint_gaussian_2d

  subroutine imprint_transverse_gaussian_2d(this,sig2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: sig2
    integer :: j,l

    this%psi = 0._dl
    do l=1,this%nfld
       do j=1,this%ny
          this%psi(1:this%nx,j,1,l) = exp(-0.5_dl*this%yGrid(j)**2/sig2) / (0.5_dl*twopi*sig2)**0.25
       enddo
    enddo
  end subroutine imprint_transverse_gaussian_2d

  ! Add momentum here
  subroutine imprint_coherent_state(this,sig2,x)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:2), intent(in) :: sig2
    real(dl), intent(in) :: x
    integer :: j,l

    this%psi = 0._dl
    do j = 1,this%ny
       do l=1,this%nfld
          this%psi(1:this%nx,j,1,l) = exp( -0.5_dl*((this%xGrid-x)**2/sig2(1)+this%yGrid(j)**2/sig2(2)) ) / (0.5_dl**2*twopi**2*sig2(1)*sig2(2))**0.25
       enddo
    enddo

    ! What was this for?  Computing norm or something?
    this%tPair%realSpace = this%psi(1:this%nx,1:this%ny,1,1)**2 + this%psi(1:this%nx,1:this%ny,2,1)**2 
  end subroutine imprint_coherent_state
  
end program Evolve_GPE
