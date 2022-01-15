program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Initial_Conditions
  use gaussianRandomField
  use Yoshida
  use Equations ! remove this in Yoshida module
  implicit none

  type(Lattice) :: mySim
  integer :: i, nf
  real(dl) :: dx, dt, alpha
  real(dl), dimension(:), allocatable :: psi_i

  nf = 2
  !  call create_lattice(mySim,128,16._dl,1)
  call create_lattice(mySim,128,8._dl,nf)
!  call set_model_parameters(-1.,0.,0.,0.)
  call set_model_parameters(1.,0.,0.01,1.,nf)  
  call initialize_trap_potential(mySim, 1._dl,1)

!  call imprint_gaussian(mySim,1._dl)
  !call imprint_sech_eigen(mySim,1)
!  call imprint_bright_soliton(mySim,2.,3.)
  !call imprint_sine(mySim,2,1._dl)
  call imprint_mean_relative_phase(mySim,0.01_dl*twopi)

  dx = mySim%dx
  alpha = 8.
  dt = dx**2/alpha
  call write_lattice_data(mySim,50)
  do i=1,100
     call step_lattice(mySim,dt,50)
     call write_lattice_data(mySim,50)
  enddo
  
contains
  
end program Evolve_GPE
