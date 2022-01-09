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
  integer :: i
  real(dl) :: dx, dt, alpha
  
  call create_lattice(mySim,256,32._dl,2)
  !  call set_model_parameters(-1.,0.,0.,0.)
  call set_model_parameters(0.,0.,0.,0.)  ! No nonlinearity
  call initialize_trap_potential(mySim, 1._dl)

  !call imprint_gaussian(mySim,1._dl)
  call imprint_sech_eigen(mySim)
  !call imprint_bright_soliton(mySim,2.,3.)
  !call imprint_sine(mySim,2,1._dl)
  !call imprint_mean_relative_phase(mySim,0.1_dl*twopi)

  dx = mySim%dx
  alpha = 16.
  dt = dx**2/alpha
  call write_lattice_data(mySim,50)
  do i=1,500
     call step_lattice(mySim,dt,100)
     call write_lattice_data(mySim,50)
  enddo
  
contains
  
end program Evolve_GPE
