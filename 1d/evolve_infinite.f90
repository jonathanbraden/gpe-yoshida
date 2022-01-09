program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Initial_Conditions
!  use gaussianRandomField
  use Yoshida
  use Equations ! remove this in Yoshida module
  implicit none

  type(Lattice) :: mySim
  integer :: i
  real(dl) :: dx, dt, alpha
  
  call create_lattice(mySim,64,4._dl,1)
  call set_model_parameters(0.,0.,0.,0.)
  call initialize_trap_potential(mySim, mySim%tPair%xGrid, 2._dl)
 
  call imprint_sech_eigen(mySim,1)

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
