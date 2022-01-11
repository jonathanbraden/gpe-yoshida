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
  
  call create_lattice(mySim,128,4._dl,1)
  call set_model_parameters(1._dl,0.,0.,0.5_dl**2)
  call initialize_trap_potential(mySim, mySim%tPair%xGrid, 2._dl)
 
!  call imprint_sech_eigen(mySim,1)
!  call imprint_gaussian(mySim,1.)
  call imprint_gray_soliton(mySim,0.5_dl,0.125_dl*twopi)

  dx = mySim%dx
  alpha = 16.
  dt = dx**2/alpha
  call write_lattice_data(mySim,50)
  do i=1,500
     call step_lattice(mySim,dt,100) ! 100
     call write_lattice_data(mySim,50)
  enddo
  
contains
  
end program Evolve_GPE
