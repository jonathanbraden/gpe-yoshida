program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Initial_Conditions
  use gaussianRandomField
  use Yoshida
  use Equations ! remove this for yoshida
  implicit none

  type(Lattice) :: mySim
  integer :: i
  
  !  call create_lattice(mySim, 16, 4._dl, 2)
  call create_lattice(mySim,512,200._dl,2)
  call set_model_parameters(-1.,0.,0.,0.)
  
  call imprint_bright_soliton(mySim,0.9,1.)
  !call imprint_sine(mySim,2,1._dl)
  !call imprint_mean_relative_phase(mySim,0.1_dl*twopi)
  
  call write_lattice_data(mySim,50)
  do i=1,100
     call step_lattice(mySim,0.02,10)
     call write_lattice_data(mySim,50)
  enddo
  
contains
  
end program Evolve_GPE
