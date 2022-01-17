program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Initial_Conditions
  use gaussianRandomField
  use Yoshida
  use Equations
  use Equations_imag
  implicit none

  type(Lattice) :: mySim
  integer :: i, nf
  real(dl) :: dx, dt, alpha
  real(dl), dimension(:), allocatable :: psi_i

  nf = 1
  call create_lattice(mySim,256,10._dl,nf)
  call set_model_parameters(250.,0.,0.,0._dl,nf)  
  call initialize_trap_potential(mySim, 1._dl,3)

  call imprint_gaussian(mySim,1._dl)
  !call imprint_sech_eigen(mySim,1)
!  call imprint_bright_soliton(mySim,2.,3.)
  !call imprint_sine(mySim,2,1._dl)
!  call imprint_mean_relative_phase(mySim,0.01_dl*twopi)

  call solve_background_w_grad_flow(mySim, 1.e-15, 1.e-15)
  call set_chemical_potential(chemical_potential(mySim))
  
  dx = mySim%dx
  alpha = 8.
  dt = dx**2/alpha
  call write_lattice_data(mySim,50)
  do i=1,100
     call step_lattice(mySim,dt,25)
     call write_lattice_data(mySim,50)
  enddo
  
contains
  
end program Evolve_GPE
