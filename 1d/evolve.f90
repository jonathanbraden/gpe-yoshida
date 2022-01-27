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
  integer :: ord

  ord = 4
  nf = 2
  call create_lattice(mySim,128,5._dl,nf)
  call set_model_parameters(250.,0.,250.*0.1_dl,0._dl,nf)  
  call initialize_trap_potential(mySim, 1._dl,3)

  call imprint_gaussian(mySim,1._dl)

  call solve_background_w_grad_flow(mySim, 1.e-15, 1.e-15)
  print*,"chemical potential is ",chemical_potential(mySim)
  print*,"full one is ", chemical_potential_full(mySim)
  print*,"field norm is ",field_norm(mySim)
  
  ! With the normalization I'm using, I don't get the correct
  ! time-evolution.  But I think this is due to my field normalization
  !  Check what happens if I divide by the field normalization
  ! No, problem is I don't have nu in the def. of chemical potential
  call set_chemical_potential(chemical_potential_full(mySim)/2.)
  
  dx = mySim%dx
  alpha = 8.
  dt = dx**2/alpha
  call write_lattice_data(mySim,50)
  do i=1,100
     call step_lattice(mySim,dt,25,ord)
     call write_lattice_data(mySim,50)
  enddo
  
contains
  
end program Evolve_GPE
