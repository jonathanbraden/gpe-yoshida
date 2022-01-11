program Solve_Background
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Equations_imag
  use Initial_Conditions
  implicit none

  type(Lattice) :: mySim
  integer :: i
  real(dl) :: dt
  
  call create_lattice(mySim,128,16._dl,1)
  call set_model_parameters(50.,0.,0.,0.)
  call initialize_trap_potential(mySim, 2._dl)

  call imprint_gaussian(mySim, 1._dl)

  open(unit=99,file='grad_descent.bin',access='stream')
  write(99) mySim%psi

  dt = mySim%dx**2/8.
  do i=1,200
     call gradient_flow(mySim,dt,50)
     write(99) mySim%psi
  enddo
  close(99)
  
contains


end program Solve_Background
