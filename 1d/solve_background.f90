program Solve_Background
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Equations_imag
  use Initial_Conditions
  implicit none

  type(Lattice) :: mySim
  integer :: i, j
  real(dl) :: dt
  real(dl) :: g_cur
  real(dl) :: err
  
  call create_lattice(mySim,128,5._dl,1)
!  call create_lattice(mySim,1024,128._dl,1)
  call set_model_parameters(250.,0.,0.,0.)
!  call initialize_trap_potential(mySim, 20._dl)
  call initialize_trap_potential(mySim, 1._dl)
  
  call imprint_gaussian(mySim, 1._dl)

  open(unit=99,file='grad_descent.bin',access='stream')
  open(unit=98,file='chemical_potential.dat')

  dt = mySim%dx**2/8.

  do j=0,250
     g_cur = j*1._dl
     call set_model_parameters(g_cur,0.,0.,0.)
     
     do i=1,200
        err =  gradient_flow(mySim,dt,50)
        if (err < 1.e-15) then
           print*,"Converged in ",i," steps of 50"
           exit
        endif
     enddo
     write(99) mySim%psi
     write(98,*) g_cur, chemical_potential(mySim)
  enddo
  close(99)
  close(98)
  
contains


end program Solve_Background
