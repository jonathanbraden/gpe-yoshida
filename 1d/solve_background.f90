program Solve_Background
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Equations_imag
  use Initial_Conditions
  implicit none

  type(Lattice) :: mySim
  integer :: i
  real(dl) :: err, dt
  
  call create_lattice(mySim,128,10._dl,1)
  call set_model_parameters(250.,0.,0.,0.,1)
  call initialize_trap_potential(mySim, 1._dl,3)
  
  call imprint_gaussian(mySim, 1._dl)
  call solve_background_w_grad_flow(mySim, 1.e-15, 1.e-15)

  open(unit=99,file='bg.dat')
  do i=1,mySim%nlat
     write(99,*) mySim%psi(i,:,:)
  enddo
  close(99)
  open(unit=99,file='mu.dat')
  write(99,*) chemical_potential(mySim)
  print*,"mu is ",chemical_potential(mySim)
  print*,"energy is ", energy(mySim)
  close(99)
  
!  call scan_g_values(mySim)
  
contains

  subroutine scan_g_values(sim)
    type(Lattice), intent(inout) :: sim

    integer :: i,j
    real(dl) :: g_cur, err
    real(dl) :: dt
    
    call imprint_gaussian(sim, 1._dl)
    
    open(unit=99,file='gpe_vacua.bin',access='stream')
    open(unit=98,file='chemical_potential.dat')

    dt = mySim%dx**2/8._dl
    
    do j=0,250
       g_cur = j*1._dl
       call set_model_parameters(g_cur,0.,0.,0.,1)
       
       do i=1,200
          err =  gradient_flow(sim,dt,50)
          if (err < 1.e-15) then
             print*,"Converged in ",i," steps of 50"
             exit
          endif
       enddo
       write(99) mySim%psi
       write(98,*) g_cur, chemical_potential(mySim), energy(mySim)
    enddo
    close(99)
    close(98)
  end subroutine scan_g_values

end program Solve_Background
