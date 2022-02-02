#include "macros.h"

! TO DO : Move the imprint Gaussian etc. into a separate file
! which I can share with the time evolution program
program Solve_Ground_State
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Equations_Imag
  use Equations, only : initialize_trap
  implicit none

  integer :: i, nf
  type(Lattice) :: mySim
  integer :: u, v
  real(dl) :: error
  integer :: ierror
  real(dl) :: g_step

  g_step = 5._dl
  nf = 1
  call create_lattice_rectangle(mySim,(/128,128/),(/32.,32./),nf)
  call initialize_model_parameters(nf)
  call set_model_parameters(0.,0.,0.,0.)
  call initialize_trap(mySim,(/16./),3)

  call imprint_gaussian_2d(mySim,(/1._dl,0.25_dl/))
  open(unit=newunit(u),file='ground-states.bin',access='stream')
  write(u) mySim%psi

  open(unit=newunit(v),file='chem_pot.dat')
  write(v,*) "# Chemical Potential for NLSE in asymmetric harmonic trap"
  write(v,*) "# V = 1/2*(x^2+16*y^2)"
  write(v,*) "# g/(hbar*omega)  mu   energy"
  write(v,*) "# Reference:  Gaussian initial state"
  write(v,*) "# mu = ", chemical_potential_full(mySim)," E = ", energy(mySim)
  do i=0,50
     call set_model_parameters(i*g_step,0.,0.,0.)
     call solve_background_w_grad_flow(mySim, 1.e-15, 1.e-15, error, ierror)
     print*,"chemical potential is ",chemical_potential_full(mySim)
     write(u) mySim%psi
     write(v,*) i*g_step, chemical_potential_full(mySim), energy(mySim)
  enddo
  
contains

    subroutine imprint_gaussian_2d(this,sig2)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:2), intent(in) :: sig2
    integer :: j,l

    this%psi = 0._dl
    do j = 1,this%ny
       do l=1,this%nfld
          this%psi(1:this%nx,j,1,l) = exp( -0.5_dl*(this%xGrid**2/sig2(1)+this%yGrid(j)**2/sig2(2)) ) / (0.5_dl**2*twopi**2*sig2(1)*sig2(2))**0.25
       enddo
    enddo

    this%tPair%realSpace = this%psi(1:this%nx,1:this%ny,1,1)**2 + this%psi(1:this%nx,1:this%ny,2,1)**2
  end subroutine imprint_gaussian_2d
  
end program Solve_Ground_State
