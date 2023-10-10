#define XIND 1:this%nx,1:this%ny

module Initial_Conditions_2d
  use constants, only : dl
  use Simulation, only : Lattice
  implicit none

contains

  subroutine imprint_gaussian_state(this,x0,y0,sig)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: x0, y0, sig

    
  end subroutine imprint_gaussian_state
  
  !>@brief
  !> Imprint approximate vortex with winding number 1
  !
  ! \psi = r/sqrt(2+r**2)e^{i\phi} for wind = 1
  !      = (x+iy)/sqrt(2+r**2)
  !
  ! To Do : Need to add factors for different condensate normalization
  !         and different g normalizations
  subroutine imprint_vortex(this,x0,y0)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: x0, y0

    real(dl) :: rho2, y, x
    integer :: i,j, l

    do l=1,this%nfld
       do j=1,this%ny
          y = this%yGrid(j) - y0
          do i=1,this%nx
             x = this%xGrid(i) - x0
             rho2 = x**2 + y**2
             this%psi(i,j,1,l) = x/sqrt(2._dl+rho2)
             this%psi(i,j,2,l) = y/sqrt(2._dl+rho2)
          enddo
       enddo
    enddo  
  end subroutine imprint_vortex

  subroutine imprint_mean_relative_phase(this,rho,phi)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: rho, phi
    
    integer :: i,j,l

    this%psi(XIND,1,1) = sqrt(rho)
    this%psi(XIND,2,1) = 0._dl
    this%psi(XIND,1,2) = this%psi(XIND,1,1)*cos(phi)
    this%psi(XIND,2,2) = this%psi(XIND,1,1)*sin(phi)
  end subroutine imprint_mean_relative_phase
    
  !>@brief
  !> Add white noise to the condensate fields
  !> The RMS fluctuation amplitude of each field is set by rms.  The RMS of each real part is thus rms/sqrt(2)
  subroutine add_white_noise(this,rms)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: rms

    integer :: i,l
    real(dl), dimension(1:this%nx) :: phase, amp

!    call initialize_rand(42,13)

    do l=1,this%nfld
       do i=1,this%ny
          call random_number(phase); call random_number(amp)
          this%psi(:,i,1,l) = this%psi(:,i,1,l) + rms*sqrt(-log(amp))*cos(twopi*phase)
          this%psi(:,i,2,l) = this%psi(:,i,2,l) + rms*sqrt(-log(amp))*sin(twopi*phase)
       enddo
    enddo
  end subroutine add_white_noise
    
end module Initial_Conditions_2d
