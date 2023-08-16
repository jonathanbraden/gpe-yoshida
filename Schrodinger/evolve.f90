#include "macros.h"

program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Yoshida
  use Equations 
  implicit none

  integer :: nf, nLat
  real(dl) :: lSize
  type(Lattice) :: mySim

  real(dl) :: dt
  integer :: i
  real(dl) :: trap_param, fld_norm, grad_ratio, m2_
  
  ! k_eff = 2/pi k_nyq = 2/dx => w2_eff = m2 + 4./dx**2
  
  nf = 1
  fld_norm = 1.
  nLat = 128

  !call evolve_vacuum_decay(nLat, 4, fld_norm)
  call evolve_preheating(256, 2, 2., 0.125*twopi)
  
contains

  subroutine test_sho(n, lSize, keff, diag)
    integer, intent(in) :: n
    real(dl), intent(in) :: lSize
    real(dl), intent(in) :: keff
    ! Not Yet Implemented
    logical, intent(in) :: diag  ! Whether or not to diagonalize the matrix

    type(Lattice) :: mySim
    integer :: nf
    real(dl) :: energy
    real(dl) :: grad_ratio ! convert to keff

    real(dl) :: dt
    integer :: i
    
    nf = 1
    call create_lattice_rectangle(mySim, (/n,n/), (/lSize,lSize/), nf)

    grad_ratio = 2.
    energy = 0.5_dl*(1._dl + sqrt(1.+grad_ratio))
    call initialize_trap(mySim,(/grad_ratio/),2,1._dl,grad_ratio)
    mySim%v_trap = mySim%v_trap - energy
    
    call imprint_gaussian_2d_diag( mySim, (/1., 1./sqrt(1.+grad_ratio)/) )

    call write_lattice_data(mySim, 50)
    
    dt = 1./512.
    do i=1,200
       call step_lattice(mySim,dt,64)
       call write_lattice_data(mySim,50)
       print*,"step ",i," time = ",mySim%time
    enddo
  end subroutine test_sho

  subroutine evolve_vacuum_decay(nLat, n_min, fld_norm)
    integer, intent(in) :: nLat, n_min
    real(dl), intent(in) :: fld_norm

    type(Lattice) :: mySim
    integer :: nf
    real(dl) :: trap_param, grad_ratio ! Move to inputs
    real(dl) :: m2_
    integer :: u
    real(dl) :: dt
    integer :: i
    
    ! Add n_min in here
    nf = 1
    call create_lattice_rectangle(mySim, (/nLat,nLat/),  &
         (/2.**0.5*fld_norm*4.*twopi,2.**0.5*fld_norm*4.*twopi/), nf)

    trap_param = 1.4
    grad_ratio = 0.2_dl !0.5 ! / 10. for mid
    m2_ = 1.-1._dl/trap_param**2
    !call initialize_trap(mySim,(/trap_param/),6, fld_norm, grad_ratio)
    !call imprint_gaussian_2d_diag( mySim, (/1./sqrt(m2_), 1./sqrt(m2_+0.25*grad_ratio)/) )

    call initialize_trap(mySim,(/trap_param/),7, fld_norm, grad_ratio)
    call imprint_gaussian_2d( mySim, (/1./sqrt(m2_), 1./sqrt(m2_+0.25*grad_ratio)/) )

    mySim%v_trap = mySim%v_trap - 0.5_dl*(sqrt(m2_) + sqrt(m2_+0.25*grad_ratio))

    u = 50
    call write_lattice_data(mySim, u)

    dt = 1./512.
    do i=1,400
       call step_lattice(mySim,dt,32)
       call write_lattice_data(mySim,u)
       print*,"step ",i," time = ",mySim%time
    enddo
 
  end subroutine evolve_vacuum_decay
  
  ! Evolve coherent state in sine-Gordon potential
  !
  ! Defaults until I work out better ones
  !  - nLat = 256
  !  - mean_fld = 0.125*twopi
  subroutine evolve_preheating(nLat, n_min, fld_norm, mean_fld)
    integer, intent(in) :: nLat, n_min
    real(dl), intent(in) :: fld_norm, mean_fld

    type(Lattice) :: mySim
    integer :: nf
    real(dl) :: grad_ratio  ! Make this an input parameter
    real(dl) :: m2_mean, m2_

    real(dl) :: dt
    integer :: i
    integer :: u
    
    nf = 1
    call create_lattice_rectangle(mySim, (/nLat,nLat/), (/2.**0.5*fld_norm*n_min*twopi,2.**0.5*fld_norm*n_min*twopi/), nf)

    grad_ratio = mean_fld**2 / 16.
    !grad_ratio = mean_fld**2 / 2.
    
    call initialize_trap(mySim,(/1._dl/), 3, fld_norm, grad_ratio)
    m2_mean = cos(mean_fld/fld_norm) 
    m2_ = m2_mean + grad_ratio
    
    call imprint_coherent_state( mySim, (/1./sqrt(m2_mean),1./sqrt(m2_)/), mean_fld*2.**0.5*fld_norm )
    !call imprint_gaussian_2d( mySim, (/1., 1./sqrt(m2_)/) )
    !call imprint_gaussian_2d_diag( mySim, (/1., 1./sqrt(m2_)/) )

    mySim%v_trap = mySim%v_trap - 0.5_dl*(1.+sqrt(m2_))

    u = 50  ! Fix this to be automated
    call write_lattice_data(mySim, u)

    dt = 1./512.
    do i=1,500
       call step_lattice(mySim,dt,128)
       call write_lattice_data(mySim,u)
       print*,"step ",i," time = ",mySim%time
    enddo
  end subroutine evolve_preheating

  
  ! Simulate tunneling in BEC potential
  subroutine evolve_bec_tunneling(lVal, keff, fld_norm, nLat )
    real(dl), intent(in) :: lVal, keff, fld_norm
    integer, intent(in) :: nLat

    type(Lattice) :: mySim
    real(dl) :: m2_
    
    call create_lattice_rectangle(mySim, (/nLat,nLat/), (/fld_norm*4.*twopi,fld_norm*4.*twopi/), nf)

    ! Fix this to use keff
    grad_ratio = 0.2_dl !0.5 ! / 10. for mid
    m2_ = 1.-1._dl/lVal**2
    call initialize_trap(mySim,(/lVal, grad_ratio/), 6, fld_norm, grad_ratio)
    mySim%v_trap = mySim%v_trap - 0.5_dl*(sqrt(m2_) + sqrt(m2_+grad_ratio))
    
    call imprint_gaussian_2d_diag( mySim, (/1./sqrt(m2_), 1./sqrt(m2_+grad_ratio)/) )
    
  end subroutine evolve_bec_tunneling
  
  ! Move these to another file
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

  
  ! Imprint Gaussian initial state in the diagonalized quadratic basis
  ! Need to check all coefficients, etc. in here
  subroutine imprint_gaussian_2d_diag(this,sig2)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:2), intent(in) :: sig2
    real(dl), dimension(1:this%nx,this%nfld,2) :: coord
    real(dl) :: norm
    integer :: j,l

    norm = 1._dl / (0.5_dl**2*twopi**2*sig2(1)*sig2(2))**0.25
    this%psi = 0._dl
    do j = 1,this%ny
       do l=1,this%nfld
          coord(:,l,1) = (this%xGrid(:) + this%yGrid(j))/sqrt(2.)
          coord(:,l,2) = (this%yGrid(j) - this%xGrid(:))/sqrt(2.)
          
          this%psi(1:this%nx,j,1,l) = norm*exp( -0.5_dl*(coord(:,l,1)**2/sig2(1)+coord(:,l,2)**2/sig2(2)) )
       enddo
    enddo
    
  end subroutine imprint_gaussian_2d_diag
  
  subroutine imprint_transverse_gaussian_2d(this,sig2)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: sig2
    integer :: j,l

    this%psi = 0._dl
    do l=1,this%nfld
       do j=1,this%ny
          this%psi(1:this%nx,j,1,l) = exp(-0.5_dl*this%yGrid(j)**2/sig2) / (0.5_dl*twopi*sig2)**0.25
       enddo
    enddo
  end subroutine imprint_transverse_gaussian_2d

  ! Add momentum here
  subroutine imprint_coherent_state(this,sig2,x)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:2), intent(in) :: sig2
    real(dl), intent(in) :: x
    integer :: j,l

    this%psi = 0._dl
    do j = 1,this%ny
       do l=1,this%nfld
          this%psi(1:this%nx,j,1,l) = exp( -0.5_dl*((this%xGrid-x)**2/sig2(1)+this%yGrid(j)**2/sig2(2)) ) / (0.5_dl**2*twopi**2*sig2(1)*sig2(2))**0.25
       enddo
    enddo

    ! What was this for?  Computing norm or something?
    this%tPair%realSpace = this%psi(1:this%nx,1:this%ny,1,1)**2 + this%psi(1:this%nx,1:this%ny,2,1)**2 
  end subroutine imprint_coherent_state
  
end program Evolve_GPE
