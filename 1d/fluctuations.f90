module Fluctuations
  use constants, only : dl, twopi
  use Simulations
  use gaussianRandomField

  type SpectralParams
     integer :: n_cut
     real(dl) :: dk, rho
     integer :: type
  end type SpectralParams
  
contains

  !>@brief
  !> Given a pre-initialised set of fields, add vacuum fluctuations on top
  subroutine add_vacuum_fluctuations_homogeneous(this, rho)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: rho

    real(dl) :: df(1:this%nlat), spec(1:this%nlat/2+1)
    integer :: i,j
    integer :: ncut

    ncut = this%nlat/2
    spec = 0._dl

    do i=2,ncut
       spec(i) = 1._dl/sqrt(this%lSize*rho)
    enddo
    do i=1,this%nfld; do j=1,2
       call generate_1dGRF(df,spec,.false.)
       this%psi(1:this%nlat,i,j) = this%psi(1:this%nlat,i,j) + df
    enddo; enddo
  end subroutine add_vacuum_fluctuations_homogeneous

  !>@brief
  !> Given an uninitialized psi, generate vacuum fluctuations
  subroutine initialize_vacuum_fluctuations_homogeneous(this,rho)
    type(Lattice), intent(inout) :: this
    real(dl), intent(in) :: rho

    real(dl) :: df(1:this%nlat), spec(1:this%nlat/2+1)
    integer :: i,j

    spec = 0._dl
    do i=2,this%nlat/2
       spec(i) = 1._dl/sqrt(this%lSize*rho)
    enddo
    do i=1,this%nfld; do j=1,2
       call generate_1dGRF(df,spec,.false.)
       this%psi(1:this%nlat,i,j) = df
    enddo; enddo
  end subroutine initialize_vacuum_fluctuations_homogeneous
  
end module Fluctuations
