#define PERIODIC
!#define INFINITE

module Simulation
  use constants, only : dl, twopi
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby_2D
#endif
  implicit none
  
  type Lattice
     real(dl), dimension(:,:,:,:), allocatable :: psi
     real(dl), dimension(:), allocatable :: xGrid, yGrid
     real(dl) :: time
     integer :: nlat, nFld
     real(dl) :: dx, lSize, dk
!#ifdef SPECTRAL
     type(transformPair2D) :: tPair
!#endif
  end type Lattice
  
contains

  subroutine create_lattice(this,n,len,nf)
    type(Lattice), intent(out) :: this
    integer, intent(in) :: n, nf
    real(dl), intent(in) :: len

    integer :: i_
    
    this%time = 0._dl
    this%nlat = n; this%lSize = len
    this%nFld = nf
    allocate(this%xGrid(1:n), this%yGrid(1:n))

#if defined(PERIODIC)
    this%dx = len/dble(n); this%dk = twopi/len
    call initialize_transform_2d(this%tPair, (/n,n/) )
    this%xGrid = [ ( (i_-1)*this%dx, i_=1,n) ]
    this%yGrid = [ ( (i_-1)*this%dx, i_=1,n) ]
#elif defined(INFINITE)
    call initialize_transform_cheby_2d(this%tPair, (/n,n/), 1)
    call transform_rational_cheby(this%tPair, len)
    this%xGrid = this%tPair%xGrid
    this%yGrid = this%tPair%yGrid
#endif
!#ifdef SPECTRAL
    allocate( this%psi(1:n,1:n,1:2,1:nf) )
!#endif
  end subroutine create_lattice

  subroutine destroy_lattice(this)
    type(Lattice), intent(inout) :: this

    this%time = -1._dl; this%nlat = -1; this%lSize = -1.
    this%dx = -1.; this%dk = -1.; this%nFld = -1
    if (allocated(this%psi)) deallocate(this%psi)
#if defined(PERIODIC)
    call destroy_transform_2d(this%tPair)
#elif defined(INFINITE)
    call destroy_transform_cheby_2d(this%tPair)
#endif
  end subroutine destroy_lattice

  subroutine write_lattice_header(this,fNum)
    type(Lattice), intent(in) :: this
    integer, intent(in) :: fNum
    logical :: o

    inquire(opened=o,unit=fNum)
    if (o) then
       print*,"File header not head implemented"
    endif
  end subroutine write_lattice_header

  subroutine write_lattice_data(this,fNum)
    type(Lattice), intent(in) :: this
    integer, intent(in) :: fNum
    logical :: o

    inquire(opened=o,unit=fNum)
    if (.not.o) open(unit=fNum,file='fields.bin',access='stream')

    write(fNum) this%psi(1:this%nlat,1:this%nlat,:,:)
    
  end subroutine write_lattice_data
  
  real(dl) function gradient_energy_spectral(this,direct) result(ge)
    type(Lattice), intent(inout) :: this
    logical, intent(in) :: direct
    integer :: n

    n = this%nlat
  end function gradient_energy_spectral
  
end module Simulation
