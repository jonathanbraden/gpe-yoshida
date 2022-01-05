module Simulation
  use constants, only : dl, twopi
  use fftw3
  implicit none
  
  type Lattice
     real(dl), dimension(:,:,:), allocatable :: psi
     real(dl) :: time
     integer :: nlat, nFld
     real(dl) :: dx, lSize, dk
!#ifdef SPECTRAL
!#ifdef PERIODIC     
     type(transformPair1D) :: tPair
!#endif
#ifdef INFINITE
     type(chebyshevPair1D) :: tPair
#endif
!#endif
  end type Lattice
  
contains

  subroutine create_lattice(this,n,len,nf)
    type(Lattice), intent(out) :: this
    integer, intent(in) :: n, nf
    real(dl), intent(in) :: len

    this%time = 0._dl
    this%nlat = n; this%lSize = len
    this%nFld = nf
    this%dx = len/dble(n); this%dk = twopi/len
!#ifdef SPECTRAL
    call initialize_transform_1d(this%tPair, n)
!#else
#ifdef INFINITE
    call initialize_transform_cheby_1d(this%tPair,n,1)
#endif
!#ifdef SPECTRAL
    allocate( this%psi(1:n,1:2,1:nf) )
!#else
!    allocate( this%psi(0:n+1,1:2,1:nf) )
!#endif
!#endif
  end subroutine create_lattice

  subroutine destroy_lattice(this)
    type(Lattice), intent(inout) :: this

    this%time = -1._dl; this%nlat = -1; this%lSize = -1.
    this%dx = -1.; this%dk = -1.; this%nFld = -1
    if (allocated(this%psi)) deallocate(this%psi)
!#ifdef SPECTRAL
!#ifdef PERIODIC    
    call destroy_transform_1d(this%tPair)
!#endif
#ifdef INFINITE
    call destroy_transform_cheby_1d(this%tPair)
#endif
!#endif
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

    write(fNum) this%psi(1:this%nlat,:,:)    
  end subroutine write_lattice_data
  
  real(dl) function gradient_energy_spectral(this,direct) result(ge)
    type(Lattice), intent(inout) :: this
    logical, intent(in) :: direct
    integer :: n

    n = this%nlat
    ! Write this
  end function gradient_energy_spectral
  
end module Simulation