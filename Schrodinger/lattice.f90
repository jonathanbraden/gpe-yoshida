#include "macros.h"

module Simulation
  use constants, only : dl, twopi
  use utils, only : newunit
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby_2D
#endif
  implicit none
  
  type Lattice
     real(dl), dimension(:,:,:,:), allocatable :: psi
     real(dl), dimension(:,:), allocatable :: v_trap
     real(dl), dimension(:), allocatable :: xGrid, yGrid
     real(dl), dimension(1:2) :: dx, lSize, dk
     real(dl) :: time, time_half
     integer :: nx, ny
     integer :: nFld
#if defined(PERIODIC)
     type(transformPair2D) :: tPair
#elif defined(INFINITE)
     type(chebyshevPair2D) :: tPair
#endif
  end type Lattice
  
contains

  subroutine create_lattice(this,n,len,nf)
    type(Lattice), intent(out) :: this
    integer, intent(in) :: n, nf
    real(dl), intent(in) :: len

    integer :: u
    integer :: i_
    
    this%time = 0._dl
    this%nx = n; this%ny = n
    this%lSize(1:2) = len
    this%nFld = nf

    ! To ensure vectorization, see what happens if I assign these
    ! using FFTW's initialization
    allocate( this%psi(1:n,1:n,1:2,1:nf) )
    allocate( this%v_trap(1:this%nx,1:this%ny) )    
    allocate( this%xGrid(1:this%nx), this%yGrid(1:this%ny) )

#if defined(PERIODIC)
    this%dx(1:2) = len/dble(n); this%dk(1:2) = twopi/len
    call initialize_transform_2d(this%tPair, (/n,n/) )
    this%xGrid = [ ( (i_-1)*this%dx(1), i_=1,n) ] - 0.5_dl*len
    this%yGrid = [ ( (i_-1)*this%dx(2), i_=1,n) ] - 0.5_dl*len
#elif defined(INFINITE)
    call initialize_transform_cheby_2d(this%tPair, (/n,n/), 1)
    call transform_rational_cheby(this%tPair, (/ len, len /))
    this%xGrid = this%tPair%xGrid
    this%yGrid = this%tPair%yGrid
    this%dx(1:2) = minval(abs(this%xGrid(2:)-this%xGrid(:n-1)))
#endif

    open(unit=newunit(u),file='xgrid.dat')
    do i_=1,this%nx
       write(u,*) this%xGrid(i_)
    enddo
    close(u)
    open(unit=newunit(u),file='ygrid.dat')
    do i_=1,this%ny
       write(u,*) this%yGrid(i_)
    enddo
    close(u)
  end subroutine create_lattice

  subroutine create_lattice_rectangle(this,n,len_,nf)
    type(Lattice), intent(out) :: this
    integer, dimension(1:2), intent(in) :: n
    real(dl), dimension(1:2), intent(in) :: len_
    integer, intent(in) :: nf

    integer :: u
    integer :: i_
    
    this%time = 0._dl
    this%nx = n(1); this%ny = n(2)
    this%lSize(1:2) = len_(1:2)
    this%nFld = nf

    allocate( this%psi(1:n(1),1:n(2),1:2,1:nf) )
    allocate( this%v_trap(1:this%nx,1:this%ny) )
    allocate( this%xGrid(1:this%nx), this%yGrid(1:this%ny) )

#if defined(PERIODIC)
    this%dx(1:2) = len_(1:2)/dble(n(1:2)); this%dk(1:2) = twopi/len_(1:2)
    call initialize_transform_2d(this%tPair, (/ this%nx,this%ny /) )
    this%xGrid = [ ( (i_-1)*this%dx(1), i_=1,this%nx ) ] - 0.5_dl*this%lSize(1)
    this%yGrid = [ ( (i_-1)*this%dx(2), i_=1,this%ny ) ] - 0.5_dl*this%lSize(2)
#elif defined(INFINITE)
    call initialize_transform_cheby_2d(this%tPair, (/this%nx,this%ny/), 1)
    call transform_rational_cheby(this%tPair, (/ len_(1), len_(2) /))
    this%xGrid = this%tPair%xGrid
    this%yGrid = this%tPair%yGrid
    this%dx(1) = minval(abs(this%xGrid(2:)-this%xGrid(:this%nx-1)))
    this%dx(2) = minval(abs(this%yGrid(2:)-this%yGrid(:this%ny-1)))
#endif

    open(unit=newunit(u),file='xgrid.dat')
    do i_=1,this%nx
       write(u,*) this%xGrid(i_)
    enddo
    close(u)
    open(unit=99,file='ygrid.dat')
    do i_=1,this%ny
       write(99,*) this%yGrid(i_)
    enddo
    close(99)
  end subroutine create_lattice_rectangle
  
  subroutine destroy_lattice(this)
    type(Lattice), intent(inout) :: this

    this%time = -1._dl; this%nx = -1; this%ny = -1
    this%dx = -1.; this%dk = -1.; this%lSize = -1.
    this%nFld = -1
    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%v_trap)) deallocate(this%v_trap)
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
    if (.not.o) open(unit=fNum,file='fields.bin',access='stream',status='replace')
    write(fNum) this%psi(XIND,:,:)
    
  end subroutine write_lattice_data

  real(dl) function field_norm(this) result(amp)
    type(Lattice), intent(inout) :: this

    real(dl), dimension(1:this%nx) :: rho2
    integer :: l
    integer :: j

    amp = 0._dl
    this%tPair%realSpace = 0._dl
    do l=1,this%nfld
       this%tPair%realSpace = this%tPair%realSpace + this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2
    enddo

#if defined(PERIODIC)
    amp = this%dx(1)*this%dx(2)*sum(this%tPair%realSpace(XIND))
#elif defined(INFINITE)
    amp = quadrature_cheby_2d(this%tPair)
#endif
  end function field_norm
  
end module Simulation
