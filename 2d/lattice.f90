#include "macros.h"

module Simulation
  use constants, only : dl, twopi
#if defined(PERIODIC)
  use fftw3
#elif defined(INFINITE)
  use Fast_Cheby_2D
#endif
  use Model_Params
  implicit none
  
  type Lattice
     real(dl), dimension(:,:,:,:), allocatable :: psi
     real(dl), dimension(:,:), allocatable :: v_trap
     real(dl), dimension(:), allocatable :: xGrid, yGrid
     real(dl), dimension(1:2) :: dx, lSize, dk
     real(dl) :: time
     integer :: nx, ny
     integer :: nFld
#if defined(PERIODIC)
     type(transformPair2D) :: tPair
#elif defined(INFINITE)
     type(chebyshevPair2D) :: tPair
#endif
     type(Model) :: model
  end type Lattice
  
contains

  subroutine create_lattice(this,n,len,nf)
    type(Lattice), intent(out) :: this
    integer, intent(in) :: n, nf
    real(dl), intent(in) :: len

    integer :: i_
    
    this%time = 0._dl
    this%nx = n; this%ny = n
    this%lSize(1:2) = len
    this%nFld = nf

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

    open(unit=99,file='xgrid.dat')
    do i_=1,this%nx
       write(99,*) this%xGrid(i_)
    enddo
    close(99)
    open(unit=99,file='ygrid.dat')
    do i_=1,this%ny
       write(99,*) this%yGrid(i_)
    enddo
    close(99)
  end subroutine create_lattice

  subroutine create_lattice_rectangle(this,n,len_,nf)
    type(Lattice), intent(out) :: this
    integer, dimension(1:2), intent(in) :: n
    real(dl), dimension(1:2), intent(in) :: len_
    integer, intent(in) :: nf

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

    open(unit=99,file='xgrid.dat')
    do i_=1,this%nx
       write(99,*) this%xGrid(i_)
    enddo
    close(99)
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
    if (.not.o) open(unit=fNum,file='fields.bin',access='stream')

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
  
  ! TO DO: Reduce memory footprint here and benchmark vector vs parallel
  real(dl) function chemical_potential(this) result(chemP)
    type(Lattice), intent(inout) :: this

    real(dl), dimension(1:this%nx,1:this%ny) :: rho2, mu_loc
    real(dl) :: g_loc
    real(dl) :: norm
    integer :: i,j, l
    
    chemP = 0._dl
    mu_loc = 0._dl

    do l=1,this%nfld
       g_loc = g_self(l)
       rho2 = this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2

       this%tPair%realSpace = this%psi(XIND,1,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,1,l)

       this%tPair%realSpace = this%psi(XIND,2,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,2,l)

       mu_loc = mu_loc + this%v_trap*rho2 + g_loc*rho2**2
    enddo

    norm = field_norm(this)
#if defined(PERIODIC)
    chemP = this%dx(1)*this%dx(2)*sum(mu_loc)
#elif defined(INFINITE)
    this%tPair%realSpace = mu_loc
    chemP = quadrature_cheby_2d(this%tPair)
#endif
    chemP = chemP / norm
  end function chemical_potential

  ! TO DO : To improve performance, I should really add some accumulators
  !         for various energy components in the loop.
  !         I can also accumulate the field normalization in here.
  real(dl) function chemical_potential_full(this) result(chemP)
    type(Lattice), intent(inout) :: this

    real(dl), dimension(1:this%nx,1:this%ny) :: rho2, mu_loc
    real(dl), dimension(1:this%nfld) :: nu_loc
    real(dl) :: g_loc
    real(dl) :: psi_norm
    integer :: l,m

    chemP = 0._dl
    mu_loc = 0._dl
    
    do l=1,this%nfld
       g_loc = g_self(l)
       nu_loc = nu(:,l); nu_loc(l) = 0._dl
       rho2 = this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2

       this%tPair%realSpace = this%psi(XIND,1,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,1,l)

       this%tPair%realSpace = this%psi(XIND,2,l)
#if defined(PERIODIC)
       call laplacian_2d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
       call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
       mu_loc = mu_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,2,l)
       
       mu_loc = mu_loc + this%v_trap*rho2 + g_loc*rho2**2

       do m=1,this%nfld
          mu_loc = mu_loc - nu_loc(m) * ( this%psi(XIND,1,l)*this%psi(XIND,1,m) + this%psi(XIND,2,l)*this%psi(XIND,2,m) )
       enddo
    enddo
! Need to include the field normalization here

    psi_norm = field_norm(this)
#if defined(PERIODIC)
    chemP = this%dx(1)*this%dx(2)*sum(mu_loc)
#elif defined(INFINITE)
    print*,"Chemical potential calculation not tested on Chebyshev"
    this%tPair%realSpace = mu_loc
    chemP = quadrature_cheby_2d(this%tPair)
#endif
    chemP = chemP / psi_norm
  end function chemical_potential_full
  
  real(dl) function energy(this) result(en)
    type(Lattice), intent(inout) :: this

    real(dl), dimension(1:this%nx,1:this%ny) :: rho2, en_loc
    integer :: i, l
    real(dl) :: g_loc

    en = 0._dl; en_loc = 0._dl

    do l=1,this%nfld
       g_loc = g_self(l)
       rho2 = this%psi(XIND,1,l)**2 + this%psi(XIND,2,l)**2

       do i=1,2
          this%tPair%realSpace = this%psi(XIND,i,l)
#if defined(PERIODIC)
          call laplacian_2d_wtype(this%tPair,this%dk)
#elif defined(INFINITE)
          call laplacian_cheby_2d_chain_mapped(this%tPair)
#endif
          en_loc = en_loc - 0.5_dl*this%tPair%realSpace*this%psi(XIND,i,l)
       enddo
       en_loc = en_loc + this%v_trap*rho2 + 0.5_dl*g_loc*rho2**2
    enddo
#if defined(PERIODIC)
    en = sum(en_loc)*this%dx(1)*this%dx(2)
#elif defined(INFINITE)
    this%tPair%realSpace = en_loc
    en = quadrature_cheby_2d(this%tPair)
#endif
  end function energy

end module Simulation
