module Model_Params
  use constants, only : dl
  use Simulation
  
  type Model
     real(dl), dimension(:,:), allocatable :: g_cross, nu
     real(dl), dimension(:), allocatable :: g_self
     real(dl) :: delta, omega
     integer :: nfld
     integer :: type
  end type Model

  real(dl), dimension(:), allocatable :: v_trap
  real(dl), dimension(:,:), allocatable :: g_cross, nu
  real(dl), dimension(:), allocatable :: g_self
  real(dl) :: mu
  real(dl) :: delta, omega
  integer :: nfld
  integer :: type

contains

  subroutine set_model_parameters(g_,gc_,nu_,mu_,nf)
    real(dl), intent(in) :: g_, gc_, nu_, mu_
    integer, intent(in) :: nf

    integer :: i_,j_
    
    if (allocated(g_self))  deallocate(g_self)
    if (allocated(g_cross)) deallocate(g_cross)
    if (allocated(nu))      deallocate(nu)
    
    allocate(g_self(1:nf))
    allocate( g_cross(1:nf,1:nf), nu(1:nf,1:nf) )

    g_self = g_
    g_cross = gc_
    nu = nu_
    mu = mu_
    
    do i_=1,nf
       g_cross(i_,i_) = 0._dl
       nu(i_,i_) = 0._dl
    enddo
    
    t_loc = 0.
  end subroutine set_model_parameters

  subroutine initialize_trap_potential(this, amp, type)
    type(Lattice), intent(inout) :: this
    real(dl), dimension(1:this%nlat) :: xGrid
    real(dl), intent(in) :: amp
    integer, intent(in), optional :: type

    integer :: i
    integer :: type_

    type_ = 1; if (present(type)) type_ = type
    allocate( v_trap(1:this%nlat) )

    select case (type_)
    case (1)
       v_trap = 0._dl
    case(2)
       v_trap = -0.5_dl*amp*(amp+1._dl) / cosh(this%xGrid)**2
    case(3)
       do i=1,this%nlat
          v_trap(i) = min(0.5_dl*this%xGrid(i)**2,32.)
       enddo
    end select
       
    open(unit=99,file='trap.dat')
    do i=1,this%nlat
       write(99,*) this%xGrid(i), v_trap(i)
    enddo
    close(99)
  end subroutine initialize_trap_potential
  
end module model_params
