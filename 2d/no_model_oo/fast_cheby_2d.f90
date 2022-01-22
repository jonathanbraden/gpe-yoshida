! TO DO:
!  - Update to allow nonsquare grids
!  - Test and debug for infinite grids
!  - Include mapping parameters (so we aren't necessarily on the [-1:1] interval
!  - Understand the blowup at the two endpoints
!  - Add calculation of minimal dx

module Fast_Cheby_2D
  use, intrinsic :: iso_c_binding
#ifdef USEOMP
  use omp_lib
#endif
  use constants
  implicit none
  include 'fftw3.f03'

  type chebyshevPair2D
     integer :: nx, ny, nx_logical, ny_logical, order
     integer :: nvol
     integer, dimension(1:2) :: type
     real(C_DOUBLE), dimension(:), allocatable :: xGrid, dTheta_dx, d2Theta_dx2
     real(C_DOUBLE), dimension(:), allocatable :: yGrid, dTheta_dy, d2Theta_dy2
     real(C_DOUBLE), dimension(:), allocatable :: Theta_x, Theta_y
     real(C_DOUBLE), pointer :: realSpace(:,:), specSpace(:,:)
     type(C_PTR) :: plan_cos_f, plan_cos_b, plan_sin_b_x1, plan_sin_b_x2
     type(C_PTR), private :: rPtr, sPtr
  end type chebyshevPair2D
  
contains

  subroutine initialize_transform_cheby_2d(this,n,type)
    type(chebyshevPair2D), intent(out) :: this
    integer, dimension(1:2), intent(in) :: n
    integer, intent(in) :: type

    integer :: i
    
    this%type(:) = type
    this%nx = n(1); this%ny = n(2)

    if (allocated(this%xGrid)) deallocate(this%xGrid)
    if (allocated(this%dTheta_dx)) deallocate(this%dTheta_dx)
    if (allocated(this%d2Theta_dx2)) deallocate(this%d2Theta_dx2)
    if (allocated(this%yGrid)) deallocate(this%yGrid)
    if (allocated(this%dTheta_dy)) deallocate(this%dTheta_dy)
    if (allocated(this%d2Theta_dy2)) deallocate(this%d2Theta_dy2)

    allocate(this%xGrid(1:this%nx), this%yGrid(1:this%ny))
    allocate(this%dTheta_dx(1:this%nx), this%dTheta_dy(1:this%ny))
    allocate(this%d2Theta_dx2(1:this%nx), this%d2Theta_dy2(1:this%ny))

    select case (type)
    case (1)  ! Gauss abscissa
       this%nx_logical = 2*this%nx; this%ny_logical = 2*this%ny
       this%type(:) = 1
    case(2)   ! Gauss-Lobatto abscissa
       this%nx_logical = 2*(this%nx+1); this%ny_logical = 2*(this%ny+1)
       this%type(:) = 2
    case default
       print*,"Need to implement default case"
       !!!! Write this
    end select

    this%nvol = this%nx_logical*this%ny_logical
    
    select case (type)
    case (1)  ! Gauss abscissa
       this%Theta_x = (/ (0.5_dl*twopi*(i-0.5_dl)/dble(this%nx), i=1,this%nx) /)
       this%xGrid = cos(this%Theta_x)
       this%dTheta_dx = -1._dl/sqrt(1._dl-this%xGrid**2) 
       this%d2Theta_dx2 = this%xGrid*this%dTheta_dx**3

       this%Theta_y = (/ (0.5_dl*twopi*(i-0.5_dl)/dble(this%ny), i=1,this%ny) /)
       this%yGrid = cos(this%Theta_y)
       this%dTheta_dy = -1._dl/sqrt(1._dl-this%yGrid**2)
       this%d2Theta_dy2 = this%yGrid*this%dTheta_dy**3
       
    case(2)   ! Gauss-Lobatto abscissa
       print*,"Need to implement Lobatto grid"
    case default
       print*,"Need to implement default case"
       !!!! Write this
    end select

    call allocate_2d_array_cheby(n(1), n(2), this%realSpace, this%specSpace, this%rPtr, this%sPtr)
    this%plan_cos_f = fftw_plan_r2r_2d( n(2), n(1), this%realSpace, this%specSpace, FFTW_REDFT10, FFTW_REDFT10, FFTW_MEASURE)
    this%plan_cos_b = fftw_plan_r2r_2d( n(2), n(1), this%specSpace, this%realSpace, FFTW_REDFT01, FFTW_REDFT01, FFTW_MEASURE)
    ! Note reverse array orderings in specifying the type of transform
    this%plan_sin_b_x1 = fftw_plan_r2r_2d( n(2), n(1), this%specSpace, this%realSpace, FFTW_REDFT01, FFTW_RODFT01, FFTW_MEASURE)
    this%plan_sin_b_x2 = fftw_plan_r2r_2d( n(2), n(1), this%specSpace, this%realSpace, FFTW_RODFT01, FFTW_REDFT01, FFTW_MEASURE)
  end subroutine initialize_transform_cheby_2d

  subroutine allocate_2d_array_cheby(n1, n2, arr, ck, fptr, sptr)
    integer, intent(in) :: n1, n2
    real(C_DOUBLE), pointer :: arr(:,:)
    real(C_DOUBLE), pointer :: ck(:,:)
    type(C_PTR) :: fptr, sptr

    fptr = fftw_alloc_real(int(n1*n2, C_SIZE_T))
    call c_f_pointer(fptr, arr, [n1,n2])
    sptr = fftw_alloc_real(int(n1*n2, C_SIZE_T))
    call c_f_pointer(sptr, ck, [n1,n2])
  end subroutine allocate_2d_array_cheby
  
  subroutine destroy_transform_cheby_2d(this)
    type(chebyshevPair2D), intent(inout) :: this

    call fftw_destroy_plan(this%plan_cos_f); call fftw_destroy_plan(this%plan_cos_b)
    call fftw_destroy_plan(this%plan_sin_b_x1); call fftw_destroy_plan(this%plan_sin_b_x2)
    call fftw_free(this%rPtr); call fftw_free(this%sPtr)
    
    if (allocated(this%xGrid)) deallocate(this%xGrid)
    if (allocated(this%dTheta_dx)) deallocate(this%dTheta_dx)
    if (allocated(this%d2Theta_dx2)) deallocate(this%d2Theta_dx2)
    if (allocated(this%yGrid)) deallocate(this%yGrid)
    if (allocated(this%dTheta_dy)) deallocate(this%dTheta_dy)
    if (allocated(this%d2Theta_dy2)) deallocate(this%d2Theta_dy2)
  end subroutine destroy_transform_cheby_2d

  ! TO DO : Fix this for the Lobatto endpoints
  subroutine transform_rational_cheby(this, lPar)
    type(chebyshevPair2D), intent(inout) :: this
    real(dl), dimension(1:2), intent(in) :: lPar

    ! Check these, and include appropriate L's
    this%dTheta_dx = (this%xGrid**2-1._dl) / lPar(1)
    this%d2Theta_dx2 = 2._dl*this%xGrid*(1._dl-this%xGrid**2)**1.5 / lPar(1)**2
    this%xGrid = lPar(1)*this%xGrid / sqrt(1._dl-this%xGrid**2)

    this%dTheta_dy = -(1._dl-this%yGrid**2) / lPar(2)
    this%d2Theta_dy2 = 2._dl*this%yGrid*(1._dl-this%yGrid**2)**1.5 / lPar(2)**2
    this%yGrid = lPar(2)*this%yGrid / sqrt(1._dl-this%yGrid**2)
  end subroutine transform_rational_cheby
  
  subroutine forward_transform_cheby_2d(this)
    type(chebyshevPair2D), intent(inout) :: this
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
  end subroutine forward_transform_cheby_2d

  subroutine backward_transform_cheby_2d(this)
    type(chebyshevPair2D), intent(inout) :: this
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    this%realSpace = this%realSpace / dble(this%nvol)
  end subroutine backward_transform_cheby_2d

  ! Write gradient chain rule here
  ! This is now minimally tested
  subroutine derivative_cheby_2d_chain_x(this)
    type(chebyshevPair2D), intent(inout) :: this

    integer :: i
    real(dl), dimension(1:this%nx,1:this%ny) :: ck
    
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    ck = this%specSpace
    do i=1,this%nx-1
       this%specSpace(i,:) = -dble(i)*ck(i+1,:)/dble(this%nvol)
    enddo
    this%specSpace(this%nx,:) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b_x1, this%specSpace, this%realSpace)
    do i=1,this%ny
       this%realSpace(:,i) = -this%realSpace(:,i) / sqrt(1._dl-this%xGrid**2)
    enddo
  end subroutine derivative_cheby_2d_chain_x
  
  ! Write gradient recurrence here

  ! This needs to be properly debugged
  subroutine laplacian_cheby_2d_recurrence(this)
    type(chebyshevPair2D), intent(inout) :: this

    integer :: i
    real(dl), dimension(1:this%nx,1:this%ny) :: spec_coeff, d2f
    real(dl), dimension(1:this%nx,1:this%ny,1:2) :: der_k
    
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    spec_coeff = this%specSpace
    ! To Do : Zero out coefficient below some threshold

    ! Start with x-direction
    der_k(this%nx,:,1) = 0._dl
    der_k(this%nx-1,:,1) = 2._dl*dble(this%nx-1)*spec_coeff(this%nx,:)
    do i=this%nx-2,1,-1
       der_k(i,:,1) = der_k(i+2,:,1) + 2._dl*dble(i)*spec_coeff(i+1,:)
    enddo

    this%specSpace(this%nx,:) = 0._dl
    this%specSpace(this%nx-1,:) = 2._dl*dble(this%nx-1)*der_k(this%nx,:,1) ! Should be zero
    do i=this%nx-2,1,-1
       this%specSpace(i,:) = this%specSpace(i+2,:) + 2._dl*dble(i)*der_k(i+1,:,1)
    enddo
    der_k(:,:,1) = this%specSpace(:,:)

    ! Derivative along y-direction
    der_k(:,this%ny,2) = 0._dl
    der_k(:,this%ny-1,2) = 2._dl*dble(this%ny-1)*spec_coeff(:,this%ny)
    do i=this%ny-2,1,-1
       der_k(:,i,2) = der_k(:,i+2,2) + 2._dl*dble(i)*spec_coeff(:,i+1)
    enddo
    
    this%specSpace(:,this%ny) = 0._dl
    this%specSpace(:,this%ny-1) = 2._dl*dble(this%ny-1)*der_k(:,this%ny,1)
    do i=this%ny-2,1,-1
       this%specSpace(:,i) = this%specSpace(:,i+2) + 2._dl*dble(i)*der_k(:,i+1,2)
    enddo

    this%specSpace = this%specSpace + der_k(:,:,1)
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    this%realSpace = this%realSpace / dble(this%nx_logical*this%ny_logical)
  end subroutine laplacian_cheby_2d_recurrence

  ! This is all wrong an a huge mess currently
  !
  ! TO DO : Need to do some performance testing when the vectorization isn't on continuous
  !  terms in memory.  I should also test doing a transpose and doing a double loop.
  !  If I want to vectorize the double loop, I could also prestore the i-index values (in 1D) and multiply

  ! TO DO : Check ordering of x vs y indices when defining the DCT vs DST in FFTW
  !   This determines this%plan_sin_b_x1 vs this%plan_sin_b_x2
  subroutine laplacian_cheby_2d_chain(this)
    type(chebyshevPair2D), intent(inout) :: this

    integer :: i, j
    real(dl), dimension(1:this%nx,1:this%ny) :: ck, d2f
    
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    ck = this%specSpace / dble(this%nvol)

    ! This does some ugly cache thrashing
    ! Test if a double loop is faster
    do i=1,this%nx-1
       this%specSpace(i,:) = -dble(i)*ck(i+1,:)
    enddo
    this%specSpace(this%nx,:) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b_x1, this%specSpace, this%realSpace)
    do i=1,this%ny
       d2f(:,i) = -this%realSpace(:,i)*this%xGrid/sqrt(1._dl-this%xGrid**2)**3 ! -cos(this%Theta_x)/sin(this%Theta)**3 = -cos(this%Theta_x)/sqrt(1._dl-cos(this%Theta_x)**2)**3
    enddo
    
    do i=1,this%nx
       this%specSpace(i,:) = -dble(i-1)**2*ck(i,:)
    enddo
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    do i=1,this%ny
       d2f(:,i) = d2f(:,i) + this%realSpace(:,i) / (1._dl-this%xGrid**2) ! *(dTheta_x)**2
    enddo

    ! Debug these lines
    do i=1,this%ny-1
       this%specSpace(:,i) = -dble(i)*ck(:,i+1)
    enddo
    this%specSpace(:,this%ny) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b_x2, this%specSpace, this%realSpace)
    do i=1,this%nx
       d2f(i,:) = d2f(i,:) - this%realSpace(i,:) * this%yGrid(:) / sqrt(1._dl-this%yGrid(:)**2)**3 !*(d2Theta_y)
    enddo
    
    do i=1,this%ny
       this%specSpace(:,i) = -dble(i-1)**2*ck(:,i)
    enddo
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    do i=1,this%nx
       d2f(i,:) = d2f(i,:) + this%realSpace(i,:)  / (1._dl-this%yGrid(:)**2) !*(dTheta_y)**2
    enddo

    this%realSpace = d2f
  end subroutine laplacian_cheby_2d_chain

  subroutine laplacian_cheby_2d_chain_infinite(this)
    type(chebyshevPair2D), intent(inout) :: this

    integer :: i
    real(dl), dimension(1:this%nx,1:this%ny) :: ck, d2f

    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    ck = this%specSpace / dble(this%nvol)

    ! Now write all of this code
  end subroutine laplacian_cheby_2d_chain_infinite
  
  subroutine laplacian_cheby_2d_chain_mapped(this)
    type(chebyshevPair2D), intent(inout) :: this

    integer :: i
    real(dl), dimension(1:this%nx,1:this%ny) :: ck, d2f
    
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    ck = this%specSpace / dble(this%nvol)

    do i=1,this%nx-1
       this%specSpace(i,:) = -dble(i)*ck(i+1,:)
    enddo
    this%specSpace(this%nx,:) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b_x1, this%specSpace, this%realSpace)
    do i=1,this%ny
       this%realSpace(:,i) = this%realSpace(:,i) * this%d2Theta_dx2
    enddo
    d2f = this%realSpace
    
    do i=1,this%nx
       this%specSpace(i,:) = -dble(i-1)**2*ck(i,:)
    enddo
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    do i=1,this%ny
       this%realSpace(:,i) =  this%realSpace(:,i) * this%dTheta_dx**2
    enddo
    d2f = d2f + this%realSpace
    
    do i=1,this%ny-1
       this%specSpace(:,i) = -dble(i)*ck(:,i+1)
    enddo
    this%specSpace(:,this%ny) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b_x2, this%specSpace, this%realSpace)
    do i=1,this%nx
       this%realSpace(i,:) = this%realSpace(i,:) * this%d2Theta_dy2
    enddo
    d2f = d2f + this%realSpace
    
    do i=1,this%ny
       this%specSpace(:,i) = -dble(i-1)**2*ck(:,i)
    enddo
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    do i=1,this%nx
       this%realSpace(i,:) = this%realSpace(i,:) * this%dTheta_dy**2
    enddo

    this%realSpace = d2f + this%realSpace
  end subroutine laplacian_cheby_2d_chain_mapped
  
  subroutine laplacian_cheby_2d_chain_mapped_(this)
    type(chebyshevPair2D), intent(inout) :: this

    integer :: i
    real(dl), dimension(1:this%nx,1:this%ny) :: ck, d2f

    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    ck = this%specSpace / dble(this%nvol)

    do i=1,this%ny-1
       this%specSpace(:,i) = -dble(i)*ck(:,i+1)
    enddo
    this%specSpace(:,this%ny) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b_x2, this%specSpace, this%realSpace)
    do i=1,this%nx
       this%realSpace(i,:) = this%realSpace(i,:) * this%d2theta_dy2 ! Check if this works better than directly saving to d2f
    enddo
    d2f = this%realSpace
    
    do i=1,this%ny
       this%specSpace(:,i) = -dble(i-1)**2*ck(:,i)
    enddo
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    do i=1,this%nx
       this%realSpace(i,:) = this%realSpace(i,:) * this%dtheta_dy**2  ! fix this
    enddo
    d2f = d2f + this%realSpace

    do i=1,this%nx-1
       this%specSpace(i,:) = -dble(i)*ck(i+1,:)
    enddo
    this%specSpace(this%nx,:) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b_x1, this%specSpace, this%realSpace)
    do i=1,this%ny
       this%realSpace(:,i) = this%realSpace(:,i) * this%d2theta_dx2
    enddo
    d2f = d2f + this%realSpace
    
    do i=1,this%nx
       this%specSpace(i,:) = -dble(i-1)**2*ck(i,:)
    enddo
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    do i=1,this%ny
       this%realSpace(:,i) = this%realSpace(:,i)*this%dTheta_dx**2
    enddo

    this%realSpace = this%realSpace + d2f
  end subroutine laplacian_cheby_2d_chain_mapped_
  
end module Fast_Cheby_2D
