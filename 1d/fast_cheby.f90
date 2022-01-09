module Fast_Cheby
  use, intrinsic :: iso_c_binding
#ifdef USEOMP
  use omp_lib
#endif
  use constants
  implicit none
  include 'fftw3.f03'

  type chebyshevPair1D
     integer :: nx, n_logical, order
     integer :: type
     real(C_DOUBLE), dimension(:), allocatable :: xGrid, dTheta, d2Theta
     real(C_DOUBLE), pointer :: realSpace(:)
     real(C_DOUBLE), pointer :: specSpace(:)
     type(C_PTR) :: plan_cos_f, plan_cos_b, plan_sin_b
     type(C_PTR), private :: rPtr, sPtr
  end type chebyshevPair1D
 
contains

  subroutine boot_openmp(nThread)
    integer, optional, intent(in) :: nThread
    integer :: errorOMP
#ifdef USEOMP
    errorOMP = fftw_init_threads()
    if (errorOMP == 0) then
       print*,"Error initializing OpenMP threading for FFTW"
       stop
    endif
    if (present(nThread)) then
       errorOMP = nThread
    else
       print*,"Defaulting to OMP_NUM_THREADS environment variable"
       errorOMP = omp_get_max_threads()
    endif
    call fftw_plan_with_nthreads(errorOMP)
    print*,"FFTW booted using ",errorOMP," threads"
#endif
  end subroutine boot_openmp

  ! Currently assumes Gauss grid.  Adjust to allow for Lobatto and Radau grids
  ! This requires different choices of the FFTW_REDFT and FFTW_RODFT types below
  subroutine initialize_transform_cheby_1d(this,n,type)
    type(chebyshevPair1D), intent(out) :: this
    integer, intent(in) :: n, type

    integer :: i

    this%nx = n
    
    if (allocated(this%xGrid))   deallocate(this%xGrid)
    if (allocated(this%dTheta))  deallocate(this%dTheta)
    if (allocated(this%d2Theta)) deallocate(this%d2Theta)
    allocate(this%xGrid(1:this%nx), this%dTheta(1:this%nx), this%d2Theta(1:this%nx))
    
    select case (type)
    case (1)  ! Gauss abscissa
       this%n_logical = 2*n
       this%type = 1
       this%order = n-1
    case (2)  ! Gauss-Lobatto abscissa
       this%n_logical = 2*(n+1)  ! Check this
       this%type = 2
       this%order = n-1
    case default
       this%nx = n; this%n_logical = 2*n
       this%type = 1
       this%order = n-1
    end select
    
    select case (type)
    case (1)  ! Gaussian Nodes
       this%xGrid = (/ ( cos(0.5_dl*twopi*(i-0.5_dl)/dble(this%nx)), i=1,this%nx ) /)
       this%dTheta = 1._dl/sqrt(1._dl-this%xGrid**2) ! check this
       this%d2Theta = this%xGrid/(1._dl-this%xGrid**2)**1.5  ! check this

    case (2)  ! Lobatto Nodes
       this%xGrid = (/ ( cos(0.5_dl*twopi*(i-1)/dble(this%order)), i=1,this%nx ) /)
       this%dTheta(1) = 0._dl; this%dTheta(this%nx) = 0._dl
       this%dTheta(2:this%nx-1) = 1./sqrt(1._dl-this%xGrid(2:this%nx-1)**2)
       this%d2Theta = 0._dl  ! Fill this in

    case default
       this%xGrid = (/ ( cos(0.5_dl*twopi*(i-0.5_dl)/dble(this%nx)), i=1,this%nx ) /)
       this%dTheta = 1._dl/sqrt(1._dl-this%xGrid**2)
       this%d2Theta = 0._dl ! Fill this in
       
    end select

    call allocate_1d_array_cheby(n, this%realSpace, this%specSpace, this%rPtr, this%sPtr)
    this%plan_cos_f = fftw_plan_r2r_1d(n, this%realSpace, this%specSpace, FFTW_REDFT10, FFTW_MEASURE)
    this%plan_cos_b = fftw_plan_r2r_1d(n, this%specSpace, this%realSpace, FFTW_REDFT01, FFTW_MEASURE)
    this%plan_sin_b = fftw_plan_r2r_1d(n, this%specSpace, this%realSpace, FFTW_RODFT01, FFTW_MEASURE)
  end subroutine initialize_transform_cheby_1d

  subroutine destroy_transform_cheby_1d(this)
    type(chebyshevPair1D), intent(inout) :: this
    call fftw_destroy_plan(this%plan_cos_f); call fftw_destroy_plan(this%plan_cos_b)
    call fftw_destroy_plan(this%plan_sin_b)
    call fftw_free(this%rPtr); call fftw_free(this%sPtr)

    if (allocated(this%xGrid))   deallocate(this%xGrid)
    if (allocated(this%dTheta))  deallocate(this%dTheta)
    if (allocated(this%d2Theta)) deallocate(this%d2Theta)
  end subroutine destroy_transform_cheby_1d

  subroutine transform_rational_cheby(this,lenPar)
    type(chebyshevPair1D), intent(inout) :: this
    real(dl), intent(in) :: lenPar

    ! Use x-coordinates
    this%dTheta = (this%xGrid**2-1._dl) / lenPar
    this%d2Theta = 2._dl*this%xGrid*(1._dl-this%xGrid**2)**1.5 / lenPar**2
    this%xGrid = lenPar*this%xGrid / sqrt(1._dl-this%xGrid**2)

    ! Use r coordinates
!    this%xGrid = this%xGrid / sqrt(1._dl-this%xGrid**2)
!    this%dTheta = -1._dl / (1._dl+this%xGrid**2) / lenPar
!    this%d2Theta = 2._dl*this%xGrid*this%dTheta**2
!    this%xGrid = lenPar*this%xGrid
    
    ! Add a correction on Lobatto grid so endpoints are exactly zero
  end subroutine transform_rational_cheby
  
  subroutine allocate_1d_array_cheby(L, arr, ck, fptr, sptr)
    integer, intent(in) :: L
    real(C_DOUBLE), pointer :: arr(:)
    real(C_DOUBLE), pointer :: ck(:)
    type(C_PTR) :: fptr, sptr

    integer :: n_real, n_spec

    n_real = L; n_spec = L  ! Fix this
    fptr = fftw_alloc_real(int(L, C_SIZE_T)); call c_f_pointer(fptr, arr, [n_real])
    sptr = fftw_alloc_real(int(L, C_SIZE_T)); call c_f_pointer(sptr, ck, [n_spec])
  end subroutine allocate_1d_array_cheby
  
  subroutine forward_transform_cheby_1d(this)
    type(chebyshevPair1D), intent(inout) :: this
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
  end subroutine forward_transform_cheby_1d

  subroutine backward_transform_cheby_1d(this)
    type(chebyshevPair1D), intent(inout) :: this
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    this%realSpace = this%realSpace / dble(this%n_logical)
  end subroutine backward_transform_cheby_1d

  ! TO DO : Need to include basis function normalization in here !!!!
  subroutine derivative_cheby_1d_recurrence(this)
    type(chebyshevPair1D), intent(inout) :: this

    integer :: i
    real(dl), dimension(1:this%nx) :: ck  ! Temporary for debugging.  Remove for performance
!    real(dl) :: max_coeff
    
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)

    ! Zero out machine precision limited modes.  Fix to compare to maxval
!    max_coeff = maxval(abs(this%specSpace))
!    where ( abs(this%specSpace)/max_coeff < 1.e-14)
!       this%specSpace = 0._dl
!    end where
    
    ck = this%specSpace  ! Figure out how to remove this
    
    this%specSpace(this%nx-1) = 2._dl*(this%nx-1._dl)*this%specSpace(this%nx)
    this%specSpace(this%nx) = 0._dl
    do i=this%nx-2, 1, -1
       this%specSpace(i) = this%specSpace(i+2) + 2._dl*dble(i)*ck(i+1)
    enddo
    ! Include relevant dk and normalization somewhere
    ! Flip sign if positive increasing order of x
    
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    this%realSpace = this%realSpace / dble(this%n_logical)
    
    ! Finally, if I transform to the infinite grid, I also need to multiply by a local factor
  end subroutine derivative_cheby_1d_recurrence

  subroutine derivative_cheby_1d_chain(this)
    type(chebyshevPair1D), intent(inout) :: this

    integer :: i

    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    
    do i=1,this%nx-1
       this%specSpace(i) = -dble(i)*this%specSpace(i+1)/dble(this%n_logical) ! Try vectorizing this, check indexing
    enddo
    this%specSpace(this%nx) = 0._dl  ! Zero out the highest mode which isn't in expansion
    call fftw_execute_r2r(this%plan_sin_b, this%specSpace, this%realSpace)

    this%realSpace = -this%realSpace / sqrt(1._dl - this%xGrid**2)  ! Fix with a where statement, add the theta_prime to type
    ! Add boundary conditions if necessary
    !if (type == LOBATTO_TYPE) then
    !   this%realSpace(1) = 0._dl        ! write this
    !   this%realSpace(this%nx) = 0._dl  ! write this
    !endif
  end subroutine derivative_cheby_1d_chain

  subroutine derivative_cheby_1d_infinite(this)
    type(chebyshevPair1D), intent(inout) :: this

    integer :: i

    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    do i=1,this%nx-1
       this%specSpace(i) = -dble(i)*this%specSpace(i+1)/dble(this%n_logical)
    enddo
    this%specSpace(this%nx) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b, this%specSpace, this%realSpace)

    ! Fix this for different L's
    !this%realSpace = -this%realSpace / (1._dl + this%xGrid**2)
    this%realSpace = this%realSpace*this%dTheta
  end subroutine derivative_cheby_1d_infinite


  subroutine derivative_cheby_1d_mapped(this)
    type(chebyshevPair1D), intent(inout) :: this

    integer :: i

    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    do i=1,this%nx-1
       this%specSpace(i) = -dble(i)*this%specSpace(i+1)/dble(this%n_logical)
    enddo
    this%specSpace(this%nx) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b, this%specSpace, this%realSpace)

    ! Fix this for different L's
    !this%realSpace = -this%realSpace / (1._dl + this%xGrid**2)
    this%realSpace = this%realSpace*this%dTheta
  end subroutine derivative_cheby_1d_mapped
  
  subroutine laplacian_cheby_1d_recurrence(this)
    type(chebyshevPair1D), intent(inout) :: this

    integer :: i
    real(dl), dimension(1:this%nx) :: spec_coeff, der_k
    
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    spec_coeff = this%specSpace
    ! To Do: Zero out coefficients below some threshold
    
    der_k(this%nx) = 0._dl
    der_k(this%nx-1) = 2._dl*dble(this%nx-1)*spec_coeff(this%nx)
    do i=this%nx-2,1,-1
       der_k(i) = der_k(i+2) + 2._dl*dble(i)*spec_coeff(i+1)
    enddo

    this%specSpace(this%nx) = 0._dl
    this%specSpace(this%nx-1) = 2._dl*dble(this%nx-1)*der_k(this%nx)
    do i=this%nx-2,1,-1
       this%specSpace(i) = this%specSpace(i+2) + 2._dl*dble(i)*der_k(i+1)
    enddo

    this%specSpace = this%specSpace / dble(this%n_logical)
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
  end subroutine laplacian_cheby_1d_recurrence

  ! Unmerged operations are working, but haven't checked about super high order convergence
  subroutine laplacian_cheby_1d_chain(this)
    type(chebyshevPair1D), intent(inout) :: this
    
    integer :: i
    real(dl), dimension(1:this%nx) :: ck, d2f
    
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    ck = this%specSpace / dble(this%n_logical)

    do i=1,this%nx-1
       this%specSpace(i) = -dble(i)*ck(i+1)
    enddo
    this%specSpace(this%nx) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b, this%specSpace, this%realSpace)
    d2f = -this%realSpace*this%xGrid/sqrt(1._dl-this%xGrid**2)**3

    do i=1,this%nx
       this%specSpace(i) = -dble(i-1)**2*ck(i)
    enddo
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)

    this%realSpace = this%realSpace / (1._dl-this%xGrid**2)
    this%realSpace = this%realSpace + d2f

    ! Should combine some of the above operations into one
    ! this%realSpace = (d2f*this%xGrid/sqrt(1._dl-this%sGrid**2) + this%realSpace)/(1._dl-this%xGrid**2)
  end subroutine laplacian_cheby_1d_chain

  subroutine laplacian_cheby_1d_mapped(this)
    type(chebyshevPair1D), intent(inout) :: this

    integer :: i
    real(dl), dimension(1:this%nx) :: ck, d2f
    
    call fftw_execute_r2r(this%plan_cos_f, this%realSpace, this%specSpace)
    ck = this%specSpace / dble(this%n_logical)

    do i=1,this%nx-1
       this%specSpace(i) = -dble(i)*ck(i+1)
    enddo
    this%specSpace(this%nx) = 0._dl
    call fftw_execute_r2r(this%plan_sin_b, this%specSpace, this%realSpace)
    d2f = this%realSpace * this%d2Theta ! Check this one

    do i=1,this%nx
       this%specSpace(i) = -dble(i-1)**2*ck(i)
    enddo
    call fftw_execute_r2r(this%plan_cos_b, this%specSpace, this%realSpace)
    this%realSpace = this%realSpace*this%dTheta**2 + d2f
  end subroutine laplacian_cheby_1d_mapped
  
end module Fast_Cheby
