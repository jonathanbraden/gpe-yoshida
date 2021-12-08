!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for implementing 1-dimensional derivatives and antiderivatives using pseudospectral methods
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> @author
!> Jonathan Braden, University College London
!>
! DESCRIPTION
!> @brief
!> An implementation of pseudospectral based one-dimensional derivatives and antiderivatives
!>
!> This module is used to compute derivatives and antiderivatives of one-dimensional fields
!> using pseudospectral approximations.
!> The required FFT's are implemented using the open source package FFTW.
!>
!> A typical usage is as follows
!> @code{.f90}
!>  type(transformPair1D) :: tPair
!>  real(C_DOUBLE), parameter :: dk = 1._dl
!>  call initialize_transform_1d(tPair, 1024)
!>  tPair%realSpace(:) = initialField(:)
!>  call derivative_n_1d(tPair,dk)
!> @endcode
!> @warning When calling any of the derivative subroutines, the field to be differentiated must already be stored in tPair%realSpace.
!> @warning The differentiated function is stored in tPair%realSpace after the subroutine call.
!>
!> The FFT convention used in FFTW is (in 1D with N lattice sites)
!> \f[ g_k = \sum_l e^{i2\pi kl / N} f(x_l) \f]
!> so that the nth derivative multiplies the Fourier amplitudes by \f[A^{(n)} = \left(i dk \right)^n \f].
!> As well, FFTW computes an unnormalised transform,
!> \f[ f(x_l) = \sum_k e^{-i2\pi kl /N} g_k \f]
!> so that the final answer must me multiplied by \f$ N_{\mathrm{lat}}^{-1} \f$
!> where \f$ N_{lat} \f$ is the total number of lattice sites.
!> In the actual implementation of the derivative routines below, these are implemented as an normalisation factor \f[ \mathrm{norm} = N_{lat}^{-1}\left( i dk \right)^n \f]
!> multiplying each Fourier amplitude.
!> Additionally, for inverse derivatives (which constitute elliptic equations which require boundary conditions),
!> the boundary conditions are automatically assumed to be periodic since this module is based on Fourier pseudospectral methods.
!> In particular, any mean value for a field we attempt to compute and inverse derivative of will be ignored, since it will not posess continuous solutions
!> within the space of continuous functions.
!> @warning The implicit periodic boundary conditions must be kept in mind when using this module to solve elliptic equations.
!>
! TO DO
!> @todo
!> @arg Include sample use cases in overall module documentation.
!> @arg Implement higher dimensional calculations
!> @arg Include required OpenMP directives
!> @arg Implement MPI communication for parallelisation
!> @arg Use separate flags for FFTW OpenMP parallellism, and parallellism in the loops
!> @arg Allow user to specify level of tuning of FFTW through a flag the user passes in during initialisation
!> @arg Link documentation of derivative functions to a single function where I'll show Fourier conventions, etc.
!>      Or else document this in the start of the module?  Improve the description of the calculation above
!> @arg Since higher dimensions allow for vector derivatives (ie. the gradient), implement a type that allows for this to be stored easily
!>      instead of the ugly rigging that is being used now.  This will introduce additional tuning parameters based on convenience, speed and memory
!> @arg Finish code example in module documentation
!> @arg Either include in this module, or make a separate convenience module, where a transform pair type is defined along with some subroutines that call this module
!>      Then have simple calls like derivative(f, df, dk) without referencing the transform, with f the input field, df the differentiated field and dk the fourier spacing
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fftw3
  use, intrinsic :: iso_c_binding
#ifdef USEOMP
  use omp_lib
#endif
  use constants
  implicit none
  include 'fftw3.f03'

  complex(C_DOUBLE_COMPLEX), parameter :: iImag = (0._dl,1._dl)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DESCRIPTION
  !> @brief
  !> Store information necessary to perform FFT's on a 1D field
  !> @param nx           integer The number of grid points describing the field
  !> @param nnx          integer Nyquist mode for the field, computed from nx
  !> @param realSpace    (array(1:nx), C_DOUBLE) Stores the real space values of the field
  !> @param specSpace    (array(1:nnx), C_DOUBLE_COMPLEX) Stores the Fourier (spectral) amplitudes of the field
  !> @param planf        (C_PTR) pointer to FFTW plan for transforming from real to spectral space
  !> @param planb        (C_PTR) pointer to FFTW plan for transforming from spectral to real space
  !> @param rPtr         (C_PTR) pointer to memory storing the real array values (needed to free memory)
  !> @param sPtr         (C_PTR) pointer to memory storing the spectral space values (needed to free memory)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type transformPair1D
     integer :: nx, nnx
     real(C_DOUBLE), pointer :: realSpace(:)
     complex(C_DOUBLE_COMPLEX), pointer :: specSpace(:)
     type(C_PTR) :: planf, planb
     type(C_PTR), private :: rPtr, sPtr
  end type transformPair1D

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DESCRIPTION
  !> @brief
  !> Store information to perform 2D FFT pseudospec based derivatives
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type transformPair2D
     integer :: nx, ny, nnx, nny
     real(C_DOUBLE), pointer :: realSpace(:,:)
     complex(C_DOUBLE_COMPLEX), pointer :: specSpace(:,:)
     type(C_PTR) :: planf, planb
     type(C_PTR), private :: rPtr, sPtr
  end type transformPair2D

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DESCRIPTION
  !> @brief
  !> Store information to perform 3D FFT pseudospec based derivatives
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type transformPair3D
     integer :: nx, ny, nz, nnx, nny, nnz
     real(C_DOUBLE), pointer :: realSpace(:,:,:)
     complex(C_DOUBLE_COMPLEX), pointer :: specSpace(:,:,:)
     type(C_PTR) :: planf, planb
     type(C_PTR), private :: rPtr, sPtr
  end type transformPair3D

contains
   
!****************************!
! Setup Transformation Pairs !
!****************************!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DESCRIPTION
  !> @brief
  !> Set up OpenMP threading for FFTW.
  !>
  !> Set OpenMP threading with nThreads.
  !> The functionality is set by the preprocessor flag USEOMP during compilation
  !> If USEOMP is not defined at compile time, this subroutine does nothing.
  !>
  !> @param[in] nThread (integer, optional) The number of OpenMP threads to use in FFTW
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

#define RSPACE1D create_transform_1d%realSpace
#define SSPACE1D create_transform_1d%specSPACE

  function create_transform_1d(n)
    type(transformPair1D) :: create_transform_1d
    integer, intent(in) :: n

    call allocate_1d_array(n, RSPACE1D, SSPACE1D, create_transform_1d%rPtr, create_transform_1d%sPtr)
    create_transform_1d%nx = n

    create_transform_1d%planf = fftw_plan_dft_r2c_1d(n, RSPACE1D, SSPACE1D, FFTW_MEASURE)
    create_transform_1d%planb = fftw_plan_dft_c2r_1d(n, SSPACE1D, RSPACE1D, FFTW_MEASURE)
  end function create_transform_1d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DESCRIPTION
  !> @brief
  !> Create a transform pair for doing 1D transforms
  !>
  !> @param[out] this (transformPair1D) The transform pair to initialise
  !> @param[in]  n (integer)            The number of lattice sites of the field
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_transform_1d(this, n)
    type(transformPair1D), intent(out) :: this
    integer, intent(in) :: n

    this%nx = n; this%nnx = n/2 + 1
    call allocate_1d_array(n, this%realSpace, this%specSpace, this%rPtr , this%sPtr)

    this%planf = fftw_plan_dft_r2c_1d(n, this%realSpace, this%specSpace, FFTW_MEASURE)
    this%planb = fftw_plan_dft_c2r_1d(n, this%specSpace, this%realSpace, FFTW_MEASURE)
  end subroutine initialize_transform_1d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DESCRIPTION
  !> @brief
  !> Destroy the 1d transform pair and free the memory
  !>
  !> @param[in,out] this (transformPair1D) The transform to destroy
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine destroy_transform_1d(this)
    type(transformPair1D), intent(inout) :: this

    call fftw_destroy_plan(this%planf); call fftw_destroy_plan(this%planb)
    call fftw_free(this%rPtr); call fftw_free(this%sPtr)
  end subroutine destroy_transform_1d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief
  !> Create a 2d transform pair
  !>
  !> Create the specified 2D scalar transform.
  !> This subroutine will automatically create real and Fourier space arrays
  !> of the correct dimentsions.
  !> It then initialises the forward and backward FFT plans in FFTW.
  !>
  !> @param[out] this (transformPair2D) The 2D transform to initialise
  !> @param[in]  n (integer, dim(1:2))  The number of grid sites for the field
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_transform_2d(this, n)
    type(transformPair2D), intent(out) :: this
    integer, intent(in), dimension(1:2) :: n
    
    this%nx = n(1); this%ny = n(2)
    this%nnx = n(1)/2+1; this%nny = n(2)/2+1

    call allocate_2d_array(n(1),n(2),this%realSpace,this%specSpace,this%rPtr,this%sPtr)

    this%planf = fftw_plan_dft_r2c_2d(n(2),n(1),this%realSpace,this%specSpace,FFTW_MEASURE)
    this%planb = fftw_plan_dft_c2r_2d(n(2),n(1),this%specSpace,this%realSpace,FFTW_MEASURE)
  end subroutine initialize_transform_2d

  !>@brief
  !> Destroy the given 2D transform and free all allocated memory.
  !> After calling this, the transform is unusable
  !>
  !>@param[in,out] this (transformPair2D) The transform to destroy.
  subroutine destroy_transform_2d(this)
    type(transformPair2D), intent(inout) :: this
    
    call fftw_destroy_plan(this%planf); call fftw_destroy_plan(this%planb)
    call fftw_free(this%rPtr); call fftw_free(this%sPtr)
  end subroutine destroy_transform_2d

  !> @brief
  !> Create a 3D transform pair
  !>
  !> @param[out] this (transformPair3D)  The transform to create
  !> @param[in]  n    (integer, dim(1:3) The number of grid points in each dimension
  !>
  !> @todo
  !> @arg finish the creation of the plans
  subroutine initialize_transform_3d(this, n)
    type(transformPair3D), intent(out) :: this
    integer, intent(in), dimension(1:3) :: n

    this%nx = n(1); this%ny = n(2); this%nz = n(3)
    this%nnx = n(1)/2+1; this%nny = n(2)/2+1; this%nnz = n(3)/2+1
    call allocate_3d_array(n(1),n(2),n(3),this%realSpace,this%specSpace)

!      this%planf = fftw_plan_dft_r2c_3d( this%realSpace, this%specSpace, FFTW_MEASURE)
!      this%planb = fftw_plan_dft_c2r_3d( this%specSpace, this%realSpace, FFTW_MEASURE)
  end subroutine initialize_transform_3d

  !> @brief
  !> Destroy a 3D transform and free allocated memory.
  !> The transform is no longer usable after calling this subroutine.
  !>
  !> @param[in,out] this (transformPair3D) The transform to destroy
  subroutine destroy_transform_3d(this)
    type(transformPair3D), intent(inout) :: this

    call fftw_destroy_plan(this%planf); call fftw_destroy_plan(this%planb)
    call fftw_free(this%rPtr); call fftw_free(this%sPtr)
  end subroutine destroy_transform_3d

  !********************!
  ! Helper Subroutines !
  !********************!
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DESCRIPTION
  !> @brief
  !> Allocate the real and spectral space arrays for use with FFTW
  !>
  !> Initialise the 1D arrays required to do one-dimensional FFTs using the FFTW package.
  !> This subroutine automates the process of
  !> @arg determining the size of the Fourier space array for a given real array
  !> @arg allocating SIMD aligned memory
  !> @arg referencing the allocated memory with the designated Fortran pointers
  !>
  !> Memory is allocated in the C_PTRs using fftw_alloc_*, so to deallocate the memory the
  !> same C_PTRs must be used in fftw_free.
  !> This is automatically implemented in destroy_arrays_1d
  !>
  !> @param[in]  L     The number of grid sites
  !> @param[out] arr   (real(C_DOUBLE) array)            Fortran pointer to the real space array
  !> @param[out] Fk    (complex(C_DOUBLE_COMPLEX) array) Fortran pointer to the spectral space array
  !> @param[out] fptr  (C_PTR)                           Pointer to the array storing the real space field
  !> @param[out] fkptr (C_PTR)                           Pointer to the array storing the spectral space field
  !>
  !> @todo
  !> @arg Add intent for the 4 array variables in the argument
  !> @arg Add a link to destroy_arrays_1d in the documentation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine allocate_1d_array(L, arr, Fk, fptr, fkptr)
    integer, intent(in) :: L
    real(C_DOUBLE), pointer :: arr(:)
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
    type(C_PTR) :: fptr, fkptr
    
    integer :: LL
    LL = L/2+1

    fptr = fftw_alloc_real(int(L, C_SIZE_T)); call c_f_pointer(fptr, arr, [L])
    fkptr = fftw_alloc_complex(int(L, C_SIZE_T)); call c_f_pointer(fkptr, Fk, [LL])
  end subroutine allocate_1d_array

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DESCRIPTION
  !> @brief
  !> Destroy the real/spectral space data and the FFTW plans.
  !>
  !> Deallocate the memory previously allocated to the arrays and FFTW plans for a transform pair.
  !> These are conveniently created using allocate_1d_array.
  !> Since SIMD aligned memory is allocated using fftw_alloc_*, the memory should
  !> be deallocated using fftw_free referencing the C_PTRs (not the Fortran pointers to the arrays).
  !> This is automatically done in this subroutine.
  !>
  !> @param[in,out] planf (C_PTR) FFTW plan for forward transform
  !> @param[in,out] planb (C_PTR) FFTW plan for backward transform
  !> @param[in,out] fptr  (C_PTR) Pointer to memory storing real array allocated using fftw_alloc_real
  !> @param[in,out] fkptr (C_PTR) Pointer to memory storing Fourier array allocated using fftw_alloc_complex
  !>
  !> @todo
  !> @arg Add a link to allocate_1d_array in the documentation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine destroy_arrays_1d(planf,planb,fptr,fkptr)
    type(C_PTR) :: planf, planb
    type(C_PTR) :: fptr, fkptr
    call fftw_destroy_plan(planf); call fftw_destroy_plan(planb)
    call fftw_free(fptr); call fftw_free(fkptr)
  end subroutine destroy_arrays_1d

  subroutine allocate_2d_array(L,M, arr, Fk, fptr, fkptr)
    integer :: L,M
    real(C_DOUBLE), pointer :: arr(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:)
    type(C_PTR) :: fptr, fkptr

    integer :: LL

    LL = L/2+1
! Check orderings in here      
    fptr = fftw_alloc_real(int(L*M, C_SIZE_T))
    call c_f_pointer(fptr, arr, [L,M])
    fkptr = fftw_alloc_complex(int(LL*M, C_SIZE_T))
    call c_f_pointer(fkptr, Fk, [LL,M])
  end subroutine allocate_2d_array

  subroutine allocate_3d_array(L,M,N, arr, Fk)
    integer :: L,M,N
    real(C_DOUBLE), pointer :: arr(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
    
    type(C_PTR) :: fptr, fkptr
    integer :: LL

    LL = L/2+1

    fptr = fftw_alloc_real(int(L*M*N, C_SIZE_T)); call c_f_pointer(fptr, arr, [L,M,N])
    fkptr = fftw_alloc_complex(int(LL*M*N, C_SIZE_T)); call c_f_pointer(fkptr, Fk, [LL,M,N])
  end subroutine allocate_3d_array

!********************************************************!
! Mathematical subroutines using the transform pair type !
!********************************************************!

  !> @brief
  !> Compute the inverse derivative under periodic boundary conditions
  !>
  !> Solve \f[ \partial_x\phi = f(x) \f] for the field f(x) stored in tPair.
  !> The boundary conditions are assumed periodic (since this is a Fourier based inversion).
  !> The mean value of \f$\phi\f$ is assumed to be zero.
  !> As well, if \f$ f(x) \f$ has a nonzero mean value, it is ignored since this would result in a
  !> discontinuous solution \f$\phi\f$ in the space of periodic functions.
  !> Since the transform pairs don't know about the side length of the box, the wavenumber of the fundamental mode must also be provided
  !> (and will simply provide a overall normalising factor).
  !>
  !> @param[in,out] this (transformPair1D) The transform pair to use.  tPair%realSpace should hold the field \f$ f(x) \f$.
  !> @param[in] dk (double) Wavenumber of the fundamental mode \f$ dk = 2\pi L^{-1}\f$.  Equivalently spacing of modes in Fourier space.
  !>
  !> @todo
  !> @arg Make sure that having a mean constant derivate doesn't screw everything up.  Why could I get a sawtooth before?
  !> @arg Implement OpenMP parallelisation of the loop.
  !> @arg Test storing the grid numbers as an array and using vectorisation vs the loop approach
  subroutine inverse_derivative_1d_wtype(tPair, dk)
    type(transformPair1D), intent(inout) :: tPair
    real(dl), intent(in) :: dk
    integer :: i
    complex(C_DOUBLE_COMPLEX) :: norm

    norm = -iImag / dk / dble(tPair%nx)
    call fftw_execute_dft_r2c(tPair%planf, tPair%realSpace, tPair%specSpace)
    do i=2,tPair%nnx
       tPair%specSpace(i) = norm*tPair%specSpace(i) / dble(i-1)
    enddo
    tPair%specSpace(1) = 0._dl
    call fftw_execute_dft_c2r(tPair%planb, tPair%specSpace, tPair%realSpace)
  end subroutine inverse_derivative_1d_wtype

  subroutine derivative_1d_wtype(tPair, dk)
    type(transformPair1D), intent(inout) :: tPair
    real(dl), intent(in) :: dk
    integer :: i
    complex(C_DOUBLE_COMPLEX) :: norm

    norm = dk*iImag / dble(tPair%nx)
    call fftw_execute_dft_r2c(tPair%planf, tPair%realSpace, tPair%specSpace)
    do i=1,tPair%nnx
       tPair%specSpace(i) = dble(i-1)*norm*tPair%specSpace(i)
    enddo
    call fftw_execute_dft_c2r(tPair%planb, tPair%specSpace, tPair%realSpace)
  end subroutine derivative_1d_wtype

  subroutine laplacian_1d_wtype(tPair, dk)
    type(transformPair1D), intent(inout) :: tPair
    real(dl), intent(in) :: dk
    integer :: i
    complex(C_DOUBLE_COMPLEX) :: norm

    norm = -dk**2 / dble(tPair%nx)
    call fftw_execute_dft_r2c(tPair%planf, tPair%realSpace, tPair%specSpace)
    do i=1,tPair%nnx
       tPair%specSpace(i) = (i-1)**2*norm * tPair%specSpace(i)
    enddo
    call fftw_execute_dft_c2r(tPair%planb, tPair%specSpace, tPair%realSpace)
  end subroutine laplacian_1d_wtype

  !> @brief
  !> Compute the nth derivative of a field \f$ f(x) \f$.
  !>
  !> Computes the nth derivative of the field \f$ f(x) \f$ stored in tPair%realSpace.
  !> This implementation assumes n is a positive integer.
  !> This subroutine can be used in place of the other explicitly named derivative subroutines.
  !> derivative_1d_wtype and laplacian_1d_wtype
  !> \warning The field whose derivative we want to compute must be previously stored in tPair%realSpace
  !>
  !> @param[in,out] tPair (transformPair1D) The transform pair to use
  !> @param[in]     dk (double) 
  !> 
  !> @todo
  !> @arg Extend to allow negative derivatives
  !> @arg Extend to allow non integer derivatives (ie. pseudo-differential equations)
  !> @arg Link documentation to derivative and laplacian subroutines
  subroutine derivative_n_1d_wtype(tPair, dk, order)
    type(transformPair1D), intent(inout) :: tPair
    real(dl), intent(in) :: dk
    integer, intent(in) :: order
    
    integer :: i
    complex(C_DOUBLE_COMPLEX) :: norm

    norm = ( dk * iImag )**order / dble(tPair%nx)
    call fftw_execute_dft_r2c(tPair%planf, tPair%realSpace, tPair%specSpace)
    do i=1,tPair%nnx
       tPair%specSpace(i) = (i-1)**order*norm *tPair%specSpace(i)
    enddo
    call fftw_execute_dft_c2r(tPair%planb, tPair%specSpace, tPair%realSpace)
  end subroutine derivative_n_1d_wtype

  !> @brief
  !> Compute the squared gradient of a 1D field.
  !>
  !> @warning This will lead to aliasing if the field is not sufficiently resolved.
  !>          For possible solutions, see for example the 2/3 rule
  !>
  !> @param[in,out] tPair (transformPair1D) Transform pair to use.  Field to transform stored in tPair%realSpace
  !> @param[in]     dk    (double) Wavenumber of the fundamental mode
  !>
  !> @todo
  !> @arg Implement 2/3 dealiasing as an option (or a separate subroutine)
  !> @arg Fill in documentation of the subroutine
  subroutine gradsquared_1d_wtype(tPair, dk)
    type(transformPair1D), intent(inout) :: tPair
    real(dl), intent(in) :: dk
    integer :: i
    complex(C_DOUBLE_COMPLEX) :: norm

    norm = iImag*dk / dble(tPair%nx)
    call fftw_execute_dft_r2c(tPair%planf, tPair%realSpace, tPair%specSpace)
    do i=1,tPair%nnx
       tPair%specSpace(i) = (i-1)*norm * tPair%specSpace(i)
    enddo
    call fftw_execute_dft_c2r(tPair%planb, tPair%specSpace, tPair%realSpace)
    
    tPair%realSpace = tPair%realSpace**2
  end subroutine gradsquared_1d_wtype

  subroutine grad_adotb_1d_wtype(n,f1,f2,tPair,dk)
    integer, intent(in) :: n
    real(C_DOUBLE), intent(in) :: f1(1:n),f2(1:n)
    type(transformPair1D), intent(inout) :: tPair
    real(dl), intent(in) :: dk

    real(C_DOUBLE) :: ftmp(1:n)
    integer :: i, nn
    complex(C_DOUBLE_COMPLEX) :: norm

    ! Why is this the right norm?  Why is the denominator squared?
    ! This seems incorrect
    norm = dk*iImag / dble(tPair%nx)**2

    nn = tPair%nnx
    tPair%realSpace = f1
    call fftw_execute_dft_r2c(tPair%planf,tPair%realSpace,tPair%specSpace)
    do i=1,nn
       tPair%specSpace(i) = (i-1)*norm*tPair%specSpace(i)
    enddo
    call fftw_execute_dft_c2r(tPair%planb,tPair%specSpace,tPair%realSpace)
    ftmp=tPair%realSpace
    
    tPair%realSpace = f2
    call fftw_execute_dft_r2c(tPair%planf,tPair%realSpace,tPair%specSpace)
    do i=1,nn
       tPair%specSpace(i) = (i-1)*norm*tPair%specSpace(i)
    enddo
    call fftw_execute_dft_c2r(tPair%planb,tPair%specSpace,tPair%realSpace)
    
    tPair%realSpace = ftmp*tPair%realSpace
  end subroutine grad_adotb_1d_wtype

! To add : allocate and initialize inplace/outofplace FFTW plans
! To add : extend this to work with MPI

!!!!!!!!!!!!!!
! Mathematical Subroutines, compute various derivatives
!!!!!!!!!!!!!!
    subroutine inverse_derivative_1d(n,f,Fk,planf,planb, dk)
      integer, intent(in) :: n
      real(dl), intent(in) :: dk
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      type(C_PTR) :: planf, planb

      integer :: i, nn
      complex(C_DOUBLE_COMPLEX) :: norm

      nn = n/2+1
      norm = -iImag / dk / dble(n)
      call fftw_execute_dft_r2c(planf, f, Fk)
      do i=2,nn
         Fk(i) = norm*Fk(i) / dble(i-1)
      enddo
      Fk(1) = 0._dl
      call fftw_execute_dft_c2r(planb, Fk, f)
    end subroutine inverse_derivative_1d

    subroutine derivative_1d(n,f,Fk,planf,planb,dk)
      integer, intent(in) :: n
      real(dl), intent(in) :: dk
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      type(C_PTR), intent(in) :: planf, planb

      integer :: i,nn

      nn = n/2+1
      call fftw_execute_dft_r2c(planf, f, Fk)
      do i=1,nn
         Fk(i) = (i-1)*dk*iImag*Fk(i)
      enddo
      call fftw_execute_dft_c2r(planb, Fk, f)
      f = f / dble(n)
    end subroutine derivative_1d

    subroutine laplacian_1d(n, f, Fk, planf, planb, dk)
      integer, intent(in) :: n
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      real(dl), intent(in) :: dk
      type(C_PTR), intent(in) :: planf, planb

      integer :: i,nn

      nn = n/2+1
      call fftw_execute_dft_r2c(planf, f, Fk)
      do i=1,nn
         Fk(i) = -((i-1)*dk)**2*Fk(i)
      enddo
      call fftw_execute_dft_c2r(planb, Fk, f)
      f = f / dble(n)
    end subroutine laplacian_1d

    subroutine gradsquared_1d(n,f,Fk,planf,planb,dk)
      integer, intent(in) :: n
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      real(dl), intent(in) :: dk
      type(C_PTR), intent(in) :: planf, planb

      integer :: i,nn
      
      nn = n/2+1
      call fftw_execute_dft_r2c(planf,f,Fk)
      do i=1,nn
         Fk(i) = (i-1)*dk*iImag*Fk(i)
      enddo
      call fftw_execute_dft_c2r(planb,Fk,f)
      f = f**2 / dble(n)**2
    end subroutine gradsquared_1d

    subroutine grad4_1d(n, f, Fk, planf, planb, dk)
      integer, intent(in) :: n
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      real(dl), intent(in) :: dk
      type(C_PTR), intent(in) :: planf, planb

      integer :: i,nn

      nn = n/2+1
      call fftw_execute_dft_r2c(planf, f, Fk)
      do i=1,nn
         Fk(i) = ((i-1)*dk)**4*Fk(i)
      enddo
      call fftw_execute_dft_c2r(planb, Fk, f)
      f = f / dble(n)
    end subroutine grad4_1d

    subroutine grad_adotb_1d(n,f1,f2,f,Fk,planf,planb,dk)
      integer, intent(in) :: n
      real(C_DOUBLE), intent(in) :: f1(1:n),f2(1:n)
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      real(dl), intent(in) :: dk
      type(C_PTR), intent(in) :: planf, planb

      real(C_DOUBLE) :: ftmp(1:n)
      integer :: i,nn

      nn = n/2+1
      f = f1
      call fftw_execute_dft_r2c(planf,f,Fk)
      do i=1,nn
         Fk(i) = (i-1)*dk*iImag*Fk(i)
      enddo
      call fftw_execute_dft_c2r(planb,Fk,f)
      ftmp=f

      f = f2
      call fftw_execute_dft_r2c(planf,f,Fk)
      do i=1,nn
         Fk(i) = (i-1)*dk*iImag*Fk(i)
      enddo
      call fftw_execute_dft_c2r(planb,Fk,f)

      f = ftmp*f / dble(n)**2
    end subroutine grad_adotb_1d

    subroutine laplacian_2d(n1, n2, f, Fk, dk)
      integer :: n1, n2
      real(C_DOUBLE), pointer :: f(:,:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:)
      real(dl) :: dk

      real(dl) :: rad2
      integer :: nn1, nn2
      integer :: i,j,jj

      type(C_PTR) :: planf, planb

      nn1 = n1/2+1; nn2=n2/2+1
      planf = fftw_plan_dft_r2c_2d(n2, n1, f, Fk, FFTW_ESTIMATE)
      planb = fftw_plan_dft_c2r_2d(n2, n1, Fk, f, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(planf, f, Fk)
      do j=1,n2; if (j>nn2) then; jj = n2+1-j; else; jj=j-1; endif
         do i=1,nn1
            rad2 = dble((i-1)**2 + jj**2)
            Fk(i,j) = -rad2*dk**2*Fk(i,j)
         enddo
      enddo

      call fftw_execute_dft_c2r(planb, Fk, f)
      call fftw_destroy_plan(planf)
      call fftw_destroy_plan(planb)
      f = f / dble(n1) / dble(n2)
    end subroutine laplacian_2d

    subroutine laplacian_3d(n1, n2, n3, f, Fk, dk)
      integer :: n1, n2, n3
      real(C_DOUBLE), pointer :: f(:,:,:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
      real(dl) :: dk

      real(dl) :: rad2
      integer :: nn1, nn2, nn3
      integer :: i,j,k,ii,jj,kk

      type(C_PTR) :: planf, planb

      nn1 = n1/2+1; nn2=n2/2+1; nn3=n3/2+1
      planf = fftw_plan_dft_r2c_3d(n3, n2, n1, f, Fk, FFTW_ESTIMATE)
      planb = fftw_plan_dft_c2r_3d(n3, n2, n1, Fk, f, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(planf, f, Fk)
      do k=1,n3; if (k>nn3) then; kk = n3+1-k; else; kk=k-1; endif
      do j=1,n2; if (j>nn2) then; jj = n2+1-j; else; jj=j-1; endif
         do i=1,nn1
            rad2 = dble((i-1)**2) + dble(jj**2) + dble(kk**2)
            Fk(i,j,k) = -rad2*dk**2*Fk(i,j,k)
         enddo
      enddo
      enddo

      call fftw_execute_dft_c2r(planb, Fk, f)
      call fftw_destroy_plan(planf)
      call fftw_destroy_plan(planb)
      f = f / dble(n1) / dble(n2) / dble(n3)
    end subroutine laplacian_3d

    real(dl) function grad_energy_1d(n, f, Fk, dk)
      integer :: n
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      real(dl) :: dk

      integer :: i,ii,nn
      real(dl) :: GE

      type(C_PTR) :: planf, planb

      nn = n/2+1
      planf = fftw_plan_dft_r2c_1d(n, f, Fk, FFTW_ESTIMATE)
      planb = fftw_plan_dft_c2r_1d(n, Fk, f, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(planf, f, Fk)

      GE = 0.
      do i=1,n
         if (i<=nn) then; ii = i-1; else; ii=n+1-i; endif
         GE = GE + Fk(ii+1)*conjg(Fk(ii+1)) * (ii*dk)**2
      enddo
      GE = GE / n / n  ! add normalizations to average and correct unnormalized inverse DFT

      grad_energy_1d = 0.5*GE
    end function grad_energy_1d

  end module fftw3
