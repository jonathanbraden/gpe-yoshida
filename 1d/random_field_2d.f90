module gaussianRandomField_2D
  use, intrinsic :: iso_c_binding
  use constants
  use fftw3

  implicit none

contains

  subroutine generate_2dGRF(field, spectrum, convolve)
    real(C_DOUBLE), dimension(:,:), intent(inout) :: field

  end subroutine generate_2dGRF

  subroutine generate_2dGRF_complex(field, spectrum, convolve)
    real(C_DOUBLE_COMPLEX), dimension(:,:), intent(inout) :: field
    real(dl), dimension(:), intent(in) :: spectrum
    logical, intent(in) :: convolve

    complex(C_DOUBLE_COMPLEX), dimension(1:size(field,dim=1) :: deviate
    complex(C_DOUBLE_COMPLEX), allocatable :: Fk(:,:)
    real(C_DOUBLE), allocatable, dimension(1:size(field,dim=1)) :: amp, phase
    integer :: nx, nnx, ny, nny
    integer :: i,j, jj
    type(C_PTR) :: fft_plan

    nx = size(field,dim=1); ny = size(field,dim=2)
    nnx = nx/2+1; nny = ny/2+1

    !initialize random numbers

    ! Create plan
    !fft_plan = fftw_plan_dft_c2c_2d(  Fk, field, FFTW_ESTIMATE )
    
    do j=1,ny; if (j>nny) then; jj = ; else; jj=j-1
       do i=1,nnx
          
       enddo
       do i=nnx+1,nx
       enddo
    enddo
       
  end subroutine generate_2dGRF_complex
