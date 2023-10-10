module Generic_Derivs
  use FFTW
  use Fast_Cheby

  implicit none

  interface laplacian_1d
     subroutine laplacian_1d_wtype(tPair,dk)
       type(transformPair1D), intent(inout) :: tPair
       real(dl), intent(in) :: dk
     end subroutine laplacian_1d_wtype

     subroutine laplacian_cheby_1d_mapped(tPair)
       type(chebyshevPair1D), intent(inout) :: tPair
     end subroutine laplacian_cheby_1d_mapped
  end interface laplacian_1d
  
end module Generic_Derivs
