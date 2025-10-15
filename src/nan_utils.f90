module nan_utils
  implicit none
contains

  pure logical function is_nan_r8(x)
    implicit none
    real(8), intent(in) :: x
    ! Fortran 95 portable NaN check: NaN is not equal to itself
    is_nan_r8 = (x /= x)
  end function is_nan_r8

  pure logical function is_finite_r8(x)
    implicit none
    real(8), intent(in) :: x
    ! finite <=> not NaN and not infinite
    is_finite_r8 = .not.(x /= x) .and. (abs(x) < huge(0.0d0))
  end function is_finite_r8

end module nan_utils
