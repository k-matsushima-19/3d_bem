module math
  
contains
  !> ベクトルの外積
  !! \param vec1 3次元ベクトル\f$ a \f$
  !! \param vec2 3次元ベクトル\f$ b \f$
  !! \return \f$a \times b\f$
  function vector_product(vec1, vec2) result(out)
    real(8),intent(in) :: vec1(3)
    real(8),intent(in) :: vec2(3)
    real(8) :: out(3)

    out(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    out(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    out(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  end function vector_product

  !> ベクトルの長さ
  !! \param vec 3次元ベクトル\f$ a \f$
  !! \return \f$ |a| \f$
  real(8) function vector_length(vec)
    real(8),intent(in) :: vec(3)

    vector_length = sqrt(dot_product(vec,vec))
  end function vector_length
end module math
