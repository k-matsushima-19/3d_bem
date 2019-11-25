module math
  use misc
  implicit none
  private :: rdot_product1, rdot_product2, rdot_product3, rdot_product4
  
  complex(8),parameter :: zero = (0.d0,0.d0)
  complex(8),parameter :: ione = (0.d0,1.d0)
  complex(8),parameter ::  one = (1.d0,0.d0)
  real(8),parameter :: pi=3.1415926535897932d0

  !> 内積でない (共役を取らない) ドット積
  interface rdot_product
     module procedure rdot_product1
     module procedure rdot_product2
     module procedure rdot_product3
     module procedure rdot_product4
  end interface rdot_product

  interface
     ! gsl/gsl_sf_bessel.h
     ! int gsl_sf_bessel_jl_array(int lmax, double x, double result_array[])
     function gsl_sf_bessel_jl_array(lmax, x, result_array) result(out) bind(c,name="gsl_sf_bessel_jl_array")
       use iso_c_binding
       integer(c_int),value :: lmax
       real(c_double),value :: x
       real(c_double) :: result_array(0:lmax)
       integer(c_int) :: out
     end function gsl_sf_bessel_jl_array

     ! gsl/gsl_sf_bessel.h
     ! int gsl_sf_bessel_yl_array(int lmax, double x, double result_array[])
     function gsl_sf_bessel_yl_array(lmax, x, result_array) result(out) bind(c,name="gsl_sf_bessel_yl_array")
       use iso_c_binding
       integer(c_int),value :: lmax
       real(c_double),value :: x
       real(c_double) :: result_array(0:lmax)
       integer(c_int) :: out
     end function gsl_sf_bessel_yl_array
  end interface
  
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

  !> 三角形上のGauss積分の積分点の面積座標と対応する重みを返す
  !! \details See https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620070316
  !! \param ng 積分点の数
  !! \param gauss_xi 積分点の面積座標
  !! \param gauss_w 積分点の重み
  subroutine generate_gauss_triangle(ng, gauss_xi, gauss_w)
    integer,intent(in) :: ng
    real(8),intent(out) :: gauss_xi(3,ng)
    real(8),intent(out) :: gauss_w(ng)
    
    call assert(ng == 6, "ng=6以外は未対応")

    if(ng == 6) then
       ! 面積座標
       gauss_xi(1:3,1) = [0.659027622374092d0, 0.231933368553031d0, 0.109039009072877d0]
       gauss_xi(1:3,2) = [0.659027622374092d0, 0.109039009072877d0, 0.231933368553031d0]
       gauss_xi(1:3,3) = [0.231933368553031d0, 0.659027622374092d0, 0.109039009072877d0]
       gauss_xi(1:3,4) = [0.231933368553031d0, 0.109039009072877d0, 0.659027622374092d0]
       gauss_xi(1:3,5) = [0.109039009072877d0, 0.659027622374092d0, 0.231933368553031d0]
       gauss_xi(1:3,6) = [0.109039009072877d0, 0.231933368553031d0, 0.659027622374092d0]

       ! 重み
       gauss_w(1) = 1.d0 / 6.d0
       gauss_w(2) = 1.d0 / 6.d0
       gauss_w(3) = 1.d0 / 6.d0
       gauss_w(4) = 1.d0 / 6.d0
       gauss_w(5) = 1.d0 / 6.d0
       gauss_w(6) = 1.d0 / 6.d0
    end if
  end subroutine generate_gauss_triangle

  !> 三角形の面積座標からCartesian座標へ変換
  !! \param x1 三角形の頂点1の3次元Cartesian座標
  !! \param x2 三角形の頂点2の3次元Cartesian座標
  !! \param x3 三角形の頂点3の3次元Cartesian座標
  !! \param a 面積座標
  !! \return aに対応する3次元Cartesian座標
  function convert_area_to_Cartesian(x1, x2, x3, a) result(out)
    real(8),intent(in) :: x1(3), x2(3), x3(3)
    real(8),intent(in) :: a(3)
    real(8) :: out(3)

#ifdef DEBUG
    ! 面積座標の3つの成分は足して1のはず
    call assert(abs(sum(a) - 1.d0) < 1d-14, "Invalid area coordinate")
#endif

    out(1) = x1(1)*a(1) + x2(1)*a(2) + x3(1)*a(3)
    out(2) = x1(2)*a(1) + x2(2)*a(2) + x3(2)*a(3)
    out(3) = x1(3)*a(1) + x2(3)*a(2) + x3(3)*a(3)
    
  end function convert_area_to_Cartesian

  ! 共役を取らない内積 (的なもの)
  complex(8) function rdot_product1(a, b)
    complex(8),intent(in) :: a(:)
    real(8),intent(in) :: b(:)

    rdot_product1 = dot_product(b, a)
  end function rdot_product1

  complex(8) function rdot_product2(a, b)
    real(8),intent(in) :: a(:)
    complex(8),intent(in) :: b(:)

    rdot_product2 = dot_product(a, b)
  end function rdot_product2

  complex(8) function rdot_product3(a, b)
    complex(8),intent(in) :: a(:)
    complex(8),intent(in) :: b(:)

    rdot_product3 = dot_product(conjg(a), b)
  end function rdot_product3

  complex(8) function rdot_product4(a, b)
    real(8),intent(in) :: a(:)
    real(8),intent(in) :: b(:)

    rdot_product4 = dot_product(a, b)
  end function rdot_product4

  !> 球面座標からデカルト座標を計算
  !! \param r 半径
  !! \param theta 方位角 (azimuthal angle) [rad]
  !! \param phi 極角 (polar angle) [rad]
  !! \return 3次元デカルト座標
  function spherical_to_Cartesian(r, theta, phi) result(out)
    real(8),intent(in) :: r
    real(8),intent(in) :: theta
    real(8),intent(in) :: phi
    real(8) :: out(3)

    real(8) :: sin_theta

    sin_theta = sin(theta)

    out(1) = r*sin_theta*cos(phi)
    out(2) = r*sin_theta*sin(phi)
    out(3) = r*cos(theta)
    
  end function spherical_to_Cartesian

  !> デカルト座標から球面座標を計算
  !! \details \f$ x_1=x_2=0 \f$のとき\f$ \phi \f$は不定であり，さらに\f$x_3=0 \f$であるとき\f$ \theta \f$も不定であるが，ここでは\f$ \phi=\pi/2,\ \theta=0 \f$を返す．
  !! \param x 3次元デカルト座標
  !! \param r 半径 \f$ r>0 \f$
  !! \param theta 方位角 (azimuthal angle) [rad] \f$ \theta\in[0,\pi] \f$
  !! \param phi 極角 (polar angle) [rad] \f$ \phi\in(-\pi,\pi] \f$
  subroutine Cartesian_to_spherical(x, r, theta, phi)
    real(8),intent(in) :: x(3)
    real(8),intent(out) :: r
    real(8),intent(out) :: theta
    real(8),intent(out) :: phi

    r = vector_length(x)
    if(r > 1d-32*abs(x(3))) then
       theta = acos(x(3)/r)
    else
       theta = 0.d0
    end if
    phi = atan2(x(2),x(1))
    
  end subroutine Cartesian_to_spherical

  !> 球面調和関数\f$ Y_n^m(x) \f$のインデックス\f$ m,n \f$を1から始まる1つのindexに置き換える
  !! \param n 球面調和関数のインデックス
  !! \param m 球面調和関数のインデックス
  integer function loct_sph(n,m)
    integer,intent(in) :: n, m
    loct_sph=n*n+n+1+m
  end function loct_sph
  
  !> 球面調和関数\f$ x\in R^3\mapsto Y_n^m(x)=Y_n^m(x/|x|) \f$を\f$ n=0 \f$から与えられた\f$ n \f$まで計算する。ynm(loct_sph(n,m))に格納される
  !! \param n インデックスの最大値 \f$ \geq 0 \f$
  !! \param x 球面調和関数の引数 (単位ベクトルでなくても良い)
  !! \param ynm 結果の配列 (長さ(n+1)**2)
  subroutine sph(n,x,ynm)
    integer,intent(in) :: n
    real(8),intent(in) :: x(3)
    complex(8),intent(out) :: ynm((n+1)**2)
    
    integer :: m, i, k, j
    complex(8) :: cc
    real(8) :: v(3), ctmp, xlen

    xlen=sqrt(dot_product(x,x))

    do i=1,3
       v(i)=x(i)/xlen
    enddo

    cc=cmplx(v(1),v(2),kind(1.d0))

    ynm(loct_sph(0,0))=1d0

    do m=0,n-1
       ynm(loct_sph(m+1,m+1))=cc*ynm(loct_sph(m,m))*sqrt((2d0*m+1d0)/(2d0*m+2d0))
    enddo


    do m=0,n-1
       ynm(loct_sph(m+1,m))=ynm(loct_sph(m,m))*v(3)*sqrt(2d0*m+1d0)
       do i=m+2,n
          ctmp=sqrt(1d0*(i+m)*(i-m))
          ynm(loct_sph(i,m))=(v(3)*(2d0*i-1d0)*ynm(loct_sph(i-1,m))-ynm(loct_sph(i-2,m))*sqrt(1d0*(i-1+m)*(i-1-m)))/ctmp
       enddo
    enddo

    do k=-1,-n,-1
       do j=-k,n
          ynm(loct_sph(j,k))=(-1.0)**(-k)*conjg(ynm(loct_sph(j,-k)))
       enddo
    enddo

    return
  end subroutine sph

  !> 0からnmax次までの球Bessel関数\f$ j_n(x) \f$を計算
  !! \param nmax 次数の最大値 \f$ \geq 0 \f$
  !! \param x 引数
  !! \param out 結果の配列
  subroutine bessel_sph_jn_array(nmax, x, out)
    integer,intent(in) :: nmax
    real(8),intent(in) :: x
    real(8),intent(out) :: out(0:nmax)

    integer :: ret

    ret = gsl_sf_bessel_jl_array(nmax, x, out)

    call assert(ret == 0)
  end subroutine bessel_sph_jn_array

  !> 0からnmax次までの球Bessel関数\f$ y_n(x) \f$を計算
  !! \param nmax 次数の最大値 \f$ \geq 0 \f$
  !! \param x 引数
  !! \param out 結果の配列
  subroutine bessel_sph_yn_array(nmax, x, out)
    integer,intent(in) :: nmax
    real(8),intent(in) :: x
    real(8),intent(out) :: out(0:nmax)

    integer :: ret

    ret = gsl_sf_bessel_yl_array(nmax, x, out)

    call assert(ret == 0)
  end subroutine bessel_sph_yn_array

  !> 0からnmax次までの第1種球Hankel関数\f$ h_n^{(1)}(x) \f$を計算
  !! \param nmax 次数の最大値 \f$ \geq 0 \f$
  !! \param x 引数
  !! \param out 結果の配列
  subroutine bessel_sph_hn_array(nmax, x, out)
    integer,intent(in) :: nmax
    real(8),intent(in) :: x
    complex(8),intent(out) :: out(0:nmax)

    integer :: ret
    real(8) :: j(0:nmax), y(0:nmax)

    ret = gsl_sf_bessel_jl_array(nmax, x, j)
    call assert(ret == 0)

    ret = gsl_sf_bessel_yl_array(nmax, x, y)
    call assert(ret == 0)

    out = dcmplx(j,y)
  end subroutine bessel_sph_hn_array

end module math
