!> 3D-Helmholtz方程式を満たす球面波\f$ I_n^m(x-c)=j_n(k|x-c|)Y_n^m((x-c)/|x-c|) \f$
module incwave3d_helmholtz_spherical_class
  use misc
  use math
  use incwave3d_helmholtz_class
  implicit none
  private  

  type,public,extends(incwave3d_helmholtz) :: incwave3d_helmholtz_spherical
     !> 複素振幅
     complex(8) :: amp
     !> 球面波の次数 
     integer :: n
     !> 球面波の次数
     integer :: m
     !> 中心
     real(8) :: origin(3)
   contains
     !> finalizer
     procedure :: destructor
     !> 初期化
     procedure :: init
     !> 点xでの入射波の値を計算
     procedure :: calc
     !> 点xでの入射波の勾配の値を計算
     procedure :: calc_derivative
  end type incwave3d_helmholtz_spherical

 

contains
  !> 初期化
  !! \param self
  !! \param w 角周波数
  !! \param c 位相速度
  !! \param amp 複素振幅
  !! \param n 球面波の次数
  !! \param m 球面波の次数
  !! \param origin 球面波の中心
  subroutine init(self, w, c, amp, n, m, origin)
    class(incwave3d_helmholtz_spherical),intent(out) :: self
    real(8),intent(in) :: w
    real(8),intent(in) :: c
    complex(8),intent(in) :: amp
    integer,intent(in) :: n
    integer,intent(in) :: m
    real(8),intent(in) :: origin(3)
    
    self%w = w
    self%amp = amp
    self%c = c

    self%n = n
    self%m = m
    self%origin = origin

    self%k = self%w / self%c
    
  end subroutine init

  !> finalizer
  !! \param self
  elemental subroutine destructor(self)
    class(incwave3d_helmholtz_spherical),intent(inout) :: self
  end subroutine destructor

  !> 点xでの入射波の値を計算
  !! \param self
  !! \param x 座標
  !! \return 入射波の値
  function calc(self, x)
    class(incwave3d_helmholtz_spherical),intent(in) :: self
    real(8),intent(in) :: x(3)
    complex(8) :: calc

    complex(8) :: Inm((self%n+1)**2)

    call calc_Inm(self%k, x, self%n, Inm)

    calc = Inm(loct(self%n,self%m))

  end function calc

  !> 点xでの入射波の勾配の値を計算
  !! \param self
  !! \param x 座標
  !! \return 入射波の勾配
  function calc_derivative(self, x)
    class(incwave3d_helmholtz_spherical),intent(in) :: self
    real(8),intent(in) :: x(3)
    complex(8) :: calc_derivative(3)

    complex(8) :: Inm_deri(3,(self%n+1)**2)

    call calc_Inm_derivative(self%k, x, self%n, Inm_deri)

    calc_derivative(:) = Inm_deri(:,loct(self%n,self%m))

  end function calc_derivative

  
  
end module incwave3d_helmholtz_spherical_class

