!> 3D-Helmholtz方程式を満たす平面波
module incwave3d_helmholtz_planewave_class
  use misc
  use math
  use incwave3d_helmholtz_class
  implicit none
  private  

  type,public,extends(incwave3d_helmholtz) :: incwave3d_helmholtz_planewave
     !> 伝播方向の単位ベクトル
     complex(8) :: pvec(3)
     !> 複素振幅
     complex(8) :: amp
   contains
     !> finalizer
     procedure :: destructor
     !> 初期化
     procedure :: init
     !> 点xでの入射波の値を計算
     procedure :: calc
     !> 点xでの入射波の勾配の値を計算
     procedure :: calc_derivative
  end type incwave3d_helmholtz_planewave

 

contains
  !> 初期化
  !! \param self
  !! \param w 角周波数
  !! \param c 位相速度
  !! \param amp 複素振幅
  !! \param pvec 伝播方向を表すベクトル\f$ p\in C^3 \f$ \f$ (p_1^2+p_2^2+p_3^2=1) \f$
  subroutine init(self, w, c, amp, pvec)
    class(incwave3d_helmholtz_planewave),intent(out) :: self
    real(8),intent(in) :: w
    real(8),intent(in) :: c
    complex(8),intent(in) :: amp
    complex(8),intent(in) :: pvec(3)    
    
    self%w = w
    self%amp = amp
    self%pvec(:) = pvec(:)
    self%c = c

    self%k = self%w / self%c
    
  end subroutine init

  !> finalizer
  !! \param self
  elemental subroutine destructor(self)
    class(incwave3d_helmholtz_planewave),intent(inout) :: self
  end subroutine destructor

  !> 点xでの入射波の値を計算
  !! \param self
  !! \param x 座標
  !! \return 入射波の値
  function calc(self, x)
    class(incwave3d_helmholtz_planewave),intent(in) :: self
    real(8),intent(in) :: x(3)
    complex(8) :: calc

    calc = self%amp * exp(ione*self%k*rdot_product(self%pvec,x))

  end function calc

  !> 点xでの入射波の勾配の値を計算
  !! \param self
  !! \param x 座標
  !! \return 入射波の勾配
  function calc_derivative(self, x)
    class(incwave3d_helmholtz_planewave),intent(in) :: self
    real(8),intent(in) :: x(3)
    complex(8) :: calc_derivative(3)

    calc_derivative(:) = self%amp * ione*self%k*self%pvec(:)*exp(ione*self%k*rdot_product(self%pvec,x))

  end function calc_derivative

  
  
end module incwave3d_helmholtz_planewave_class

