!> 3D-Helmholtzの入射波の抽象クラス
module incwave3d_helmholtz_class
  implicit none
  private

  type,public,abstract :: incwave3d_helmholtz
     !> 角周波数
     real(8) :: w
     !> 位相速度
     real(8) :: c
     !> 波数
     real(8) :: k

   contains
     !> 点xでの入射波の値を計算
     procedure(calc_def),deferred :: calc
     !> 点xでの入射波の勾配の値を計算
     procedure(calc_derivative_def),deferred :: calc_derivative
  end type incwave3d_helmholtz

  interface
     function calc_def(self, x)
       import incwave3d_helmholtz 
       class(incwave3d_helmholtz),intent(in) :: self
       real(8),intent(in) :: x(3)
       complex(8) :: calc_def
     end function calc_def

     function calc_derivative_def(self, x)
       import incwave3d_helmholtz
       class(incwave3d_helmholtz),intent(in) :: self
       real(8),intent(in) :: x(3)
       complex(8) :: calc_derivative_def(3)
     end function calc_derivative_def
  end interface

end module incwave3d_helmholtz_class
