! 3D-HelmholtzのBEMの影響係数行列の成分を計算するためのclass
module layerpotential3d_helmholtz_class
  use math
  use misc
  implicit none
  private

  !> 3D-HelmholtzのBEMの影響係数行列の成分を計算するためのclass
  type,public :: layerpotential3d_helmholtz
     !> 角周波数
     real(8) :: w
     !> 位相速度
     real(8) :: c
     !> 波数
     real(8) :: k

     !> [-1,1]上の1次元Gauss積分の積分点の数
     integer :: ng_1d
     !> [-1,1]上の1次元Gauss積分の積分点の座標 (ng)
     real(8),allocatable :: gauss_1d_xi(:)
     !> [-1,1]上の1次元Gauss積分の積分点の重み (ng)
     real(8),allocatable :: gauss_1d_w(:)

     !> 三角形上のGauss積分の積分点の数
     integer :: ng_tri
     !> 三角形上のGauss積分の積分点の面積座標 (3,ng_tri)
     real(8),allocatable :: gauss_tri_xi(:,:)
     !> 三角形上のGauss積分の積分点の重み (ng_tri)
     real(8),allocatable :: gauss_tri_w(:)
   contains
     !> finalizer
     procedure :: destructor
     !> 初期化
     procedure :: init
     !> 1重層,2重層,2重層随伴、2重層微分の三角形上の積分を計算
     procedure :: calc
  end type layerpotential3d_helmholtz

  ! 積分点を計算する関数をCから呼ぶ
  interface
     subroutine generate_gauss(m, x, w) bind(c)
       use iso_c_binding
       integer(c_int),value :: m
       real(c_double) :: x(m)
       real(c_double) :: w(m)
     end subroutine generate_gauss
  end interface

  

contains
  !> finalizer
  !! \param self
  elemental subroutine destructor(self)
    class(layerpotential3d_helmholtz),intent(inout) :: self

    if(allocated(self%gauss_1d_xi)) deallocate(self%gauss_1d_xi)
    if(allocated(self%gauss_1d_w)) deallocate(self%gauss_1d_w)
    if(allocated(self%gauss_tri_xi)) deallocate(self%gauss_tri_xi)
    if(allocated(self%gauss_tri_w)) deallocate(self%gauss_tri_w)
  end subroutine destructor
  
  !> 初期化
  !! \param self
  !! \param ng_1d 1次元Gauss積分の点数
  !! \param w 角周波数
  !! \param c 位相速度
  subroutine init(self, ng_1d, w, c)
    class(layerpotential3d_helmholtz),intent(out) :: self
    integer,intent(in) :: ng_1d
    real(8),intent(in) :: w
    real(8),intent(in) :: c

    self%ng_1d = ng_1d
    self%w = w
    self%c = c

    self%k = self%w / self%c

    ! 積分点の座標と重みを計算
    allocate(self%gauss_1d_xi(ng_1d), self%gauss_1d_w(ng_1d))
    call generate_gauss(self%ng_1d, self%gauss_1d_xi, self%gauss_1d_w)

    ! 三角形Gauss積分の点数はとりあえず6とする
    self%ng_tri = 6
    allocate(self%gauss_tri_xi(3,self%ng_tri), self%gauss_tri_w(self%ng_tri))
    call generate_gauss_triangle(self%ng_tri, self%gauss_tri_xi, self%gauss_tri_w)
    
  end subroutine init
  
  !> 1重層,2重層,2重層随伴、2重層微分の三角形上の積分を計算
  !! \details y1, y2, y3の順番は法線n_yの向きと矛盾していない必要あり
  !! \param self
  !! \param x 選点
  !! \param y1 積分する要素の頂点 (節点) 1の座標
  !! \param y2 積分する要素の頂点 (節点) 2の座標
  !! \param y3 積分する要素の頂点 (節点) 3の座標
  !! \param area 積分する要素の面積
  !! \param n_x 点xでの外向き単位法線ベクトル
  !! \param n_y 積分する要素の外向き単位法線ベクトル
  !! \param S 1重層の積分
  !! \param D 2重層の積分
  !! \param Ds 2重層随伴の積分
  !! \param N 2重層微分の積分
  subroutine calc(self, x, y1, y2, y3, area, n_x, n_y, S, D, Ds, N)
    class(layerpotential3d_helmholtz),intent(in) :: self
    real(8),intent(in) :: x(3)
    real(8),intent(in) :: y1(3), y2(3), y3(3)
    real(8),intent(in) :: area
    real(8),intent(in) :: n_x(3)
    real(8),intent(in) :: n_y(3)
    complex(8),intent(out) :: S, D, Ds, N

    integer :: ig
    real(8) :: y_vec(3), r_vec(3), r
    complex(8) :: G, H, e

    ! Gauss積分
    S = zero
    D = zero
    Ds = zero
    N = zero
    do ig=1,self%ng_tri
       ! 積分点 (Cartesian)
       y_vec = convert_area_to_Cartesian(y1, y2, y3, self%gauss_tri_xi(:,ig))

       r_vec = x - y_vec
       r = vector_length(r_vec)

       ! expはあらかじめ計算しておく
       e = exp(ione*self%k*r)
       
       ! Single layer
       G = e / (4.d0*pi*r)
       S = S + self%gauss_tri_w(ig) * G
       ! Double layer
       H = (1.d0-ione*self%k*r)/(4.d0*pi*r**2) * e * dot_product(r_vec,n_y)/r
       D = D + self%gauss_tri_w(ig) * H
    end do

    ! 面積を掛ける
    S = S * area
    D = D * area
    
  end subroutine calc
  
  !> 基本解
  complex(8) function func_G(self, x, y1, y2, y3, n_x, n_y)
    class(layerpotential3d_helmholtz),intent(in) :: self
    real(8),intent(in) :: x(3)
    real(8),intent(in) :: y1(3), y2(3), y3(3)
    real(8),intent(in) :: n_x(3)
    real(8),intent(in) :: n_y(3)
    
  end function func_G
end module layerpotential3d_helmholtz_class
