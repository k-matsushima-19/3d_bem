!> 3D-Helmholtzを解くためのBEM
module bem3d_helmholtz_class
  use math
  use misc
  use mesh3d_class
  use bem3d_class
  use layerpotential3d_helmholtz_class
  use incwave3d_helmholtz_planewave_class
  use solution_helmholtz_class

  !> 3D-Helmholtzを解くためのBEM
  type,public,extends(bem3d) :: bem3d_helmholtz
     !> 角周波数
     real(8) :: w
     !> 位相速度 1:外部 2:内部
     real(8) :: cs(2)
   contains
     !> 初期化
     procedure :: init
     !> 係数行列を作る
     procedure :: gen_amat
     !> 右辺ベクトルを作る
     procedure :: gen_bvec
     !> 解を返す
     procedure :: gen_solution
  end type bem3d_helmholtz

contains
  !> meshと各パラメータを渡してオブジェクトを初期化
  !! \param self
  !! \param mesh BEMに使う境界mesh
  !! \param b_condition 境界条件 "neumann"
  !! \param w 角周波数
  !! \param cs 各領域の位相速度 1:外部 2:内部
  subroutine init(self, mesh, b_condition, w, cs)
    class(bem3d_helmholtz),intent(out) :: self
    type(mesh3d),intent(inout),target :: mesh
    character(*),intent(in) :: b_condition
    real(8),intent(in) :: w
    real(8),intent(in) :: cs(2)

    self%mesh => mesh
    self%b_condition = b_condition

    self%w = w
    self%cs = cs

    self%size = self%mesh%ne

    call assert(self%b_condition == "neumann")
    
  end subroutine init

  !> 係数行列を作る
  !! \param self
  subroutine gen_amat(self)
    class(bem3d_helmholtz),intent(inout) :: self

    integer :: i, j
    complex(8) :: S, D, Ds, N
    type(layerpotential3d_helmholtz),allocatable :: lps(:)

    allocate(lps(2))
    call lps(1)%init(10, self%w, self%cs(1))
    call lps(2)%init(10, self%w, self%cs(2))

    allocate(self%amat(self%size,self%size))
    self%amat(:,:) = zero

    if(self%b_condition == "neumann") then       
       do j=1,self%mesh%ne
          do i=1,self%mesh%ne
             ! layer potentialを計算
             call lps(1)%calc(&
                  self%mesh%c(:,i), &
                  self%mesh%p(:,self%mesh%nd(1,j)),&
                  self%mesh%p(:,self%mesh%nd(2,j)),&
                  self%mesh%p(:,self%mesh%nd(3,j)),&
                  self%mesh%ar(j),&
                  self%mesh%an(:,i),&
                  self%mesh%an(:,j),&
                  S, D, Ds, N)

             self%amat(i,j) = D
          end do
       end do

       ! 対角に1/2を足す
       do i=1,self%mesh%ne
          self%amat(i,i) = self%amat(i,i) + 0.5d0
       end do
    end if
  end subroutine gen_amat

  !> 右辺ベクトルを作る
  !! \param self
  !! \param incw 入射波
  subroutine gen_bvec(self, incw)
    class(bem3d_helmholtz),intent(inout) :: self
    type(incwave3d_helmholtz_planewave),intent(in) :: incw

    integer :: i
    complex(8) :: u

    allocate(self%bvec(self%size))

    if(self%b_condition == "neumann") then

       do i=1,self%mesh%ne
          u = incw%calc(self%mesh%c(:,i))
          self%bvec(i) = u
       end do
       
    end if
  end subroutine gen_bvec

  !> 解を返す
  !! \param self
  !! \param sol 解
  subroutine gen_solution(self, sol)
    class(bem3d_helmholtz),intent(in) :: self
    type(solution_helmholtz),intent(out) :: sol

    call sol%init(self%mesh)

    if(self%b_condition == "neumann") then
       sol%u_bndry(:,1) = self%bvec(:)
       sol%u_bndry(:,2) = zero
    end if
  end subroutine gen_solution
  
end module bem3d_helmholtz_class
