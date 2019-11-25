!> BEMの抽象クラス
module bem3d_class
  use math
  use misc
  use mesh3d_class
  implicit none

  !> bemの抽象クラス
  type,public,abstract :: bem3d
     !> mesh
     type(mesh3d),pointer :: mesh => null()

     !> 係数行列のサイズ
     integer :: size

     !> 係数行列 (size,size)
     complex(8),allocatable :: amat(:,:)

     !> 右辺ベクトル (size)
     complex(8),allocatable :: bvec(:)

     !> ipiv (size)
     integer,allocatable :: ipivot(:)

     !> 境界条件
     character(20) :: b_condition

     !> 係数行列amatがLU分解されているか
     logical :: is_factorized_amat = .false.

   contains
     !> 係数行列amatをLU分解する
     procedure :: lu_decompose
     !> LU分解した係数行列amatとbvecを用いて連立方程式を解く
     procedure :: lu_solve
  end type bem3d

contains
  !> 係数行列amatをLU分解する
  !! \param self
  subroutine lu_decompose(self)
    class(bem3d),intent(inout) :: self

    integer :: info


    allocate(self%ipivot(self%size))    
    call zgetrf(self%size, self%size, self%amat, self%size, self%ipivot, info)
    call assert(info == 0)

    self%is_factorized_amat = .true.
    
  end subroutine lu_decompose

  !> LU分解した係数行列amatとbvecを用いて連立方程式を解く
  !! \param self
  subroutine lu_solve(self)
    class(bem3d),intent(inout) :: self
    integer :: info

    call assert(self%is_factorized_amat)


    call zgetrs('N', self%size, 1, self%amat, self%size, self%ipivot, self%bvec, self%size, info)    
    call assert(info == 0)
  end subroutine lu_solve
  
end module bem3d_class
