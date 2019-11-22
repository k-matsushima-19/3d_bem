!> 境界メッシュ
module mesh_class
  implicit none
  private

  !> 三角形メッシュ
  type :: mesh
     !> 要素数
     integer :: ne
     !> 節点数
     integer :: np

   contains
     !> 初期化
     procedure :: init
  end type mesh

contains
  !> \details 初期化する
  subroutine init(self)
    class(mesh),intent(inout) :: self
  end subroutine init

end module mesh_class
