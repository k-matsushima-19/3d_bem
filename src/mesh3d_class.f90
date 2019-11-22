!> 境界メッシュ
module mesh3d_class
  use math
  use misc
  implicit none
  private

  !> 三角形メッシュ
  type,public :: mesh3d
     !> 要素数 (三角形の数)
     integer :: ne
     !> 節点数
     integer :: np

     !> 節点座標の配列 (3,np)
     real(8),allocatable :: p(:,:)
     !> 各要素が有する節点の番号 (3,ne)
     integer,allocatable :: nd(:,:)

     !> 各要素の重心 (3,ne)
     real(8),allocatable :: c(:,:)
     !> 外向き単位法線ベクトル (3,ne)
     real(8),allocatable :: an(:,:)
     !> 各要素の面積 (ne)
     real(8),allocatable :: ar(:)
   contains
     !> 初期化
     procedure :: init
     !> .offファイルを読み込む
     procedure :: read_off
     !> 外向き単位法線ベクトルをファイルに出力する
     procedure :: output_an
     !> 要素の向き (法線の方向) を逆転する
     procedure :: invert
  end type mesh3d

contains
  !> \details ne, np, p, ndを与えた後，その他のメンバを初期化する
  subroutine init(self)
    class(mesh3d),intent(inout) :: self

    integer :: i

    call assert(allocated(self%p))
    call assert(allocated(self%nd))

    ! 重心
    if(allocated(self%c)) deallocate(self%c)
    allocate(self%c(3,self%ne))
    do i=1,self%ne
       self%c(:,i) = (self%p(:,self%nd(1,i))+self%p(:,self%nd(2,i))+self%p(:,self%nd(3,i))) / 3.d0
    end do

    ! 法線と要素面積
    if(allocated(self%an)) deallocate(self%an)
    if(allocated(self%ar)) deallocate(self%ar)
    allocate(self%an(3,self%ne))
    allocate(self%ar(self%ne))
    do i=1,self%ne
       ! (節点1->節点2) x (節点1->節点3)を計算
       self%an(:,i) = vector_product(self%p(:,self%nd(2,i))-self%p(:,self%nd(1,i)), self%p(:,self%nd(3,i))-self%p(:,self%nd(1,i)))
       ! 面積 = ↑の外積の長さ / 2
       self%ar(i) = vector_length(self%an(:,i)) * 0.5d0
       ! 単位ベクトルにする
       self%an(:,i) = self%an(:,i) / (self%ar(i)*2.d0)       
    end do
    
  end subroutine init

  !> \details .off (Object File Format) ファイルを読み込む
  !! \details フォーマットは https://shape.cs.princeton.edu/benchmark/documentation/off_format.html を参照
  !! \param self
  !! \param filename ファイル名
  subroutine read_off(self, filename)
    class(mesh3d),intent(inout) :: self
    character(*),intent(in) :: filename

    character(1000) :: ctmp
    integer :: itmp, i

    open(10,file=filename,status="old")

    ! 一行目は無視
    read(10,*) ctmp
    ! 節点数，要素数，(辺の数)をread
    read(10,*) self%np, self%ne, itmp
    ! 各節点の座標をread
    allocate(self%p(3,self%np))
    do i=1,self%np
       read(10,*) self%p(1:3,i)
    end do
    ! 各要素の節点番号をread
    allocate(self%nd(3,self%ne))
    do i=1,self%ne
       ! 1列目は節点の数
       read(10,*) itmp, self%nd(1:3,i)
       ! 三角形メッシュのみ対応
       call assert(itmp == 3)
    end do
    ! .offの節点番号は0から始まるので，1足す
    self%nd(:,:) = self%nd(:,:) + 1
    close(10)

    call self%init
    
  end subroutine read_off

  !> 外向き単位法線ベクトルをファイルに出力する
  !! \param self
  !! \param filename ファイル名
  subroutine output_an(self, filename)
    class(mesh3d),intent(in) :: self
    character(*),intent(in) :: filename

    integer :: i

    open(10,file=filename,status="replace")

    do i=1,self%ne
       write(10,'(6e24.16)') self%c(1:3,i), self%an(1:3,i)
    end do

    close(10)
  end subroutine output_an

  !> 要素の向き (法線の方向) を逆転する
  subroutine invert(self)
    class(mesh3d),intent(inout) :: self

    integer :: i, itmp

    ! 各要素の節点2と節点3を入れ替え
    do i=1,self%ne
       itmp = self%nd(3,i)
       self%nd(3,i) = self%nd(2,i)
       self%nd(2,i) = itmp
    end do

    ! 法線を反転
    self%an(:,:) = -self%an(:,:)
  end subroutine invert

end module mesh3d_class
