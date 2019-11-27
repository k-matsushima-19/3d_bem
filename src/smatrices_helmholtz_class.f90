!> 3D-HelmholtzのS行列を複数用いて，multiple scatteringを解くための構造体
module smatrices_helmholtz_class
  use misc
  use math
  use smatrix_helmholtz_class
  use solution_helmholtz_class
  use mesh3d_class
  implicit none
  private

  !> 3D-HelmholtzのS行列を複数用いて，multiple scatteringを解くための構造体
  type,public :: smatrices_helmholtz
     !> 角周波数
     real(8) :: w
     !> 位相速度 1:外部 2:内部
     real(8) :: cs(2)
     !> 外部領域 (領域1) の波数
     real(8) :: k

     !> 次数nの最大値
     integer :: nmax
     !> S行列のサイズ
     integer :: size

     !> 散乱体の数
     integer :: nobj
     !> それぞれの散乱体のS行列 (nobj)
     type(smatrix_helmholtz),allocatable :: smats(:)

     !> TとSを並べた行列
     complex(8),allocatable :: amat(:,:)
     !> amatのサイズ = size*nobj
     integer :: size_all
     !> amatのLU分解に使うipiv
     integer,allocatable :: ipivot(:)
     !> 右辺ベクトル
     complex(8),allocatable :: bvec(:)
     !> 各散乱体へ外から入射する波の係数
     complex(8),allocatable :: alphas(:,:)

     !> 係数行列amatがLU分解されているか
     logical :: is_factorized_amat = .false.
   contains
     !> 初期化
     procedure :: init
     !> 単位行列，S行列，translation matrixを並べた係数行列amatを作る
     procedure :: gen_amat
     !> amatをLU分解
     procedure :: lu_decompose
     !> LU分解した係数行列amatとbvecを用いて連立方程式を解く
     procedure :: lu_solve
     !> 平面波\f$ \mathrm{amp}*\exp(\mathrm{i}kp\cdot x) \f$が入射波するときの右辺ベクトルを生成
     procedure :: gen_bvec_planewave
     !> 現在のbvecから境界・内点の解を計算
     procedure :: gen_solution
  end type smatrices_helmholtz

contains
  !> 初期化
  !! \param self
  !! \param nobj 散乱体の数
  !! \param smats それぞれの散乱体のS行列
  subroutine init(self, nobj, smats)
    class(smatrices_helmholtz),intent(out) :: self
    integer,intent(in) :: nobj
    type(smatrix_helmholtz),intent(in) :: smats(nobj)

    integer :: iobj

    self%nobj = nobj
    allocate(self%smats(self%nobj))
    do iobj=1,nobj
       self%smats(iobj) = smats(iobj)
    end do

    self%w = self%smats(1)%w
    self%cs = self%smats(1)%cs
    self%k = self%smats(1)%k

    ! ややこしいので，すべてのS行列のnmaxは同じであるとする
    do iobj=2,self%nobj
       call assert(smats(iobj)%nmax == smats(1)%nmax)
    end do
    self%nmax = self%smats(1)%nmax
    self%size = (self%nmax+1)**2
    self%size_all = self%size * self%nobj

    allocate(self%amat(self%size_all,self%size_all))
    allocate(self%alphas(self%size,self%nobj))
    allocate(self%bvec(self%size_all))
    
  end subroutine init

  !> 単位行列，S行列，translation matrixを並べた係数行列amatを作る
  !! \param self
  subroutine gen_amat(self)
    class(smatrices_helmholtz),intent(inout) :: self

    integer :: i, k, l
    complex(8),allocatable :: T(:,:)
    integer :: ind(2)

    
    ! 対角成分を1に
    self%amat(:,:) = zero
    do i=1,self%size_all
       self%amat(i,i) = one
    end do

    ! Translation matrixをamatの対応する位置に配置
    allocate(T(self%size,self%size))
    do l=1,self%nobj
       do k=1,self%nobj
          ! 対角ブロックはなにもしない
          if(k == l) cycle
          
          ! ブロックの左上の位置
          ind(1) = self%size*(k-1)+1
          ind(2) = self%size*(l-1)+1

          ! T^(kl)の計算
          call self%smats(1)%calc_T(self%smats(k)%mesh%center-self%smats(l)%mesh%center, T)

          ! -S*Tを配置
          self%amat(ind(1):ind(1)+self%size-1,ind(2):ind(2)+self%size-1) = -matmul(self%smats(k)%S,T)
       end do
    end do
    
  end subroutine gen_amat

  !> 係数行列amatをLU分解する
  !! \param self
  subroutine lu_decompose(self)
    class(smatrices_helmholtz),intent(inout) :: self

    integer :: info


    allocate(self%ipivot(self%size_all))    
    call zgetrf(self%size_all, self%size_all, self%amat, self%size_all, self%ipivot, info)
    call assert(info == 0)

    self%is_factorized_amat = .true.
    
  end subroutine lu_decompose
  
  !> LU分解した係数行列amatとbvecを用いて連立方程式を解く．bvecに解が上書きされる．
  !! \param self
  subroutine lu_solve(self)
    class(smatrices_helmholtz),intent(inout) :: self
    integer :: info

    call assert(self%is_factorized_amat)


    call zgetrs('N', self%size_all, 1, self%amat, self%size_all, self%ipivot, self%bvec, self%size_all, info)    
    call assert(info == 0)
  end subroutine lu_solve

  !> 平面波\f$ \mathrm{amp}*\exp(\mathrm{i}kp\cdot x) \f$が入射波するときの右辺ベクトルを生成
  !! \param self
  !! \param amp 複素振幅
  !! \param pvec 伝播方向を表すベクトル\f$ p\in S^2 \f$
  subroutine gen_bvec_planewave(self, amp, pvec)
    class(smatrices_helmholtz),intent(inout) :: self
    complex(8),intent(in) :: amp
    real(8),intent(in) :: pvec(3)
    
    integer :: k
    complex(8),allocatable :: alpha(:)
    complex(8),allocatable :: Ys(:)
    integer :: n, m

    allocate(Ys(self%size))
    allocate(alpha(self%size))
    ! 右辺ベクトルの第kブロック
    do k=1,self%nobj
       ! 第k散乱体の原点まわりの入射波ベクトルを計算
       call sph(self%nmax, pvec, Ys)
       do n=0,self%nmax
          do m=-n,n
             alpha(loct(n,m)) = ione**n*(2*n+1)*conjg(Ys(loct(n,m)))
          end do
       end do
       alpha = alpha * exp(ione*self%k*rdot_product(pvec,self%smats(k)%mesh%center))
       ! selfで保持しておく
       self%alphas(:,k) = alpha(:)

       ! S(k)*alphaをbvecに配置
       self%bvec(1+self%size*(k-1):self%size*k) = matmul(self%smats(k)%S,alpha)
    end do
    
  end subroutine gen_bvec_planewave

  !> 現在のbvecから境界・内点の解を計算
  !! \param self
  !! \param sol 解
  subroutine gen_solution(self, sol)
    class(smatrices_helmholtz),intent(in) :: self
    type(solution_helmholtz),intent(out) :: sol

    integer :: k, l
    complex(8),allocatable :: A(:), T(:,:), B(:)
    type(solution_helmholtz),allocatable :: sols(:)
    integer :: n, m
    ! allocatable,targetはサブルーチンを抜けた瞬間にdeallocされることに注意 (sol%mesh=>meshしても，抜けた瞬間にsol%meshの宛先は不定になる)
    ! type(mesh3d),allocatable,target :: mesh
    type(mesh3d),pointer :: mesh
    integer :: itmp

    allocate(A(self%size))
    allocate(T(self%size,self%size))
    allocate(B(self%size))
    allocate(sols(self%nobj))
    
    do k=1,self%nobj
       ! 第k散乱体のincoming coefficients Aを計算
       ! A(k) = alpha(k) + \sum_{l\neq k} T(k,l)*B(l)
       A(:) = self%alphas(:,k)
       do l=1,self%nobj
          if(k == l) cycle

          ! B(l)
          B(:) = self%bvec(1+self%size*(l-1):self%size*l)
          ! T^(kl)
          call self%smats(1)%calc_T(self%smats(k)%mesh%center-self%smats(l)%mesh%center, T)
          A(:) = A(:) + matmul(T,B)
       end do

       ! Aとself%smats(k)%solsを用いて解を作る
       call sols(k)%init(self%smats(k)%mesh, self%w, self%cs(1))
       sols(k)%u_bndry(:,:) = zero
       ! 各基底の解を重ね合わせる
       do n=0,self%nmax
          do m=-n,n
             sols(k)%u_bndry(:,:) = sols(k)%u_bndry(:,:) + A(loct(n,m))*self%smats(k)%sols(loct(n,m))%u_bndry(:,:)
          end do
       end do
       
    end do

    ! call sols(1)%output_bndry("bndry1.dat")
    ! call sols(2)%output_bndry("bndry2.dat")

    !
    ! solsを結合して1つの構造体を作る
    !
    ! meshの作成
    allocate(mesh)
    mesh = self%smats(1)%mesh
    ! 2~nobj番散乱体を結合
    do k=2,self%nobj
       call mesh%add_mesh(self%smats(k)%mesh)
    end do
    ! ↑で作ったsolsを結合
    call sol%init(mesh, self%w, self%cs(1))
    itmp = 1
    do k=1,self%nobj
       sol%u_bndry(itmp:itmp-1+self%smats(k)%mesh%ne,1:2) = sols(k)%u_bndry(:,1:2)
       itmp = itmp + self%smats(k)%mesh%ne
    end do

    ! call sol%output_bndry("bndry.dat")
    ! stop
    
  end subroutine gen_solution

end module smatrices_helmholtz_class
