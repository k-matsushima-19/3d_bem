!> 3D-HelmholtzのS行列
module smatrix_helmholtz_class
  use misc
  use math
  use mesh3d_class
  use solution_helmholtz_class
  use bem3d_helmholtz_class
  use incwave3d_helmholtz_spherical_class
  implicit none
  private

  !> 3D-HelmholtzのS行列
  type,public :: smatrix_helmholtz
     !> 境界メッシュ
     type(mesh3d),pointer :: mesh => null()
     !> 角周波数
     real(8) :: w
     !> 位相速度 1:外部 2:内部
     real(8) :: cs(2)
     !> 外部領域 (領域1) の波数
     real(8) :: k

     !> 次数nの最大値
     integer :: nmax
     
     !> 境界条件
     character(20) :: b_condition

     !> S行列のサイズ
     integer :: size

     !> S行列 \f$ S_{kn}^{lm} \f$ (1:size,1:size)
     !! \details \f$ S_{kn}^{lm} \f$はS(loct(k,l),loct(n,m))に格納される
     complex(8),allocatable :: S(:,:)

     !> \f$ S_{kn}^{lm} \f$の(n,m)に対応する解 (1:size)
     !! \details (n,m)に対応する解はsols(loct(n,m))に格納される
     type(solution_helmholtz),allocatable :: sols(:)
   contains
     !> BEMを用いてS行列を作る
     procedure :: init
     !> 散乱体が球の場合について，解析解を用いてS行列を作る
     procedure :: init_sphere
     !> incoming coefficients Aを与えたときの解solを計算する
     procedure :: calc_sol
     !> \f$ p\in S^2 \f$方向に進む平面波を表すincoming coefficients Aを計算
     procedure :: calc_A_planewave
     !> Sの成分をファイルに出力
     procedure :: output_S
  end type smatrix_helmholtz

contains
  !> BEMを用いてS行列を作る
  !! \param self
  !! \param mesh 境界mesh
  !! \param b_condition 境界条件 "neumann"
  !! \param w 角周波数
  !! \param cs 各領域の位相速度 1:外部 2:内部
  !! \param nmax 次数nの最大値
  subroutine init(self, mesh, b_condition, w, cs, nmax)
    class(smatrix_helmholtz),intent(out) :: self
    type(mesh3d),intent(inout),target :: mesh
    character(*),intent(in) :: b_condition
    real(8),intent(in) :: w
    real(8),intent(in) :: cs(2)
    integer,intent(in) :: nmax

    type(bem3d_helmholtz),allocatable :: bem
    integer :: m, n
    type(incwave3d_helmholtz_spherical),allocatable :: incw

    self%mesh => mesh
    self%b_condition = b_condition

    self%w = w
    self%cs = cs

    self%k = self%w / self%cs(1)

    self%nmax = nmax

    self%size = (self%nmax+1)**2

    ! Neumann条件のみ対応
    call assert(self%b_condition == "neumann")

    !
    ! BEM
    !
    allocate(bem)
    call bem%init(self%mesh, self%b_condition, self%w, self%cs)

    ! 係数行列
    write(*,*) "# Constructing amat"
    call bem%gen_amat
    write(*,*) "# LU factorising"
    call bem%lu_decompose

    ! 各n,mについて右辺ベクトルを作り，解いて，その解を積分する
    allocate(incw)
    allocate(self%sols(self%size))
    allocate(self%S(self%size,self%size))
    do n=0,self%nmax
       do m=-n,n
          ! 入射波
          call incw%init(self%w, self%cs(1), one, n, m, self%mesh%center)
          ! 右辺ベクトル
          call bem%gen_bvec(incw)
          
          ! 解く
          write(*,*) "# LU solving", loct(n,m), "/", self%size
          call bem%lu_solve

          ! 解をself%solsの対応する場所に保持する
          call bem%gen_solution(self%sols(loct(n,m)))

          ! この解を境界積分してoutgoing coefficientsを計算する
          call self%sols(loct(n,m))%calc_B(self%nmax, self%S(:,loct(n,m)))
       end do
    end do
    
  end subroutine init
  
  !> 散乱体が球の場合について，解析解を用いてS行列を作る
  !! \param self
  !! \param mesh 球の境界mesh
  !! \param b_condition 境界条件 "neumann"
  !! \param w 角周波数
  !! \param cs 各領域の位相速度 1:外部 2:内部
  !! \param nmax 次数nの最大値
  subroutine init_sphere(self, mesh, b_condition, w, cs, nmax)
    class(smatrix_helmholtz),intent(out) :: self
    type(mesh3d),intent(inout),target :: mesh
    character(*),intent(in) :: b_condition
    real(8),intent(in) :: w
    real(8),intent(in) :: cs(2)
    integer,intent(in) :: nmax

    real(8) :: rad
    integer :: n, m, k, l
    integer :: i

    real(8),allocatable :: besj(:)
    complex(8),allocatable :: besh(:)

    complex(8),allocatable :: A(:), B(:)
    complex(8),allocatable :: Inm(:), Onm(:)

    real(8) :: x(3)
    
    self%mesh => mesh
    self%b_condition = b_condition

    self%w = w
    self%cs = cs

    self%k = self%w / self%cs(1)

    self%nmax = nmax

    self%size = (self%nmax+1)**2

    ! Neumann条件のみ対応
    call assert(self%b_condition == "neumann")

    ! 球の半径
    rad = mesh%radius

    ! 球Bessel,Hankel関数の計算
    allocate(besj(0:self%nmax+1))
    call bessel_sph_jn_array(self%nmax+1, self%k*rad, besj)
    allocate(besh(0:self%nmax+1))
    call bessel_sph_hn_array(self%nmax+1, self%k*rad, besh)
    
    ! Sの成分を計算
    allocate(self%S(self%size,self%size))
    self%S(:,:) = zero
    do n=0,self%nmax
       do m=-n,n
          self%S(loct(n,m), loct(n,m)) = -&
               (n/(self%k*rad)*besj(n) - besj(n+1)) / &
               (n/(self%k*rad)*besh(n) - besh(n+1))
       end do
    end do

    !
    ! 解を計算
    !
    allocate(self%sols(self%size))
    allocate(A(self%size), B(self%size))
    allocate(Inm(self%size), Onm(self%size))
    
    do n=0,self%nmax
       do m=-n,n          
          ! incoming coefficients Aは(n,m)だけ1で，他は0
          A(:) = zero
          A(loct(n,m)) = one
          ! outgoing coefficients B = S*A
          B = matmul(self%S,A)

          ! 初期化
          call self%sols(loct(n,m))%init(self%mesh, self%w, self%cs(1))
          self%sols(loct(n,m))%u_bndry(:,:) = zero
          ! 各点iで入射波と散乱波を足す
          do i=1,self%mesh%ne
             x(:) = mesh%c(:,i)
             
             ! IとOを計算
             call calc_Inm(self%k, x-self%mesh%center, self%nmax, Inm)
             call calc_Onm(self%k, x-self%mesh%center, self%nmax, Onm)

             ! 入射波
             do k=0,self%nmax
                do l=-k,k
                   ! u
                   self%sols(loct(n,m))%u_bndry(i,1) = self%sols(loct(n,m))%u_bndry(i,1) + A(loct(k,l))*Inm(loct(k,l))
                end do
             end do

             ! 散乱波
             do k=0,self%nmax
                do l=-k,k
                   ! u
                   self%sols(loct(n,m))%u_bndry(i,1) = self%sols(loct(n,m))%u_bndry(i,1) + B(loct(k,l))*Onm(loct(k,l))
                end do
             end do
          end do
          
       end do
    end do    
    
  end subroutine init_sphere
  
  !> incoming coefficients Aを与えたときの解solを計算する
  subroutine calc_sol(self, A, sol)
    class(smatrix_helmholtz),intent(in) :: self
    complex(8),intent(in) :: A(self%size)
    type(solution_helmholtz),intent(out) :: sol

    integer :: n, m

    call sol%init(self%mesh, self%w, self%cs(1))

    sol%u_bndry(:,:) = zero
    do n=0,self%nmax
       do m=-n,n
          sol%u_bndry(:,:) = sol%u_bndry(:,:) + A(loct(n,m))*self%sols(loct(n,m))%u_bndry(:,:)
       end do
    end do
  end subroutine calc_sol

  !> \f$ p\in S^2 \f$方向に進む平面波を表すincoming coefficients Aを計算
  !! \param self
  !! \param pvec 入射波の進行方向\f$ p\in S^2 \f$
  !! \param amp 入射波の複素振幅
  !! \param A incoming coefficients
  subroutine calc_A_planewave(self, pvec, amp, A)
    class(smatrix_helmholtz),intent(in) :: self
    real(8),intent(in) :: pvec(3)
    complex(8),intent(in) :: amp
    complex(8),intent(out) :: A(self%size)

    integer :: n, m
    complex(8),allocatable :: Ys(:)

    ! 球面調和関数Y(pvec)を計算
    allocate(Ys(self%size))
    call sph(self%nmax, pvec, Ys)

    do n=0,self%nmax
       do m=-n,n
          A(loct(n,m)) = ione**n*(2*n+1)*conjg(Ys(loct(n,m)))
       end do
    end do
    A = A * exp(ione*self%k*rdot_product(pvec,self%mesh%center))

    deallocate(Ys)
  end subroutine calc_A_planewave

  !> Sの成分をファイルに出力
  subroutine output_S(self, filename)
    class(smatrix_helmholtz),intent(in) :: self
    character(*),intent(in) :: filename

    integer :: n, m
    integer :: k, l

    open(10,file=filename,status="replace")
    ! 外部領域の波数
    write(10,*) self%k
    ! nmax
    write(10,*) self%nmax
    do k=0,self%nmax
       do l=-k,k
          do n=0,self%nmax
             do m=-n,n
                write(10,"(4i4,2e24.16)") k, l, n, m, self%S(loct(k,l),loct(n,m))
             end do
          end do
          ! 空行を入れる
          write(10,*) ""
       end do
    end do
    close(10)
  end subroutine output_S
end module smatrix_helmholtz_class
