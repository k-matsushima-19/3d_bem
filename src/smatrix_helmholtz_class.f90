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
  public :: output_E

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

     !> \f$ E_{n n^\prime n^{\prime\prime}}^{m m^\prime m^{\prime\prime}} \f$の次数\f$ n,n^\prime,n^{\prime\prime} \f$の最大値
     integer :: nmax_E
     !> \f$ E_{n n^\prime n^{\prime\prime}}^{m m^\prime m^{\prime\prime}} \f$の値
     complex(8),allocatable :: E(:,:,:)
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
     !> \f$ O_n^m \f$のtranslation matrix \f$ T_{nn^\prime}^{mm^\prime} \f$を計算
     procedure :: calc_T
     !> \f$ I_n^m \f$のtranslation matrix \f$ \tilde{T}_{nn^\prime}^{mm^\prime} \f$を計算
     procedure :: calc_T_tilde
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
    integer :: n, nn, nnn
    integer :: m, mm, mmm
    type(incwave3d_helmholtz_spherical),allocatable :: incw

    real(8) :: tmp1, tmp2

    self%mesh => mesh
    self%b_condition = b_condition

    self%w = w
    self%cs = cs

    self%k = self%w / self%cs(1)

    self%nmax = nmax

    self%size = (self%nmax+1)**2

    ! Neumann条件のみ対応
    call assert(self%b_condition == "neumann")

    ! Eをファイルから読み込む
    write(*,*) "# Reading E from E.dat"
    open(10,file="E.dat",status="old")
    read(10,*) self%nmax_E
    call assert(self%nmax_E >= self%nmax)
    allocate(self%E(self%size,self%size,self%size))
    do nnn=0,self%nmax_E
       do mmm=-nnn,nnn
          do nn=0,self%nmax_E
             do mm=-nn,nn
                do n=0,self%nmax_E
                   do m=-n,n
                      read(10,'(2e24.16)') tmp1, tmp2
                      ! 必要な次数の分だけ格納
                      if(n <= self%nmax) then
                         self%E(loct(n,m),loct(nn,mm),loct(nnn,mmm)) = dcmplx(tmp1,tmp2)
                      end if
                      
                   end do
                end do
             end do
          end do
       end do
    end do
    close(10)

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
    integer :: k, l
    integer :: i

    real(8),allocatable :: besj(:)
    complex(8),allocatable :: besh(:)

    complex(8),allocatable :: A(:), B(:)
    complex(8),allocatable :: Inm(:), Onm(:)

    real(8) :: x(3)

    integer :: n, nn, nnn
    integer :: m, mm, mmm
    real(8) :: tmp1, tmp2
    
    self%mesh => mesh
    self%b_condition = b_condition

    self%w = w
    self%cs = cs

    self%k = self%w / self%cs(1)

    self%nmax = nmax

    self%size = (self%nmax+1)**2

    ! Neumann条件のみ対応
    call assert(self%b_condition == "neumann")

    ! Eをファイルから読み込む
    write(*,*) "# Reading E from E.dat"
    open(10,file="E.dat",status="old")
    read(10,*) self%nmax_E
    call assert(self%nmax_E >= self%nmax)
    allocate(self%E(self%size,self%size,self%size))
    do nnn=0,self%nmax_E
       do mmm=-nnn,nnn
          do nn=0,self%nmax_E
             do mm=-nn,nn
                do n=0,self%nmax_E
                   do m=-n,n
                      read(10,'(2e24.16)') tmp1, tmp2
                      ! 必要な次数の分だけ格納
                      if(n <= self%nmax) then
                         self%E(loct(n,m),loct(nn,mm),loct(nnn,mmm)) = dcmplx(tmp1,tmp2)
                      end if
                      
                   end do
                end do
             end do
          end do
       end do
    end do
    close(10)

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
    ! 原点
    write(10,*) self%mesh%center
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

  !> \f$ E_{n n^\prime n^{\prime\prime}}^{m m^\prime m^{\prime\prime}} \f$を計算してファイルに出力する．
  !! \param filename_mesh 単位球の境界メッシュ (数値積分用) の.offファイル
  !! \param filename 出力ファイル
  !! \param nmax 次数\f$ n,n^\prime,n^{\prime\prime} \f$の最大値
  subroutine output_E(filename_mesh, filename, nmax)
    character(*),intent(in) :: filename_mesh
    character(*),intent(in) :: filename
    integer,intent(in) :: nmax

    integer :: i, ig
    integer :: n, nn, nnn
    integer :: m, mm, mmm
    complex(8),allocatable :: integ(:,:,:), E(:,:,:)
    integer :: size
    
    type(mesh3d),allocatable :: sphere

    ! 三角形上のGauss積分の積分点の数
    integer :: ng_tri
    ! 三角形上のGauss積分の積分点の面積座標 (3,ng_tri)
    real(8),allocatable :: gauss_tri_xi(:,:)
    ! 三角形上のGauss積分の積分点の重み (ng_tri)
    real(8),allocatable :: gauss_tri_w(:)

    real(8) :: y_vec(3)
    complex(8),allocatable :: ynm(:)

    ! 三角形Gauss積分の点数はとりあえず6とする
    ng_tri = 6
    allocate(gauss_tri_xi(3,ng_tri), gauss_tri_w(ng_tri))
    call generate_gauss_triangle(ng_tri, gauss_tri_xi, gauss_tri_w)

    allocate(sphere)
    call sphere%read_off(filename_mesh)
    

    ! 各要素についてGauss積分し，integに代入
    size = (nmax+1)**2
    allocate(integ(size,size,size))
    allocate(E(size,size,size))
    allocate(ynm(size))
    E(:,:,:) = zero
    do i=1,sphere%ne
       write(*,*) i, "/", sphere%ne
       
       ! Gauss積分
       integ(:,:,:) = zero
       do ig=1,ng_tri
          ! 積分点 (Cartesian)
          y_vec = convert_area_to_Cartesian(&
               sphere%p(:,sphere%nd(1,i)),&
               sphere%p(:,sphere%nd(2,i)),&
               sphere%p(:,sphere%nd(3,i)),&
               gauss_tri_xi(:,ig))

          ! 球面調和関数Yを計算
          call sph(nmax, y_vec, ynm)
          ! 各(n,m),(nn,mm),(nnn,mmm)について被積分関数を計算
          do nnn=0,nmax
             do mmm=-nnn,nnn
                do nn=0,nmax
                   do mm=-nn,nn
                      do n=0,nmax
                         do m=-n,n
                            integ(loct(n,m),loct(nn,mm),loct(nnn,mmm)) = integ(loct(n,m),loct(nn,mm),loct(nnn,mmm)) + &
                                 gauss_tri_w(ig) * ynm(loct(n,m)) * ynm(loct(nn,-mm)) * ynm(loct(nnn,-mmm))
                                 
                         end do
                      end do
                   end do
                end do
             end do
          end do                         
       end do
       ! integの計算終わり

       ! 係数を掛けてEに足す
       do nnn=0,nmax
          do mmm=-nnn,nnn
             do nn=0,nmax
                do mm=-nn,nn
                   do n=0,nmax
                      do m=-n,n
                         E(loct(n,m),loct(nn,mm),loct(nnn,mmm)) = E(loct(n,m),loct(nn,mm),loct(nnn,mmm)) + &
                              sphere%ar(i)*integ(loct(n,m),loct(nn,mm),loct(nnn,mmm)) &
                              * ione**(nn+nnn-n)/(4.d0*pi)*(-1)**(mm+mmm)
                      end do
                   end do
                end do
             end do
          end do
       end do
                 
    end do

    ! ファイルに出力
    open(10,file=filename)
    write(10,*) nmax
    do nnn=0,nmax
       do mmm=-nnn,nnn
          do nn=0,nmax
             do mm=-nn,nn
                do n=0,nmax
                   do m=-n,n
                      write(10,'(2e24.16)') E(loct(n,m),loct(nn,mm),loct(nnn,mmm))
                   end do
                end do
             end do
          end do
       end do
    end do

    close(10)
    
  end subroutine output_E

  !> \f$ O_n^m \f$のtranslation matrix \f$ T_{nn^\prime}^{mm^\prime} \f$を計算
  !! \param self
  !! \param vec \f$ T_{nn^\prime}^{mm^\prime} \f$の引数
  !! \param out \f$ T_{nn^\prime}^{mm^\prime} \f$の値を格納する配列
  subroutine calc_T(self, vec, out)
    class(smatrix_helmholtz),intent(in) :: self
    real(8),intent(in) :: vec(3)
    complex(8),intent(out) :: out(self%size,self%size)

    integer :: n, nn, nnn
    integer :: m, mm, mmm
    complex(8),allocatable :: Onm(:)

    allocate(Onm(self%size))

    ! Oを計算
    call calc_Onm(self%k, vec, self%nmax, Onm)
    
    ! 各(n,m), (nn,mm)について，Tの成分を計算
    out(:,:) = zero
    do nnn=0,self%nmax
       do mmm=-nnn,nnn          
          
          do nn=0,self%nmax
             do mm=-nn,nn
                do n=0,self%nmax
                   do m=-n,n
                      out(loct(n,m),loct(nn,mm)) = out(loct(n,m),loct(nn,mm)) + &
                           (2*n+1)*(2*nnn+1)*self%E(loct(nn,mm),loct(n,m),loct(nnn,mmm))*Onm(loct(nnn,mmm))
                      ! out(loct(n,m),loct(nn,mm)) = out(loct(n,m),loct(nn,mm)) + &
                      !      (2*nn+1)*(2*nnn+1)*self%E(loct(n,m),loct(nn,mm),loct(nnn,mmm))*Onm(loct(nnn,mmm))
                   end do
                end do
             end do
          end do

       end do
    end do
    
  end subroutine calc_T

  !> \f$ I_n^m \f$のtranslation matrix \f$ \tilde{T}_{nn^\prime}^{mm^\prime} \f$を計算
  !! \param self
  !! \param vec \f$ T_{nn^\prime}^{mm^\prime} \f$の引数
  !! \param out \f$ T_{nn^\prime}^{mm^\prime} \f$の値を格納する配列
  subroutine calc_T_tilde(self, vec, out)
    class(smatrix_helmholtz),intent(in) :: self
    real(8),intent(in) :: vec(3)
    complex(8),intent(out) :: out(self%size,self%size)

    integer :: n, nn, nnn
    integer :: m, mm, mmm
    complex(8),allocatable :: Inm(:)

    allocate(Inm(self%size))

    ! Iを計算
    call calc_Inm(self%k, vec, self%nmax, Inm)
    
    ! 各(n,m), (nn,mm)について，Tの成分を計算
    out(:,:) = zero
    do nnn=0,self%nmax
       do mmm=-nnn,nnn          
          
          do nn=0,self%nmax
             do mm=-nn,nn
                do n=0,self%nmax
                   do m=-n,n
                      out(loct(n,m),loct(nn,mm)) = out(loct(n,m),loct(nn,mm)) + &
                           (2*n+1)*(2*nnn+1)*self%E(loct(nn,mm),loct(n,m),loct(nnn,mmm))*Inm(loct(nnn,mmm))
                      ! out(loct(n,m),loct(nn,mm)) = out(loct(n,m),loct(nn,mm)) + &
                      !      (2*nn+1)*(2*nnn+1)*self%E(loct(n,m),loct(nn,mm),loct(nnn,mmm))*Inm(loct(nnn,mmm))
                   end do
                end do
             end do
          end do

       end do
    end do
    
  end subroutine calc_T_tilde
  
end module smatrix_helmholtz_class
