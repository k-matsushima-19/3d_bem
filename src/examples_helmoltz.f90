!> 3次元Helmholtzに関連する例
module examples_helmholtz
  implicit none

contains
  !> 球による平面波散乱の解析解を計算
  subroutine multi_smatrix_ana
    use misc
    use file_io
    use math
    use mesh3d_class
    use solution_helmholtz_class
    use smatrix_helmholtz_class
    use smatrices_helmholtz_class

    real(8) :: w
    real(8) :: c
    real(8) :: k
    
    integer :: nmax

    type(mesh3d),allocatable :: mesh1, mesh2
    character(:),allocatable :: fn_mesh

    type(solution_helmholtz),allocatable :: sol
    type(smatrix_helmholtz),allocatable :: smat1, smat2
    type(smatrices_helmholtz),allocatable :: smats

    complex(8),allocatable :: A(:)

    real(8) :: pvec(3)

    ! tmp
    real(8) :: x_p(3), x_q(3), x(3)
    complex(8),allocatable :: O(:), B(:), TB(:), T(:,:), I(:)
    complex(8) :: ret
    integer :: n, m
    
    !
    ! read input
    !
    open(10, file="input.conf")
    call read_config(10, "w", w, "3.0 # angular freq.")
    call read_config(10, "c", c, "1.2 # phase velocity")
    call read_config(10, "nmax", nmax, "+20 # maximum index of spherical funcs >= 0")
    close(10)

    

    ! ! Eを計算
    ! call output_E("E.off", "E.dat", 10)
    ! stop

    ! 球1のメッシュを読み込み
    allocate(mesh1)
    call mesh1%read_off("sphere.off")
    call mesh1%invert
    ! mesh1をshift
    call mesh1%shift([-1.5d0,0.d0,0.d0])
    ! S行列を作る
    allocate(smat1)
    ! call smat1%init_sphere(mesh1, "neumann", w, [c,c], nmax)
    call smat1%init(mesh1, "neumann", w, [c,c], nmax)
    call smat1%output_S("smatrix1.dat")

    ! 球2のメッシュを読み込み
    allocate(mesh2)
    call mesh2%read_off("sphere.off")
    call mesh2%invert
    ! mesh2をshift
    call mesh2%shift([+1.5d0,0.d0,0.d0])
    ! S行列を作る
    allocate(smat2)
    ! call smat2%init_sphere(mesh2, "neumann", w, [c,c], nmax)
    call smat2%init(mesh2, "neumann", w, [c,c], nmax)
    call smat1%output_S("smatrix2.dat")
    
    ! 2つのS行列からsmatrices構造体を作る
    allocate(smats)
    call smats%init(2, [smat1, smat2])

    call smats%gen_amat
    call smats%lu_decompose
    pvec = [0.d0, 1.d0, 0.d0]
    call smats%gen_bvec_planewave(one, pvec)
    call smats%lu_solve

    allocate(sol)
    call smats%gen_solution(sol)

    call sol%output_bndry("bndry.dat")
    
  end subroutine multi_smatrix_ana
  
  !> 球による平面波散乱の解析解を計算
  subroutine compute_smatrix_ana
    use misc
    use file_io
    use math
    use mesh3d_class
    use solution_helmholtz_class
    use smatrix_helmholtz_class

    real(8) :: w
    real(8) :: c
    real(8) :: k
    
    integer :: nmax

    type(mesh3d),allocatable :: mesh
    character(:),allocatable :: fn_mesh

    type(solution_helmholtz),allocatable :: sol
    type(smatrix_helmholtz),allocatable :: smat

    complex(8),allocatable :: A(:)

    real(8) :: pvec(3)
    
    !
    ! read input
    !
    open(10, file="input.conf")
    call read_config(10, "fn_mesh", fn_mesh, "elm_init.off # メッシュのファイル名")
    call read_config(10, "w", w, "3.0 # angular freq.")
    call read_config(10, "c", c, "1.2 # phase velocity")
    call read_config(10, "nmax", nmax, "+20 # maximum index of spherical funcs >= 0")
    close(10)

    ! 球のメッシュを読み込み
    allocate(mesh)
    call mesh%read_off(fn_mesh)
    call mesh%invert

    allocate(smat)
    ! call smat%init_sphere(mesh, "neumann", w, [c,c], nmax)
    call smat%init(mesh, "neumann", w, [c,c], nmax)

    call smat%output_S("smatrix.dat")

    ! 例として，[0.d0, 1.d0, 0.d0]方向の平面波が進行するときの解を計算
    allocate(A(smat%size))
    pvec = [0.d0, 1.d0, 0.d0]
    call smat%calc_A_planewave(pvec, one, A)
    allocate(sol)
    call smat%calc_sol(A, sol)
    
    call sol%output_bndry("bndry.dat")
    call smat%sols(loct(0,0))%output_bndry("bndry_0_0.dat")
    call smat%sols(loct(2,1))%output_bndry("bndry_2_1.dat")

  end subroutine compute_smatrix_ana
  
  !> 球による平面波散乱の解析解を計算
  subroutine ana_sphere_planewave
    use misc
    use file_io
    use math
    use mesh3d_class
    use solution_helmholtz_class

    real(8) :: w
    real(8) :: c
    real(8) :: k
    real(8) :: pvec(3)

    integer :: n, m
    integer :: nmax, size, index
    complex(8),allocatable :: A(:), B(:)
    complex(8),allocatable :: Ys(:)
    real(8),allocatable :: besj(:)
    complex(8),allocatable :: besh(:)

    type(mesh3d),allocatable :: mesh
    character(:),allocatable :: fn_mesh

    real(8) :: rad
    real(8) :: r, theta, phi

    integer :: i
    real(8) :: x(3)
    complex(8) :: u

    type(solution_helmholtz),allocatable :: sol
    
    !
    ! read input
    !
    open(10, file="input.conf")
    call read_config(10, "fn_mesh", fn_mesh, "elm_init.off # メッシュのファイル名")
    call read_config(10, "w", w, "3.0 # angular freq.")
    call read_config(10, "c", c, "1.2 # phase velocity")
    call read_config(10, "nmax", nmax, "+20 # maximum index of spherical funcs >= 0")
    close(10)

    k = w / c
    pvec = [0.d0, 1.d0, 0.d0]

    ! 球のメッシュを読み込み
    allocate(mesh)
    call mesh%read_off(fn_mesh)
    call mesh%invert

    ! size = sum_0^nmax sum_{m=-n}^n
    size = (nmax+1)**2

    !
    ! Aを計算
    !
    allocate(A(size))
    ! 球面調和関数Y(pvec)を計算
    allocate(Ys(size))
    call sph(nmax, pvec, Ys)
    
    do n=0,nmax
       do m=-n,n
          index = loct(n,m)
          
          A(index) = ione**n*(2*n+1)*conjg(Ys(index))
       end do
    end do

    deallocate(Ys)

    rad = mesh%radius
    write(*,*) "Radius: ", rad

    !
    ! Bを計算
    !
    ! 球Bessel,Hankel関数の計算
    allocate(besj(0:nmax+1))
    call bessel_sph_jn_array(nmax+1, k*rad, besj)
    allocate(besh(0:nmax+1))
    call bessel_sph_hn_array(nmax+1, k*rad, besh)

    ! Neumann条件からBを求める
    allocate(B(size))
    do n=0,nmax
       do m=-n,n
          index = loct(n,m)

          B(index) = -&
               (n/(k*rad)*besj(n) - besj(n+1)) / &
               (n/(k*rad)*besh(n) - besh(n+1)) &
               *A(index)
               
       end do
    end do

    deallocate(besj)
    deallocate(besh)

    !
    ! 境界上の解を計算
    !
    allocate(sol)
    call sol%init(mesh, w, c)

    allocate(besh(0:nmax))    
    allocate(Ys(size))
    
    do i=1,mesh%ne
       x(:) = mesh%c(:,i)

       ! xを球面座標に変換
       call Cartesian_to_spherical(x, r, theta, phi)

       ! Hankel関数の計算
       call bessel_sph_hn_array(nmax, k*r, besh)

       ! 球面調和関数Y(x)を計算
       call sph(nmax, x, Ys)

       ! 散乱波
       sol%u_bndry(i,:) = zero
       do n=0,nmax
          do m=-n,n
             index = loct(n,m)

             sol%u_bndry(i,1) = sol%u_bndry(i,1) + B(index)*besh(n)*Ys(index)
          end do
       end do

       ! 入射波を足す
       sol%u_bndry(i,1) = sol%u_bndry(i,1) + exp(ione*k*rdot_product(pvec,x))
    end do

    ! 出力
    call sol%output_bndry("bndry.dat")
    call sol%output_bndry_geomview("bndry.off")

    ! r = 1.d0
    ! theta = 2.d0
    ! phi = 1.2d0
    ! x(:) = [r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)]

    
    
    ! u = exp(ione*k*rdot_product(pvec,x))
    ! write(*,*) u

    ! ! 球面調和関数Y(x)を計算
    ! call sph(nmax, x, Ys)
    ! ! 球Bessel関数の計算
    ! allocate(besj(0:nmax))
    ! call bessel_sph_jn_array(nmax, k*rad, besj)
    ! u = zero
    ! do n=0,nmax
    !    do m=-n,n
    !       index = loct(n,m)
          
    !       u = u + A(index)*besj(n)*Ys(index)
    !    end do
    ! end do
    ! write(*,*) u
    
  end subroutine ana_sphere_planewave

  !> \f$ I_n^m(x) \f$入射のHelmholtz方程式の散乱問題をBEMで解く．
  subroutine solve_helmholtz3d_spherical
    use math
    use file_io
    use mesh3d_class
    use incwave3d_helmholtz_spherical_class
    use bem3d_helmholtz_class
    use solution_helmholtz_class
    implicit none

    real(8) :: w
    real(8) :: c

    type(mesh3d),allocatable :: mesh
    type(incwave3d_helmholtz_spherical),allocatable :: incw
    type(bem3d_helmholtz),allocatable :: bem
    type(solution_helmholtz),allocatable :: sol

    character(:),allocatable :: fn_mesh

    integer :: n, m
    
    !
    ! read input
    !
    open(10, file="input.conf")
    call read_config(10, "fn_mesh", fn_mesh, "elm_init.off # メッシュのファイル名")
    call read_config(10, "w", w, "3.0 # angular freq.")
    call read_config(10, "c", c, "1.2 # phase velocity")
    call read_config(10, "n", n, "3 # order of I")
    call read_config(10, "m", m, "-1 # order of I")
    close(10)

    allocate(mesh)
    call mesh%read_off(fn_mesh)

    call mesh%invert
    call mesh%output_an("normal.dat")

    allocate(bem)
    call bem%init(mesh, "neumann", w, [c,c])

    ! 係数行列
    write(*,*) "# Constructing amat"
    call bem%gen_amat
    write(*,*) "# LU factorising"
    call bem%lu_decompose

    ! 右辺ベクトル
    allocate(incw)
    call incw%init(w, c, one, n, m, mesh%center)
    call bem%gen_bvec(incw)

    ! 解く
    write(*,*) "# LU solving"
    call bem%lu_solve

    allocate(sol)
    call bem%gen_solution(sol)
    call sol%output_bndry("bndry.dat")
  end subroutine solve_helmholtz3d_spherical

  !> 平面波入射のHelmholtz方程式の散乱問題をBEMで解く．
  subroutine solve_helmholtz3d_planewave
    use math
    use file_io
    use mesh3d_class
    use incwave3d_helmholtz_planewave_class
    use bem3d_helmholtz_class
    use solution_helmholtz_class
    implicit none

    real(8) :: w
    real(8) :: c

    type(mesh3d),allocatable :: mesh
    type(incwave3d_helmholtz_planewave),allocatable :: incw
    type(bem3d_helmholtz),allocatable :: bem
    type(solution_helmholtz),allocatable :: sol

    character(:),allocatable :: fn_mesh
    
    !
    ! read input
    !
    open(10, file="input.conf")
    call read_config(10, "fn_mesh", fn_mesh, "elm_init.off # メッシュのファイル名")
    call read_config(10, "w", w, "3.0 # angular freq.")
    call read_config(10, "c", c, "1.2 # phase velocity")
    close(10)

    allocate(mesh)
    call mesh%read_off(fn_mesh)

    call mesh%invert
    call mesh%output_an("normal.dat")

    allocate(bem)
    call bem%init(mesh, "neumann", w, [c,c])

    ! 係数行列
    write(*,*) "# Constructing amat"
    call bem%gen_amat
    write(*,*) "# LU factorising"
    call bem%lu_decompose

    ! 右辺ベクトル
    allocate(incw)
    call incw%init(w, c, one, [zero, one, zero])
    call bem%gen_bvec(incw)

    ! 解く
    write(*,*) "# LU solving"
    call bem%lu_solve

    allocate(sol)
    call bem%gen_solution(sol)
    call sol%output_bndry("bndry.dat")
  end subroutine solve_helmholtz3d_planewave
  
end module examples_helmholtz
