!> 3次元Helmholtzに関連する例
module examples_helmholtz
  implicit none

contains
  !> 球による平面波散乱の解析解を計算
  subroutine ana_sphere
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
    integer :: n_max, size, index
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
    call read_config(10, "n_max", n_max, "+20 # maximum index of spherical funcs >= 0")
    close(10)

    k = w / c
    pvec = [zero, one, zero]

    ! 球のメッシュを読み込み
    allocate(mesh)
    call mesh%read_off(fn_mesh)
    call mesh%invert

    ! size = sum_0^n_max sum_{m=-n}^n
    size = (n_max+1)**2

    !
    ! Aを計算
    !
    allocate(A(size))
    ! 球面調和関数Y(pvec)を計算
    allocate(Ys(size))
    call sph(n_max, pvec, Ys)
    
    do n=0,n_max
       do m=-n,n
          index = loct_sph(n,m)
          
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
    allocate(besj(0:n_max+1))
    call bessel_sph_jn_array(n_max+1, k*rad, besj)
    allocate(besh(0:n_max+1))
    call bessel_sph_hn_array(n_max+1, k*rad, besh)

    ! Neumann条件からBを求める
    allocate(B(size))
    do n=0,n_max
       do m=-n,n
          index = loct_sph(n,m)

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
    call sol%init(mesh)

    allocate(besh(0:n_max))    
    allocate(Ys(size))
    
    do i=1,mesh%ne
       x(:) = mesh%c(:,i)

       ! xを球面座標に変換
       call Cartesian_to_spherical(x, r, theta, phi)

       ! Hankel関数の計算
       call bessel_sph_hn_array(n_max, k*r, besh)

       ! 球面調和関数Y(x)を計算
       call sph(n_max, x, Ys)

       ! 散乱波
       sol%u_bndry(i,:) = zero
       do n=0,n_max
          do m=-n,n
             index = loct_sph(n,m)

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
    ! call sph(n_max, x, Ys)
    ! ! 球Bessel関数の計算
    ! allocate(besj(0:n_max))
    ! call bessel_sph_jn_array(n_max, k*rad, besj)
    ! u = zero
    ! do n=0,n_max
    !    do m=-n,n
    !       index = loct_sph(n,m)
          
    !       u = u + A(index)*besj(n)*Ys(index)
    !    end do
    ! end do
    ! write(*,*) u
    
  end subroutine ana_sphere

  subroutine solve_helmholtz3d
    use math
    use mesh3d_class
    use incwave3d_helmholtz_planewave_class
    use bem3d_helmholtz_class
    use solution_helmholtz_class
    implicit none

    real(8),parameter :: w = 3.0d0
    real(8),parameter :: c = 1.2d0

    type(mesh3d),allocatable :: mesh
    type(incwave3d_helmholtz_planewave),allocatable :: incw
    type(bem3d_helmholtz),allocatable :: bem
    type(solution_helmholtz),allocatable :: sol

    ! allocate(mesh)
    ! call mesh%read_off("sphere.off.old")
    ! call mesh%shift(-[6.6154999999999973d0, -313.21600000000001d0, 165.90299999999999d0])
    ! call mesh%expand(1.d0 / 170.79050000000001d0)
    ! call mesh%output_off("sphere.off")

    ! stop

    allocate(mesh)
    call mesh%read_off("sphere.off")

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
    call sol%output_bndry_geomview("bndry.off")
  end subroutine solve_helmholtz3d
end module examples_helmholtz
