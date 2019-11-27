!> mainのプログラム
program main
  ! use mesh3d_class

  ! type(mesh3d) :: mesh1, mesh2

  ! call mesh1%read_off("sphere.off")
  ! call mesh2%read_off("sphere.off")

  ! call mesh1%shift([-1.5d0, 0.d0, 0.d0])
  ! call mesh2%shift([+1.5d0, 0.d0, 0.d0])

  ! call mesh1%add_mesh(mesh2)
  ! call mesh1%output_off("double.off")
  
  use examples_helmholtz
  implicit none

  ! call ana_sphere
  ! call solve_helmholtz3d_planewave
  ! call solve_helmholtz3d_spherical
  ! call compute_smatrix_ana
  call multi_smatrix_ana
  
end program main

  
