!> mainのプログラム
program main  
  use examples_helmholtz
  implicit none

  ! call ana_sphere
  ! call solve_helmholtz3d_planewave
  ! call solve_helmholtz3d_spherical
  call compute_smatrix_ana
  
end program main

  
