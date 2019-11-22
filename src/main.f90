!> mainのプログラム
program main
  use mesh3d_class
  implicit none

  type(mesh3d),allocatable :: mesh

  allocate(mesh)
  call mesh%read_off("sphere.off")
  call mesh%invert

  call mesh%output_an("normal.dat")

end program main

  
