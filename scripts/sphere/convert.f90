! merge.fが出力したファイルを.offに整形
program convert
  implicit none

  integer :: ne, np, itmp, i
  character(100) :: filename

  integer,allocatable :: nd(:,:)
  real(8),allocatable :: p(:,:)

  read(*,"(a)") filename
  open(10,file=trim(adjustl(filename)), status="old")

  read(10,*) ne, np, itmp
  ! nd
  allocate(nd(3,ne))
  do i=1,ne
     read(10,*) itmp, nd(1:3,i), itmp, itmp
  end do
  ! p
  allocate(p(3,np))
  do i=1,np
     read(10,*) itmp, p(1:3,i)
  end do
  close(10)

  ! 上書き
  open(10,file=trim(adjustl(filename)))
  write(10,'(A)') "OFF"
  write(10,'(3i7)') np, ne, 0
  do i=1,np
     write(10,'(3e24.16)') p(1,i),p(2,i),p(3,i)
  end do

  do i=1,ne              
     write(10,'(i2,3i7)') 3,nd(1,i)-1,nd(2,i)-1,nd(3,i)-1
  end do
  close(10)
  
end program convert
