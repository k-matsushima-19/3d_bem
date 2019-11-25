module file_io
  implicit none

  interface read_config
     module procedure read_config_r
     module procedure read_config_rs
     module procedure read_config_i
     module procedure read_config_c
     module procedure read_config_l
  end interface read_config
contains

  !read �ե�������ɤ߹����key�򸡺�����val�˼¿����Ǽ
  recursive subroutine read_config_r(f_n, key, val, default)
    implicit none
    integer,intent(in) ::f_n
    character(*),intent(in) :: key
    real(8),intent(out) :: val
    character(*),intent(in) :: default
    integer :: ios, index_key, index_sh
    integer,parameter :: max_length = 2000
    character(len=max_length) :: line
    logical :: has_key

    rewind(f_n)
    has_key = .false.

    ! ��Ԥ����ɤ߹���
    do
       read(f_n,'(a)',iostat=ios) line
       if(ios < 0) exit
       line = trim(adjustl(line))
       if(line(1:1) == "#") cycle ! �����ȹԤ����Ф�
       index_key = index(line, key)
       ! key�����Ĥ��ä����
       if(index_key == 1) then
          ! key�������
          line = line(index_key+len(key):len_trim(line))
          line = adjustl(line)
          ! �ǽ��ʸ����=�Ǥʤ���а㤦
          if(line(1:1) /= "=") cycle          
          line = line(2:len_trim(line)) ! =��ä�
          ! #������Ȥ��Ϥ������鱦�������
          index_sh = index(line, "#")          
          if(index_sh /= 0) line = line(1:index_sh-1)
          line = line(1:len_trim(line))
          line = trim(adjustl(line))
          read(line,*) val
          has_key = .true.
          exit
       endif
    enddo
    ! key�����Ĥ���ʤ��ä����ϥǥե���Ȥ�񤭤����ƺƵ�
    if(.not. has_key) then
       backspace(f_n)
       write(*,*) "# '" // key // "'" // " was not found in the config file."
       write(*,*) "# '" // key //" = " // default // "' is to be used."
       write(f_n,*) key//" = "//default
       call read_config_r(f_n, key, val, default)
    endif
  end subroutine read_config_r

  !read �ե�������ɤ߹����key�򸡺�����val�˼¿����Ǽ
  recursive subroutine read_config_rs(f_n, key, val, default)
    implicit none
    integer,intent(in) ::f_n
    character(*),intent(in) :: key
    real(8),intent(out) :: val(:)
    character(*),intent(in) :: default
    integer :: ios, index_key, index_sh
    integer,parameter :: max_length = 2000
    character(len=max_length) :: line
    logical :: has_key

    rewind(f_n)
    has_key = .false.

    ! ��Ԥ����ɤ߹���
    do
       read(f_n,'(a)',iostat=ios) line
       if(ios < 0) exit
       line = trim(adjustl(line))
       if(line(1:1) == "#") cycle ! �����ȹԤ����Ф�
       index_key = index(line, key)
       ! key�����Ĥ��ä����
       if(index_key == 1) then
          ! key�������
          line = line(index_key+len(key):len_trim(line))
          line = adjustl(line)
          ! �ǽ��ʸ����=�Ǥʤ���а㤦
          if(line(1:1) /= "=") cycle          
          line = line(2:len_trim(line)) ! =��ä�
          ! #������Ȥ��Ϥ������鱦�������
          index_sh = index(line, "#")          
          if(index_sh /= 0) line = line(1:index_sh-1)
          line = line(1:len_trim(line))
          line = trim(adjustl(line))
          read(line,*) val
          has_key = .true.
          exit
       endif
    enddo
    ! key�����Ĥ���ʤ��ä����ϥǥե���Ȥ�񤭤����ƺƵ�
    if(.not. has_key) then
       backspace(f_n)
       write(*,*) "# '" // key // "'" // " was not found in the config file."
       write(*,*) "# '" // key //" = " // default // "' is to be used."
       write(f_n,*) key//" = "//default
       call read_config_rs(f_n, key, val, default)
    endif
  end subroutine read_config_rs

  ! �ե�������ɤ߹����key�򸡺�����val���������Ǽ
  recursive subroutine read_config_i(f_n, key, val, default)
    implicit none
    integer,intent(in) ::f_n
    character(*),intent(in) :: key
    integer,intent(out) :: val
    character(*),intent(in) :: default
    integer :: ios, index_key, index_sh
    integer,parameter :: max_length = 2000
    character(len=max_length) :: line
    logical :: has_key

    rewind(f_n)
    has_key = .false.

    ! ��Ԥ����ɤ߹���
    do
       read(f_n,'(a)',iostat=ios) line
       if(ios < 0) exit
       line = trim(adjustl(line))
       if(line(1:1) == "#") cycle ! �����ȹԤ����Ф�
       index_key = index(line, key)
       ! key�����Ĥ��ä����
       if(index_key == 1) then
          ! key�������
          line = line(index_key+len(key):len_trim(line))
          line = adjustl(line)
          ! �ǽ��ʸ����=�Ǥʤ���а㤦
          if(line(1:1) /= "=") cycle          
          line = line(2:len_trim(line)) ! =��ä�
          ! #������Ȥ��Ϥ������鱦�������
          index_sh = index(line, "#")          
          if(index_sh /= 0) line = line(1:index_sh-1)
          line = line(1:len_trim(line))
          line = trim(adjustl(line))
          read(line,*) val
          has_key = .true.
       endif
    enddo
    ! key�����Ĥ���ʤ��ä����ϥǥե���Ȥ�񤭤����ƺƵ�
    if(.not. has_key) then
       backspace(f_n)
       write(*,*) "# '" // key // "'" // " was not found in the config file."
       write(*,*) "# '" // key //" = " // default // "' is to be used."
       write(f_n,*) key//" = "//default
       call read_config_i(f_n, key, val, default)
    endif
  end subroutine read_config_i

  ! �ե�������ɤ߹����key�򸡺�����val��ʸ������Ǽ
  ! require�Υǥե���Ȥ�true
  recursive subroutine read_config_c(f_n, key, val, default)
    implicit none
    integer,intent(in) ::f_n
    character(*),intent(in) :: key
    character(:),allocatable,intent(inout) :: val
    character(*),intent(in) :: default
    integer :: ios, index_key, index_sh
    integer,parameter :: max_length = 2000
    character(len=max_length) :: line
    logical :: has_key

    rewind(f_n)
    has_key = .false.

    ! ��Ԥ����ɤ߹���
    do
       read(f_n,'(a)',iostat=ios) line
       if(ios < 0) exit
       line = trim(adjustl(line))
       if(line(1:1) == "#") cycle ! �����ȹԤ����Ф�
       index_key = index(line, key)
       ! key�����Ĥ��ä����
       if(index_key == 1) then
          ! key�������
          line = line(index_key+len(key):len_trim(line))
          line = adjustl(line)
          ! �ǽ��ʸ����=�Ǥʤ���а㤦
          if(line(1:1) /= "=") cycle          
          line = line(2:len_trim(line)) ! =��ä�
          ! #������Ȥ��Ϥ������鱦�������
          index_sh = index(line, "#")          
          if(index_sh /= 0) line = line(1:index_sh-1)
          line = line(1:len_trim(line))
          line = trim(adjustl(line))
          !read(line,*) val
          allocate(character(len=len_trim(line))::val)
          val(:) = line(:)
          has_key = .true.
          exit
       endif
    enddo
    ! key�����Ĥ���ʤ��ä����ϥǥե���Ȥ�񤭤����ƺƵ�
    if(.not. has_key) then
       backspace(f_n)
       write(*,*) "# '" // key // "'" // " was not found in the config file."
       write(*,*) "# '" // key //" = " // default // "' is to be used."
       write(f_n,*) key//" = "//default
       call read_config_c(f_n, key, val, default)
    endif
  end subroutine read_config_c

  ! �ե�������ɤ߹����key�򸡺�����val�˿������Ǽ
  ! require�Υǥե���Ȥ�true
  recursive subroutine read_config_l(f_n, key, val, default)
    implicit none
    integer,intent(in) ::f_n
    character(*),intent(in) :: key
    logical,intent(out) :: val
    character(*),intent(in) :: default
    integer :: ios, index_key, index_sh
    integer,parameter :: max_length = 2000
    character(len=max_length) :: line, filename
    logical :: has_key

    rewind(f_n)
    has_key = .false.

    ! ��Ԥ����ɤ߹���
    do
       read(f_n,'(a)',iostat=ios) line
       if(ios < 0) exit
       line = trim(adjustl(line))
       if(line(1:1) == "#") cycle ! �����ȹԤ����Ф�
       index_key = index(line, key)
       ! key�����Ĥ��ä����
       if(index_key == 1) then
          ! key�������
          line = line(index_key+len(key):len_trim(line))
          line = adjustl(line)
          ! �ǽ��ʸ����=�Ǥʤ���а㤦
          if(line(1:1) /= "=") cycle          
          line = line(2:len_trim(line)) ! =��ä�
          ! #������Ȥ��Ϥ������鱦�������
          index_sh = index(line, "#")          
          if(index_sh /= 0) line = line(1:index_sh-1)
          line = line(1:len_trim(line))
          line = trim(adjustl(line))
          if(line == 'true' .or. line == 'True' .or. line == 't' .or. line == 'T') then
             val = .true.
          else if(line == 'false' .or. line == 'False' .or. line == 'f' .or. line == 'F') then
             val = .false.
          else
             inquire(f_n, name=filename)
             write(*,*) "Error: Invalid config file"
             write(*,*) "You must set "//"'"//key//"' "//"in "//trim(filename)
          endif
          has_key = .true.
          exit
       endif
    enddo
    ! key�����Ĥ���ʤ��ä����ϥǥե���Ȥ�񤭤����ƺƵ�
    if(.not. has_key) then
       backspace(f_n)
       write(*,*) "# '" // key // "'" // " was not found in the config file."
       write(*,*) "# '" // key //" = " // default // "' is to be used."
       write(f_n,*) key//" = "//default
       call read_config_l(f_n, key, val, default)
    endif
  end subroutine read_config_l

  subroutine progress()
    implicit none
    integer,save :: count, count_max
    integer(kind=4) :: per, n, i
    character(len=57) :: bar="???% |                                                  |"
    count = count + 1
    per = count*100/count_max
    !write(*,*) count, per
    n = per / 2
    write(unit=bar(1:3),fmt="(i3)") per
    do i=1, n
       bar(6+i:6+i)="="
    enddo
    ! print the progress bar.
    write(unit=6,fmt="(a1,a57)",advance="no") char(13),bar
    if (per/=100) then
       flush(unit=6)
    else
       write(unit=6,fmt=*)
       ! initialize
       count = 0
       bar = "???% |                                                  |"
    endif
    return
  contains
    subroutine init_progress(n)
      integer,intent(in) :: n

      count_max = n
      count = 0
    end subroutine init_progress
  end subroutine progress

end module file_io
