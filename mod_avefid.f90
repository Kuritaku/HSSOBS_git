module mod_avefid
  implicit none

  public :: MISC_get_available_fid  !--- get an available file ID
  public :: MISC_get_little_fid  !--- get an available file ID
  !## share/mod_misc.f90

  !----minimum available file-id
  integer, parameter, private :: min_fid = 7

  !----maxmum available file-id
  ! fixme
  integer, parameter, private :: max_fid = 99

contains
  !
  function  MISC_get_available_fid()  &
        result(fid)           ! -- file id
    !
    implicit none
    !
    integer  :: fid
    logical  :: i_opened

    !------------------------------------------------------

    do fid = min_fid, max_fid
      INQUIRE( fid , opened=i_opened )
      if (.not. i_opened ) return
    end do
  end function MISC_get_available_fid
  !
  function  MISC_get_little_fid()  &
        result(fid)           ! -- file id
    !
    implicit none
    !
    integer  :: fid
    logical  :: i_opened

    !------------------------------------------------------
      !
      ! Kurihana's bashrc , 9999 is little-endian IO
      !
      fid = 9999
      !INQUIRE( fid , opened=i_opened )
      !if (.not. i_opened ) return
  end function MISC_get_little_fid

end module
