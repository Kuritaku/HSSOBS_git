!------------------------------------------------------------------------------
!
! ++ Document here
!     Get ico-zlayer number from  obs-height
!
! ++ History 
!
!     Version      Date      Editer               Comments
!  --------------------------------------------------------------------------
!       0.0    2018-05-16   T. Kurihana   obs-zheight to ico zlayer num
!       0.1    2018-05-16   T. Kurihana   test - clear
!
!
!------------------------------------------------------------------------------
module mod_setup_zlayer
  implicit none

  public :: setup_icozgrid, &
            setup_icozlayer

contains
  !
  subroutine setup_icozgrid(                &
                             zall,          & !--IN
                             zgridinfo,     & !--IN
                             ico_zgrid      & !--OUT
                           )

    ! ++ module
    use mod_hssconst

    implicit none

    integer,              intent(in)  ::  zall
    character(kind_char), intent(in)  ::  zgridinfo
    real(kind_real),      intent(out) ::  ico_zgrid(zall)

    !--------------------------------------------------------------------------

    !--------------------------------------------------------
    !++ Setup ico zgrid
    !--------------------------------------------------------
    open(1, file=trim(zgridinfo), form='formatted', action='read' )
      read(1, *) ico_zgrid
    close(1)
    !
    write(*,*) '      ####  Finish  :   Ico - Zlayer    #####    '

    return
  end subroutine

  !----------------------------------------------------------------------------

  subroutine setup_icozlayer(               &
                              zall,         &  !--IN
                              ico_zgrid,    &  !--IN
                              obs_height,   &  !--IN
                              ico_znum      &  !--OUT
                            )

    ! ++ module
    use mod_hssconst

    implicit none

    real(kind_real),      parameter   :: init_val = 99999.d0

    integer,              intent(in)  :: zall
    real(kind_real),      intent(in)  :: ico_zgrid(zall)  ! ico-z height  [m]
    real(kind_real),      intent(in)  :: obs_height       ! obs-z height  [m]
    integer,              intent(out) :: ico_znum         ! ico-grid znumber

    integer               :: iz
    integer               :: i, j
    integer               :: tmp_znum
    real(kind_real)       :: dheight(0:zall)
    !--------------------------------------------------------------------------

    !initialize
    dheight(0)  = init_val
    tmp_znum    = 1
    do iz = 1, zall
      dheight(iz)   = abs( ico_zgrid(iz) - obs_height   )
      !
      ! judge z-height
      if (  dheight(iz) < dheight(iz-1)   ) then
        tmp_znum = iz
      end if

    end do

    ico_znum = tmp_znum
    return
  end subroutine

  !----------------------------------------------------------------------------

  !subroutine setup_dpr_icozlayer(               &
  !                              )

  !  implicit none

  !end subroutine


  !----------------------------------------------------------------------------

  !subroutine bubble_sort(             &
  !                      )

  !  implicit none


  !end subroutine 
end module
