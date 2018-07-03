!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
module mod_getrgnid
  implicit none 

  public   :: setup_targetrgn

  !private  ::

contains 
  !
  subroutine  setup_targetrgn(  &
                 gl, rl,           &
                 gall, lall,       &
                 edgeall,          &
                 GrdEdge,          &
                 tgdim,            &
                 tgrid,            &
                 gwflag,           &
                 RgnNumber         &
                 )

    implicit none

    integer, parameter    :: kind_real = 8
    integer, intent(in)   :: gl, rl
    integer, intent(in)   :: gall, lall
    integer, intent(in)   :: edgeall
    integer, intent(in)   :: tgdim
    real(kind_real), intent(in)   :: GrdEdge(lall,edgeall)
    real(kind_real), intent(in)   :: tgrid(tgdim)
    integer, intent(in)   :: gwflag(lall)
    integer, intent(out)  :: RgnNumber(lall)

    integer   :: irgn
    integer   :: iedge
    integer   :: flag

    integer   :: i,j,k

    real(kind_real)   :: ilon , ilat 
    real(kind_real), parameter   :: dlon = 1.0d-2
    real(kind_real), parameter   :: dlat = 1.0d-2
    real(kind_real), parameter   :: wlon = 5.0d0
    real(kind_real), parameter   :: wlat = 5.0d0

    !--------------------------------------------------------------------------

    ! replace target grid value to const
    ilon = tgrid(1)
    ilat = tgrid(2)

    RgnNumber(:) = -9999 ! undef value = -9999
    do irgn = 1, lall

      flag = gwflag(irgn)
      if ( flag == 0 ) then
        if ( ( ilon >= GrdEdge(irgn, 2)-dlon .and. ilon <= GrdEdge(irgn, 4)+dlon )  .and.  & 
             ( ilat >= GrdEdge(irgn, 1)-dlat .and. ilat <= GrdEdge(irgn, 3)+dlat)  ) then

              RgnNumber(irgn) = irgn-1
              flag = 2

        else if ( flag == 0 .and.  &
             ( ilon >= GrdEdge(irgn, 2)-wlon .and. ilon <= GrdEdge(irgn, 4)+wlon )  .and.  &
             ( ilat >= GrdEdge(irgn, 1)-wlat .and. ilat <= GrdEdge(irgn, 3)+wlat)  ) then

              RgnNumber(irgn) = irgn-1
              flag = 2
        !else
          !write(13,*) ilon, ilat
          !write(13,*) GrdEdge(1, 1),GrdEdge(1, 2),GrdEdge(1, 3),GrdEdge(1, 4)
          !write(13, *) '  '

        end if

      else if ( flag == 1 ) then ! over Greenwitch Line

        if ( ( ilon+360.d0 >= GrdEdge(irgn, 2)-dlon .and. ilon <= GrdEdge(irgn, 4)+dlon )  .and.  &
             ( ilat >= GrdEdge(irgn, 1)-dlat .and. ilat <= GrdEdge(irgn, 3)+dlat)  ) then

              RgnNumber(irgn) = irgn-1
              flag = 2

        else if ( flag == 1 .and. &
             ( ilon >= GrdEdge(irgn, 2)-dlon .and. ilon <= 360.d0 + GrdEdge(irgn, 4)+dlon )  .and.  &
             ( ilat >= GrdEdge(irgn, 1)-dlat .and. ilat <= GrdEdge(irgn, 3)+dlat)  ) then

              RgnNumber(irgn) = irgn-1
              flag = 2

        else if ( flag == 1 .and. &
             ( ilon+360.d0 >= GrdEdge(irgn, 2)-wlon .and. ilon <= GrdEdge(irgn, 4)+wlon )  .and.  &
             ( ilat >= GrdEdge(irgn, 1)-wlat .and. ilat <= GrdEdge(irgn, 3)+wlat)  ) then

              RgnNumber(irgn) = irgn-1
              flag = 2

        else if ( flag == 1 .and. &
             ( ilon >= GrdEdge(irgn, 2)-wlon .and. ilon <= 360.d0 + GrdEdge(irgn, 4)+wlon )  .and.  &
             ( ilat >= GrdEdge(irgn, 1)-wlat .and. ilat <= GrdEdge(irgn, 3)+wlat)  ) then

              RgnNumber(irgn) = irgn-1
              flag = 2
        end if

      end if

    end do ! irgn 

    return
  end subroutine 

end  module 
