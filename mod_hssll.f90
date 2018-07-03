!------------------------------------------------------------------------------
!
!++ Document
!
!++  Hisotry
!
!       Version    Date       Editer        Description
!    ----------------------------------------------------------------
!        0.0      2018-04-06  T. Kurihana   Modulize Latlon Sel. Part
!
!
!++  API
!
!
!------------------------------------------------------------------------------
module mod_hssll

  public   :: setup_latlontab

  !private

contains
  !
  subroutine  setup_latlontab(                &
                               gl,            & ! IN
                               rl,            & ! IN
                               gall,          & ! IN
                               lall,          & ! IN
                               ctlon,         & ! IN
                               ctlat,         & ! IN
                               opt_pass,      & ! IN
                               RgnCandMaxID,  & ! IN
                               rgnid_dim,     & ! OUT
                               RgnCandArray   & ! OUT
                             )

    use  mod_llselection,  only :  &  ! set latlon thresval
         latlon_selection

    use  mod_getrgnid,     only :  &  ! set rgnid
         setup_targetrgn

    implicit none
  
    integer,         parameter   :: kind_real  = 8
    integer,         parameter   :: edgeall    = 5
    integer,         parameter   :: tgdim      = 2 ! dim. of target grid
    integer,         parameter   :: grd_Xdim   = 1 ! dim.ID of longitude
    integer,         parameter   :: grd_Ydim   = 2 ! dim.ID of longitude

    integer, intent(in)         :: gl
    integer, intent(in)         :: rl
    integer, intent(in)         :: gall
    integer, intent(in)         :: lall
    integer, intent(in)         :: RgnCandMaxID

    real(kind_real), intent(in)  :: ctlon
    real(kind_real), intent(in)  :: ctlat

    logical, intent(in)          :: opt_pass

    !integer, intent(inout)         :: rgnid_dim
    integer, intent(out)         :: rgnid_dim
    !integer, allocatable,  intent(inout)   :: RgnCandArray(:)
    integer, intent(out)   :: RgnCandArray(RgnCandMaxID )

    integer, allocatable, save   :: gwflag(:)

    integer,         allocatable :: RgnNumber(:) ! tmp array to judge rgnid
    real(kind_real), allocatable :: GrdEdge(:,:)

    real(kind_real)              :: tmp_grid(tgdim)

    integer                      :: icount
    integer                      :: l
    !--------------------------------------------------------------------------

    if (  opt_pass  ) then
      ! Initialize
      !RgnCandArray(:) = 0
      RgnCandArray(:) = -1

      ! set temporaraly grid values
      tmp_grid(grd_Xdim) = ctlon
      tmp_grid(grd_Ydim) = ctlat

      ! set Grid Edge values in all rgn
      allocate( GrdEdge( lall , edgeall )  )
      allocate( gwflag(lall)  )
      call  latlon_selection(           &
                            gl,         & ! IN
                            rl,         & ! IN
                            gall,       & ! IN
                            lall,       & ! IN
                            GrdEdge,    & ! OUT
                            gwflag      & ! INOUT
                            )


      ! set rgnid
      allocate( RgnNumber( lall ) )
      call   setup_targetrgn(             &
                             gl,          & 
                             rl,          &
                             gall,        &
                             lall,        &
                             edgeall,     &
                             GrdEdge,     &
                             tgdim,       &
                             tmp_grid,    &
                             gwflag,      &
                             RgnNumber    &
                            )

      rgnid_dim = count( mask = RgnNumber >= 0) ! number of valid rgn
      !allocate( RgnCandArray(rgnid_dim) )
      icount = 1
      do l = 1, lall
         if (  RgnNumber(l) >=  0 .and. icount <= rgnid_dim ) then
           RgnCandArray( icount ) = RgnNumber(l)
           icount = icount + 1
         end if
      end do
      if ( sum(RgnNumber) == -9999*lall  ) then
           write(6,*)  ' searching fail '
           write(200,"(2f12.7)")  ctlon , ctlat
           !write(20, * )  ctlon , ctlat
      end if

      write(6, * ) '  Number of potential rgn  ', rgnid_dim
      write(6,*)  '   Rgn Number  '
      write(6, * ) RgnCandArray

      deallocate(GrdEdge , RgnNumber,  gwflag  )

      write(200,*) ' finish search  '

      return

      ! ======> until here search potential rgn-grid number
      !------------------------------------------------------------------------

    else if (  .not. opt_pass ) then
       !
       !++ Same Rgn grids  as Prior Process
       !
       return
       !
    else
       !
       !++ Exception
       !
       write( * ,*) '  Stop LL loop'
       stop
       !
     end if

  end subroutine setup_latlontab

  !----------------------------------------------------------------------------

end module
