!------------------------------------------------------------------------------
!
! ++ Document here
!     Get latlon info. (degree) from ico-grid (radian) 
!
! ++ History 
!     0.0  2018  03  13th   : T. Kurihana
!     0.1  2018  03  21th   : T. Kurihana ; I/O & Array pass test are cleared
!
!
!------------------------------------------------------------------------------
module mod_ico2latlon
  implicit none 

  public  :: setup_ico2ll , &
             MISC_get_latlon

  !private :: MISC_get_latlon !# /share/mod_misc.f90

contains
  !----------------------------------------------------------------------------
  subroutine setup_ico2ll(    &
             gl  , rl ,       & ! IN
             gall, lall ,     & ! IN
             rgnID,           & ! IN
             inputdir,        & ! IN
             !inputname,       & ! IN ! assume grid.rgnXXXXX
             outputdir,       & ! IN
             outputname,      & ! IN
             offline,         & ! IN
             opt_bin,         & ! IN
             lat, lon         & ! OUT
                         )
    !
    use mod_setupfname, only : &
        get_ifilename

    use mod_avefid,     only : &
        MISC_get_available_fid

    !
    implicit none
    integer, parameter    :: kind_real = 8

    integer,        intent(in)   :: gl, rl
    integer,        intent(in)   :: gall, lall
    integer,        intent(in)   :: rgnID
    character(256), intent(in)   :: inputdir
    character(256), intent(in)   :: outputdir
    character(256), intent(in)   :: outputname
    logical,        intent(in)   :: offline ! offline true => gen latlon files in dir 
    logical,        intent(in)   :: opt_bin ! output option  true : binary / false : txt

    real(kind_real), intent(out) :: lat(gall)
    real(kind_real), intent(out) :: lon(gall)

    character(128)               :: filenumber
    character(40)                :: prefix
    real(kind_real), allocatable :: atm(:,:) 

    real(kind_real)              :: pi
    real(kind_real)              :: coef_deg
    real(kind_real)              :: lat_val, lon_val
    real(kind_real)              :: x, y, z

    integer                      :: log_ID    = 200
    integer                      :: gall_input
    integer                      :: iogall, iolall
    integer                      :: ilat, ilon
    integer                      :: grd_Xdim, grd_Ydim, grd_Zdim
    integer                      :: igall

    ! set coefficients
    pi = 4.d0 * datan(1.d0)
    coef_deg = 180.d0/pi

    grd_Xdim = 1
    grd_Ydim = 2
    grd_Zdim = 3
    !--------------------------------------------------------------------------

    !----------------------------------
    ! ++ check gall & lall
    !----------------------------------
    iogall = ( 2**(gl - rl)+2 )**2
    iolall = ( 2**rl )**2 * 10
    if (  iogall /= gall .or. iolall /= lall ) then
        write( 6, * ) '   ! Inconsistent in Glevel or Rlevel settings  ! '
        write( 6, * ) '   Error at mod_ico2latlon.f90 [subroutine | setup_ico2ll ] '
        stop
    end if

    !---------------------------------
    ! ++ Open & Read Hgrid file
    !---------------------------------
    allocate(  atm( gall, 3))

    call get_ifilename(     &
               rgnID,       & !IN
               filenumber,  &
               prefix       &
                      )


    log_ID = MISC_get_available_fid()
    open(log_ID,file=trim(inputdir)//'/'//'grid.rgn'//trim(prefix)//trim(filenumber), &
                form='unformatted', access='sequential')

        read(log_ID) gall_input
        read(log_ID) atm( : , grd_Xdim)
        read(log_ID) atm( : , grd_Ydim)
        read(log_ID) atm( : , grd_Zdim)
    close(log_ID)

    !--------------------------------------------
    ! ++  Convert Ico to LatLon
    !  &  Store or Pass  LatLon Info.
    !--------------------------------------------
    if ( offline ) then

      log_ID = MISC_get_available_fid()
      if (    opt_bin   ) then
        open(log_ID, & 
             file=trim(outputdir)//'/bin_'//trim(outputname)//'.rgn'//trim(prefix)//trim(filenumber),&
             form='unformatted' , access='direct', &
             recl=kind_real*gall,  action='write')


        do igall = 1, gall

          x = atm(igall, grd_Xdim)
          y = atm(igall, grd_Ydim)
          z = atm(igall, grd_Zdim)

          call MISC_get_latlon(   &
                       x, y, z,   &  ! IN
                       lat_val,   &  ! OUT
                       lon_val    &  ! OUT
                       )

          lon(igall) = lon_val * coef_deg
          lat(igall) = lat_val * coef_deg
        end do 

        write(log_ID , rec=1 ) lon( : )
        write(log_ID , rec=2 ) lat( : )

        close(log_ID)


      else if ( opt_bin .ne. .true. ) then
        log_ID = MISC_get_available_fid()
        open(log_ID, & 
             file=trim(outputdir)//'/'//trim(outputname)//'.rgn'//trim(prefix)//trim(filenumber),&
             action='write')

        do igall = 1, gall

          x = atm(igall, grd_Xdim)
          y = atm(igall, grd_Ydim)
          z = atm(igall, grd_Zdim)

          call MISC_get_latlon(   &
                       x, y, z,   &  ! IN
                       lat_val,   &  ! OUT
                       lon_val    &  ! OUT
                       )

          write(log_ID, * ) lon_val*coef_deg, lat_val*coef_deg

        end do

        close(log_ID)

     else
       !
       ! ++ Exception Processing
       !
       write( 6, *)  '  !   Error Occurred in the I/O Output process  ! '
       write( 6, *)  '  Error Detection ===>  [ Suborutine | setup_ico2ll ]  '
       write( 6, *)  '      Stop at Program : mod_ico2latlon.f90        '
       stop

     end if

    else ! online

        do igall = 1, gall

          x = atm(igall, grd_Xdim)
          y = atm(igall, grd_Ydim)
          z = atm(igall, grd_Zdim)

          call MISC_get_latlon(   &
                       x, y, z,   &  ! IN
                       lat_val,   &  ! OUT
                       lon_val    &  ! OUT
                       )

          lon(igall) = lon_val*coef_deg
          lat(igall) = lat_val*coef_deg

        end do

    end if 

    deallocate( atm )

    return

  end subroutine

  !-----------------------------------------------------------------------------

  
  subroutine MISC_get_latlon( &
       x, y, z,               &  !--- IN  : Cartesian coordinate ( on the sphere )
       lat, lon               & !--- OUT : latitude and longitude
                            )
    !
    implicit none

    ! IN/OUT variables
    real(8),intent(in) :: x,y,z
    real(8),intent(out) :: lat, lon

    ! local variables
    real(8), parameter :: epsilon = 1.0D-99
    real(8) :: leng,leng_xy
    !
    leng=dsqrt(x*x+y*y+z*z)
    !
    ! --- vector length equals to zero.
    if(leng<epsilon) then
       lat=0.0D0
       lon=0.0D0
       return
    endif
    ! --- vector is parallele to z axis.
    if(z/leng>=1.0D0) then
       lat=dasin(1.0D0)
       lon=0.0D0
       return
    else if(z/leng<=-1.0D0) then
       lat=dasin(-1.0D0)
       lon=0.0D0
       return
    endif
    ! --- not parallele to z axis
    lat=dasin(z/leng)
    !
    leng_xy=dsqrt(x*x+y*y)
    if(leng_xy<epsilon) then
       lon=0.0D0
       return
    endif
    if(x/leng_xy>=1.0D0) then
       lon=dacos(1.0D0)
       if(y<0.0D0) lon=-lon
       return
    else if(x/leng_xy<=-1.0D0) then
       lon=dacos(-1.0D0)
       if(y<0.0D0) lon=-lon
       return
    endif
    lon=acos(x/leng_xy)
    if(y<0.0D0) lon=-lon
    return
  end subroutine MISC_get_latlon


end module 
