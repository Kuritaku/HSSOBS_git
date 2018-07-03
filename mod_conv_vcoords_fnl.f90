module  mod_conv_vcoords_fnl
  implicit none

  public :: get_zval ,              &
            get_pval ,              &
            get_fnltemp_boundary,   &
            deallocate_fnltemp,     &
            get_fnllatlon_ID

  ! variables in function
  integer,          parameter, private  :: kind_dble = 8
  real(kind_dble),  parameter, private  :: g = 9.80665
  real(kind_dble),  parameter, private  :: gamma_d = 0.0065
  real(kind_dble),  parameter, private  :: R_d = 287.d0
  real(kind_dble),  parameter, private  :: ps = 1000.d0
  !                                        avoide minus Z as possible

  ! ####  1000hPa temepature grid
  !
  ! preprocessing bash script only decode at 1000hPa temepature
  real(kind_dble),  allocatable, save, public   :: fnl_temp( : , : )
  real(kind_dble),  allocatable, save, public   :: fnl_lon( : )
  real(kind_dble),  allocatable, save, public   :: fnl_lat( : )

  ! common number of lon-lat
  ! shared vars for avoiding many IO
  integer, save, public   :: fnl_nlon
  integer, save, public   :: fnl_nlat

contains
  !
  function  get_zval( plev, T0 ) &
      result(zval)
    !
    implicit none
    !
    real(kind_dble)   :: plev     ! IN
    real(kind_dble)   :: T0       ! IN
    real(kind_dble)   :: zval     ! OUT
    real(kind_dble)   :: coef

    coef = R_d * gamma_d / g
    zval = T0 / gamma_d * ( 1.d0 - (plev/ps)**(coef)  )

  end function get_zval

  !--------------------------------------------------------------------
  !
  function  get_pval( zlev, T0 ) &
      result(pval)
    !
    implicit none
    !
    real(kind_dble)   :: zlev     ! IN
    real(kind_dble)   :: T0       ! IN
    real(kind_dble)   :: pval     ! OUT
    real(kind_dble)   :: coef

    coef = g / ( R_d * gamma_d )
    pval = ps * ( ( (T0 - gamma_d * zlev)/T0 )**coef )

  end function get_pval

  !--------------------------------------------------------------------


  subroutine get_fnltemp_boundary()
    !
    ! ++ Documentation
    !     Return ncep fnl reanalysis 1000hPa  temperature
    !
    ! ++  Memo
    !     FNL Data Structure [nx][ny] = [0:359][-90:90]
    !
    use mod_hssconst

    use mod_avefid,   only :  &
        MISC_get_available_fid

    implicit none

    integer                         :: fid
    integer                         :: ierr
    integer                         :: irec
    integer                         :: nlon ! fnl-reanalysis nlon
    integer                         :: nlat ! fnl-reanalysis nlat
    integer                         :: ilon, ilat
    real(kind_sngl), allocatable    :: temp(:,:)
    real(kind_real), allocatable    :: tmp_latconv(:,:)

    character(kind_char)            :: fnl_tempfiledir
    character(kind_char)            :: fnl_tempfile

    logical                         :: opt_west2east
    logical                         :: opt_south2north

    namelist / p2z_param / &
      nlon, nlat,          &
      fnl_tempfiledir,     &
      fnl_tempfile

    !------------------------------------------------------------------

    ! Initalize
    opt_west2east   = .true.
    opt_south2north = .true.

    ! open namelist
    open(1, file='hssdriver.cnf', iostat=ierr)
      read(1, nml=p2z_param)
    close(1)
    if ( ierr /= 0) then
        write(*,*) '      #### Inappropriate variable in namelist ####'
        write(*,*) '      [MOD | mod_p2z.f90] / [subroutine: get_estimated_t0]]'
        stop
    end if

    ! allocation array
    allocate( fnl_temp( nlon , nlat  )    )
    allocate( tmp_latconv( nlon , nlat  )    )
    allocate( fnl_lon( nlon ) )
    allocate( fnl_lat( nlat ) )
    allocate( temp( nlon , nlat  )    )
    !copy  coefficient
    fnl_nlon = nlon
    fnl_nlat = nlat

    ! open temperature init
    fid = MISC_get_available_fid()
    irec= 1
    open(fid, file=trim(fnl_tempfiledir)//'/'//trim(fnl_tempfile), &
              form='unformatted', access='direct', recl=4*nlon)
      do ilat = 1, nlat
        read(fid, rec=irec)  temp(:,ilat)
        irec = irec + 1
      end do
    close(fid)

    !--------------------------------------------
    !++ Direction of data insert
    !--------------------------------------------
    ! gen fnl lon-lat
    ! lon
    if ( opt_west2east  ) then
      do ilon = 1, nlon
          fnl_lon(ilon) =  dble(ilon - 1)
      end do
      !
      fnl_temp(:,:) = dble(temp(:,:))
      !
    else if ( .not. opt_west2east  ) then
      do ilon = nlon, 1, -1
          fnl_lon(ilon) =  dble(ilon - 1)
      end do
      !
      do ilon = 1, nlon
        fnl_temp(ilon,:) = dble(temp(nlon-ilon+1,:))
      end do
      !
    end if
    !
    ! lat
    if ( opt_south2north  ) then
      do ilat = 1, ilat
          fnl_lat(ilat) =  -90.d0 + dble(ilat-1)
      end do
      !
      ! pass
      !
    else if ( .not. opt_south2north  ) then
      do ilat = ilat, 1, -1
          fnl_lat(ilat) =  -90.d0 + dble(ilat-1)
      end do
      !
      tmp_latconv(:,:) = fnl_temp(:,:)
      do ilat = 1 , nlat
          fnl_temp(:, ilat) = tmp_latconv(:,nlat-ilat+1)
      end do
      !
    end if

    deallocate(temp)
    deallocate(tmp_latconv)
    return
  end subroutine  get_fnltemp_boundary
  !
  !--------------------------------------------------------------------
  !
  subroutine deallocate_fnltemp()
    implicit none
    ! subroutine for deallocate mod_p2z_fnl array

    deallocate( fnl_temp )
    deallocate( fnl_lon )
    deallocate( fnl_lat )

    return
  end subroutine
  !
  !--------------------------------------------------------------------
  !
  subroutine get_fnllatlon_ID(                  &
                               flon,            &  ! IN
                               ctlat,           &  ! IN
                               lon_id,          &  ! OUT
                               lat_id           &  ! OUT
                             )
    !
    !++ Document
    !   lon - coordination is standerlized to 0 - 360 system
    !
    use mod_hssconst

    implicit none

    ! IN/OUT
    real(kind_real),  intent(in)        :: flon
    real(kind_real),  intent(in)        :: ctlat
    integer,          intent(out)       :: lon_id
    integer,          intent(out)       :: lat_id
    integer                             :: ilon, ilat
    integer                             :: icount
    integer                             :: ierr
    real(kind_real)                     :: ctlon
    real(kind_real)                     :: delta_lon , delta_lat

    !------------------------------------------------------------------

    ! convert 0 to 360
    if (  flon < 0.d0  ) ctlon = flon + 360.d0
    !write(*,*) ctlon , ctlat

    icount  = 0
    do ilon = 1, fnl_nlon-1
        do ilat = 1, fnl_nlat -1
            if (  ctlat >= fnl_lat(ilat) .and. ctlat < fnl_lat(ilat+1)  ) then
                !
                if ( ctlon >= fnl_lon(ilon) .and. ctlon < fnl_lon(ilon+1) ) then
                  lon_id = ilon
                  lat_id = ilat
                  icount = icount + 1
                end if
            endif
        end do
    end do
    !
    ! Return to main loop
    if ( icount >= 1  ) return
    !
    !   boundary I
    !   : lon ===> 359.d0  < current targrt longitude < 360.d0
    !   : lat ===> all
    !
    do ilat = 1, fnl_nlat-1
       !
       if (  ctlat >= fnl_lat(ilat) .and. ctlat < fnl_lat(ilat+1)  ) then
           if ( ctlon >= 359.d0 .and. ctlon <= 360.d0 ) then
               lon_id = fnl_nlon
               lat_id = ilat
               return
            end if
       end if
    end do
    !
    ! Return to main loop
    if ( icount >= 1  ) return
    !
    !   boundary II
    !   : lon ===> 359.d0  < current targrt longitude < 360.d0
    !   : lat ===> +-89.d0 < current targrt latitude  < -+90.d0
    !
    if (  ctlat >= -90.d0  .and. ctlat <= -89.d0   ) then
       !
       if ( ctlon >= 359.d0 .and. ctlon <= 360.d0 ) then
         lon_id = fnl_nlon
         lat_id = 1
         return
       end if
    end if
    !
    if (  ctlat >= 89.d0  .and. ctlat <= 90.d0   ) then
       !
       if ( ctlon >= 359.d0 .and. ctlon <= 360.d0 ) then
         lon_id = fnl_nlon
         lat_id = fnl_nlat
         return
       end if
    end if
    !
    !   boundary III
    !   : lon ===> All
    !   : lat ===> +89.d0 < current targrt latitude  =< +90.d0
    !     inlcude northpole
    !
    if (  ctlat > 89.d0  .and. ctlat <= 90.d0   ) then
       !
       do ilon = 1, fnl_nlon-1
           if ( ctlon >= fnl_lon(ilon) .and. ctlon < fnl_lon(ilon+1) ) then
              lon_id = ilon
              lat_id = ilat
              icount = icount + 1
            end if
       end do
       !
    end if


    if (icount == 0 ) then
        write(*, *) '  !!! stop there is a bug  !!!  '
        write(*, *) '  [ mod | get_fnllatlon_ID ] '
        write(*,*) '  ####    ', ctlon, ctlat
        do ilon = 1, fnl_nlon
        write(90,*) fnl_lon(ilon)
        end do
        do ilat = 1, fnl_nlat
        write(91,*) fnl_lat(ilat)
        end do
        stop
    end if

    return
  end subroutine get_fnllatlon_ID
  !
end module
