!------------------------------------------------------------------------------
!
! ++ Document here
!     Get latlon info. (degree) from ico-grid (radian) 
!
! ++ History 
!
!     Version   Date    Editer   Description
!    ----------------------------------------------------------------------
!     0.0  2018  04  10th   T. Kurihana
!
!------------------------------------------------------------------------------
module  mod_icoid2ll
  implicit none 

  public  :: get_llinfo

contains
  !
  subroutine get_llinfo(                &
                        gall,           & !--IN
                        g_num,          & !--IN
                        rgnid,          & !--IN
                        lon,            & !--OUT
                        lat             & !--OUT
                       )

    !
    ! ++ Load Modules
    !
    use  mod_setupfname,     only : &
         get_ifilename

    use  mod_ico2latlon,     only : &
         MISC_get_latlon


    implicit none

    ! Parameter
    integer,         parameter     :: kind_real = 8
    integer,         parameter     :: kind_char = 256
    integer,         parameter     :: icodim    = 3
    integer,         parameter     :: grd_Xdim  = 1
    integer,         parameter     :: grd_Ydim  = 2
    integer,         parameter     :: grd_zdim  = 3
    !
    ! INOUT
    integer,         intent(in)    :: gall
    integer,         intent(in)    :: g_num
    integer,         intent(in)    :: rgnid
    real(kind_real), intent(out)   :: lon
    real(kind_real), intent(out)   :: lat

    !
    character(kind_char)           :: icodatadir
    character(kind_char)           :: icobasename
    character(128)                 :: gridnumber
    character(40)                  :: prefix

    integer,         save          :: file_id    = 200
    integer,         save          :: input_size = 8

    logical,         save          :: opt_bin = .true.
    logical,         save          :: use_icohgrid = .true.

    real(kind_real)                :: ico_grid(gall, icodim)
    real(kind_real)                :: lon_val
    real(kind_real)                :: lat_val
    real(kind_real)                :: ilon
    real(kind_real)                :: ilat
    real(kind_real)                :: iz
    real(kind_real)                :: pi
    real(kind_real)                :: coef_r2d
    integer                        :: gall_input
    integer                        :: ierr
    integer                        :: irec
    integer                        :: ij
    integer                        :: i, j , k

    namelist  / icoid_param / &
          opt_bin,            &
          use_icohgrid,       &
          input_size,         &
          icodatadir ,        &
          icobasename

    !--------------------------------------------------------------------------

    ! set Pi
    pi = 4.d0 * datan(1.d0)
    coef_r2d = 180.d0 / pi

    ! open namelist
    open(1, file='hssdriver.cnf' )
      read(1, nml=icoid_param, iostat=ierr)
    close(1)
    if (  ierr  /= 0 ) then
      ! OK
      write(*,*) '  ! Stop !   Inappropriate variables in icoid_param   '
      stop
    end if


    ! convert int-id to char-id
    call   get_ifilename(                  &
                          rgnid,           & !--IN
                          gridnumber,      & !--OUT
                          prefix           & !--OUT
                        )

    ico_grid(:,:) = 0.d0
    if (  opt_bin  ) then
      ! binary
      open(file_id, &
         file=trim(icodatadir)//'/'//trim(icobasename)//trim(prefix)//trim(gridnumber), &
         form='unformatted', access='sequential', action='read'       &
          )

          read(file_id) gall_input
          read(file_id) ico_grid(:,grd_Xdim)
          read(file_id) ico_grid(:,grd_Ydim)
          read(file_id) ico_grid(:,grd_Zdim)
     close(file_id)
   !
   else if (  .not. opt_bin  ) then
     ! text
     ! fixme
     stop
      open(file_id, &
         file=trim(icodatadir)//'/'//trim(icobasename)//trim(prefix)//trim(gridnumber), &
         form='formatted',  action='read'       &
          )

          read(file_id, * ) ico_grid(:,:)
     close(file_id)
   !
   else
     !
     !++ Exception Error Processing
     write(*, *) '   Stop : Inappropriate Bool for opt_bin    '
     write(*, *) '   [subroutine | mod_icoid2ll.f90  ]  '
     stop
     !
   end if

   !get info
   if ( use_icohgrid  ) then
     ! not converted to ico ==> ll

     ilon =  ico_grid( g_num , grd_Xdim )
     ilat =  ico_grid( g_num , grd_Ydim )
     iz   =  ico_grid( g_num , grd_Zdim )

     call MISC_get_latlon (           &
                            ilon,     &  !--IN
                            ilat,     &  !--IN
                            iz,       &  !--IN
                            lat_val,  &  !--OUT
                            lon_val   &  !--OUT
                          )

     lon = lon_val * coef_r2d
     lat = lat_val * coef_r2d
    !
   else if ( .not. use_icohgrid ) then
     ! already converted to ico ==> ll
     lon = ico_grid( g_num , grd_Xdim )
     lat = ico_grid( g_num , grd_Ydim )
   else
      !
      write(*,*) '  Error Usage of  no llconverted nicam hgird '
      stop
   end if

  return
  end subroutine

end module 
