!------------------------------------------------------------------------------
!
!++  Document 
!    Main Program of High Speed Super-OBservation System
!
!++  Hisotry
!
!       Version    Date       Editer        Description
!    ----------------------------------------------------------------
!        0.0      2018-03-21  T. Kurihana   Building Driver HSSOBS
!        1.0      2018-04-09  T. Kurihana   Coding Alpha ver.
!        1.1      2018-04-25  T. Kurihana   Bug-Fix in dacos(X) when X=+/-1
!        2.0      2018-05-12  T. Kurihana   Code Develop in Rgrid Table
!        2.1      2018-05-12  T. Kurihana   Separate IO Input
!        2.2      2018-05-12  T. Kurihana   Make Module for parameters
!        2.3      2018-05-16  T. Kurihana   Add module: get ico-zlayer number
!        3.0      2018-06-09  T. Kurihana   Modify Program for NCEP Prepbufr
!        3.1      2018-06-10  T. Kurihana   Fix/Add input/output lon data
!
!
!++  Input/Output Data Structure
!       Longitude value should be 0 <= lon  < 360
!       This program deals with data as 0to360 longitude system
!       ===> Inputdata -180 to 180 is modified in program
!
!       * Output option can be selective either 0to360 or -180to180
!
!
!++  API
!
!------------------------------------------------------------------------------
program driver_hssobs
  !
  ! ++ Modules to build the entire system
  !

  use  mod_hssconst

  use  mod_common_prepbufr

  use  mod_ico2latlon ,       only :  &
       setup_ico2ll

  use  mod_setup_zlayer,      only :  &
       setup_icozgrid,                &
       setup_icozlayer

  use  mod_conv_vcoords_fnl,  only :  &
       fnl_temp,                      &
       get_fnllatlon_ID,              &
       get_fnltemp_boundary,          &
       deallocate_fnltemp,            &
       get_pval

  use  mod_compave,           only :  &
       !setup_ave,                 &
       !get_aveval
       setup_obs_ave,                 &
       get_obs_aveval

  use  mod_icoid2ll,          only :  &
       get_llinfo
  !
  ! ++ Flowchart of Superobservation
  !

  !----------------------------------------------------------------------------

  implicit none

  ! before ico2ll
  character(kind_char)         :: obsdatadir
  character(kind_char)         :: obsinputfile
  character(kind_char)         :: inputdir
  character(kind_char)         :: outputdir
  character(kind_char)         :: outputname
  !after ico2ll
  character(kind_char)         :: rgnfilename  ='bin_hgrid'
  character(kind_char)         :: inputfiledir ='./'
  character(kind_char)         :: zgridinfo ='./ico_38zlayer.txt'
  character(kind_char)         :: logfname  ='prepbufr_'

  ! save variables
  logical,         save        :: offline    = .true.
  logical,         save        :: opt_bin    = .true.
  logical,         save        :: opt_pass   = .true.
  logical,         save        :: BuildRgnTableIs
  logical,         save        :: DoLoopIs
  !
  real(kind_real), parameter   :: RNundef    = -9999.d0
  real(kind_real), parameter   :: LUTundef   = -9999.d0
  !
  integer,         save        :: gl
  integer,         save        :: rl
  integer,         save        :: gall
  integer,         save        :: lall
  integer,         save        :: zall
  integer,         save        :: irec
  integer,         save        :: input_size
  integer,         save        :: tmp_rgnnumber

  real(kind_real)              :: garray(lut_dim, lut_elem)
  real(kind_real)              :: dist_mini(lut_elem)
  real(kind_real)              :: ctlon
  real(kind_real)              :: ctlat
  real(kind_real)              :: zval
  real(kind_real)              :: physval
  real(kind_real)              :: slon
  real(kind_real)              :: slat
  real(kind_real)              :: err
  real(kind_real)              :: pbtype
  real(kind_real)              :: plevel
  real(kind_real)              :: T0
  real(kind_real)              :: hss_dist
  real(kind_real), allocatable :: lon(:)
  real(kind_real), allocatable :: lat(:)
  real(kind_real), allocatable :: ico_zgrid(:)
  !real(kind_real), allocatable , save :: ave_box(:,:,:)   ! atm(grid,zlayer, rgn)
  !real(kind_real), allocatable , save :: ave_count(:,:,:) ! atm(grid,zlayer, rgn)
  !
  ! Modify average box for prepbufr
  real(kind_real), allocatable , save :: ave_box(:,:,:,:)   ! atm(grid,zlayer, rgn, obstype)
  real(kind_real), allocatable , save :: ave_count(:,:,:,:) ! atm(grid,zlayer, rgn, obstype)
  real(kind_real), allocatable , save :: ave_err_box(:,:,:,:)  !err(grid,zlayer, rgn, obstype)
  real(kind_real), allocatable , save :: ave_err_count(:,:,:,:)!err(grid,zlayer, rgn, obstype)

  integer                      :: g_num
  integer                      :: rgn_tabs(tab_id)
  integer                      :: icount
  integer                      :: rgnid
  integer                      :: ij
  integer                      :: irgn
  integer                      :: l
  integer                      :: ielem
  integer                      :: ierr
  integer                      :: obsall
  integer                      :: ilon, ilat
  integer                      :: pb_id
  integer                      :: vid
  integer                      :: hss_rgn
  integer                      :: hss_g
  integer                      :: hss_z
  integer                      :: lon_id
  integer                      :: lat_id

  logical                      :: opt_outlon_0to360 = .false.
  logical                      :: opt_zcoords       = .false.
  logical                      :: opt_pcoords       = .true.

  ! iteration variables
  integer                      :: i,j,k ! debug
  integer                      :: iobs  ! debug
  integer                      :: ii

  ! time cal. vars
  real(kind_real)              :: t1, t2

  ! namelist
  namelist  / hssdriver_param / &
              gl,               &
              rl,               &
              zall,             &
              obsall,           &
              offline,          &
              opt_bin,          &
              opt_outlon_0to360,&
              opt_zcoords,      &
              opt_pcoords,      &
              obsdatadir,       &
              obsinputfile,     &
              inputdir,         &
              input_size,       &
              outputdir,        &
              outputname,       &
              rgnfilename,      &
              inputfiledir,     &
              zgridinfo,        &
              logfname

  namelist  / prepbufr_param  / &
              pb_id

  !----------------------------------------------------------------------------

  !--------------------------------------------------------
  !++ Open Namelist
  !--------------------------------------------------------
  open(1, file='hssdriver.cnf' )
    read(1, nml=hssdriver_param, iostat=ierr )
    read(1, nml=prepbufr_param,  iostat=ierr )
  close(1)
  ! fixme future here module
  if (  ierr == 0 ) then
     write(*,*) '      Found Namelist    '
  else
     stop
  end if

  !--------------------------------------------------------
  !++ Setup Coefficients
  !--------------------------------------------------------
  ! gall lall
  gall = (  2**(gl - rl )+2 )**2
  lall = (  2**rl )**2 * 10


  ! setup store average matrix
  !allocate( ave_box(   gall , zall , 0:lall-1 ) )
  !allocate( ave_count( gall , zall , 0:lall-1 ) )
  allocate( ave_box(   gall , zall , 0:lall-1, nobtype ) )
  allocate( ave_count( gall , zall , 0:lall-1, nobtype ) )
  allocate( ave_err_box(   gall , zall , 0:lall-1, nobtype ) )
  allocate( ave_err_count( gall , zall , 0:lall-1, nobtype ) )

  ! initialize
  ave_box(  :,:,:,:) = 0.d0
  ave_count(:,:,:,:) = 0.d0
  ave_err_box(  :,:,:,:) = 0.d0
  ave_err_count(:,:,:,:) = 0.d0


  !--------------------------------------------------------
  !++ Setup ico 2 ll
  !--------------------------------------------------------
  allocate( lon(gall) , lat(gall) )
  do rgnid = 0, lall-1
    call setup_ico2ll(            &
                 gl  , rl ,       & ! IN
                 gall, lall ,     & ! IN
                 rgnid,           & ! IN
                 inputdir,        & ! IN
                 outputdir,       & ! IN
                 outputname,      & ! IN
                 offline,         & ! IN
                 opt_bin,         & ! IN
                 lat, lon         & ! OUT
                     )

  end do
  write(*,*) '      ####  Finish  :   Ico2ll    #####    '

  !--------------------------------------------------------
  !++ Setup ico zgrid
  !--------------------------------------------------------
  allocate( ico_zgrid(zall) )
  call setup_icozgrid(              &
                       zall,        &  !--IN
                       zgridinfo,   &  !--IN
                       ico_zgrid    &  !--OUT
                     )

  !--------------------------------------------------------
  !++ Setup temperature boundary data array
  !--------------------------------------------------------
  call get_fnltemp_boundary()


  call cpu_time(t1)

  !==================================================================
  !
  !                        Start Main Process
  !
  !==================================================================

  write(*,*) "================================================================"
  write(*,*) " "
  write(*,*) " ########              Start Main Process                #######"
  write(*,*) " "
  write(*,*) "================================================================"

  ! Intialization
  irec = 1
  BuildRgnTableIs   = .true.
  DoLoopIs          = .true.
  !
  ! Main Loop
  do iobs = 1, obsall

  !--------------------------------------------------------
  !++ Open Obs data
  !--------------------------------------------------------
    call open_prepbufr_data(                     &
                            input_size,          &  !--IN
                            obsdatadir,          &  !--IN
                            obsinputfile,        &  !--IN
                            ctlon,               &  !--OUT
                            ctlat,               &  !--OUT
                            zval,                &  !--OUT
                            physval,             &  !--OUT
                            err,                 &  !--OUT
                            pbtype,              &  !--OUT
                            irec                 &  !--INOUT
                           )

  !--------------------------------------------------------
  !++ Compute R-regions & G-grids
  !--------------------------------------------------------
    DoLoopIs =  .true.
    do while ( DoLoopIs )
      call get_rgvalue(                   &
                        gl,               &
                        rl,               &
                        gall,             &  !--IN
                        lall,             &  !--IN
                        ctlon,            &  !--IN
                        ctlat,            &  !--IN
                        offline,          &  !--IN
                        opt_bin,          &  !--IN
                        rgnfilename,      &  !--IN
                        inputfiledir,     &  !--IN
                        hss_g,            &  !--OUT
                        hss_rgn,          &  !--OUT
                        hss_dist,         &  !--OUT
                        DoLoopIs,         &  !--OUT
                        tmp_rgnnumber,    &  !--INOUT
                        BuildRgnTableIs   &  !--INOUT
                      )

    end do ! while

    ! compute z-val [m]
    call setup_icozlayer(             &
                          zall,       &  !--IN
                          ico_zgrid,  &  !--IN
                          zval,       &  !--IN
                          hss_z       &  !--OUT
                        )

    ! store data for  average computation
    !call setup_ave(                 &
    !                  hss_g,        &
    !                  hss_rgn,      &
    !                  hss_z,        &
    !                  gall,         &
    !                  lall,         &
    !                  zall,         &
    !                  physval,      &
    !                  ave_box,      &  !--OUT
    !                  ave_count     &  !--OUT
    !              )

    ! Modify below for prepbufr data
    ! store prepbufr data for average computation
    call setup_obs_ave(               &
                        hss_g,        &
                        hss_rgn,      &
                        hss_z,        &
                        gall,         &
                        lall,         &
                        zall,         &
                        physval,      &
                        err,          &
                        pbtype,       &
                        ave_box,      &  !--OUT
                        ave_count,    &  !--OUT
                        ave_err_box,  &  !--OUT
                        ave_err_count &  !--OUT
                      )


  end do ! outer obs number loop
  close(1000) ! file end
  call cpu_time(t2)
  ! ======> until here store data
  !----------------------------------------------------------------------------

  !call get_aveval(          &
  !              gall,       &
  !              lall,       &
  !              zall,       &
  !              ave_box,    &
  !              ave_count   &
  !               )

  call get_obs_aveval(                 &  !<--IN
                       gall,           &
                       lall,           &
                       zall,           &
                       ave_box,        &  !<--INOUT
                       ave_count,      &
                       ave_err_box,    &
                       ave_err_count   &
                     )


  !--------------------------------------------------------
  !++ Open Log File
  !--------------------------------------------------------
  vid  = get_prepbufr_varID( pb_id )
  open(300, file='log_'//trim(logfname)//'.txt', form='formatted', action='write'   )
  open(400, file='bin_'//trim(logfname)//'.dat', form='unformatted', action='write' )
  write(6,"(4A5,2A12)") 'G_num'  , 'Z' , 'Rgn' , 'OBS','Ave. Val', 'Ave. Err'
  do l = 1, nobtype
    do k = 0, lall-1
      do j = 1, zall
        do i = 1, gall
          if (  ave_count( i, j, k, l )  > 0.d0 ) then
            write(*,"(4I5,2f12.6)") i , j , k , l, &
                                    ave_box( i , j ,k, l), &
                                    ave_err_box( i , j ,k, l)
            !write(700,*) i , j , k
            call  get_llinfo(           &
                              gall,     & !--IN
                              i,        & !--IN g_num
                              k,        & !--IN rgnid
                              slon,     & !--OUT
                              slat      & !--OUT
                          )
            ! save log
            write(300,"(1i6, 6f14.7)") vid, slon, slat , ico_zgrid(j),    &
                                       ave_box( i , j , k , l),      &
                                       ave_err_box( i , j , k , l),  &
                                       sngl(l)


            ! output longitude option
            if ( .not. opt_outlon_0to360  )  then  ! -180 to 180 prepbufr style
                if ( slon > 180.d0 ) slon = slon - 360.d0
            else if ( opt_outlon_0to360  )  then  ! 0 to 360
                if ( slon < 0.d0 )   slon = slon + 360.d0
            end if

            ! save binary
            if ( opt_zcoords  )  then
              write(400) sngl(vid), sngl(slon), sngl(slat),       &
                         sngl(ico_zgrid(j)) ,                     &
                         sngl(ave_box( i , j , k , l)),           &
                         sngl(ave_err_box( i , j , k , l)),       &
                         sngl(l)

            else if (  opt_pcoords ) then

              ! get nearest T value
              call  get_fnllatlon_ID(         &
                                     slon,    & !--IN
                                     slat,    & !--IN
                                     lon_id,  & !--OUT
                                     lat_id   & !--OUT
                                    )

              T0= fnl_temp(lon_id, lat_id)
              plevel  =  get_pval( ico_zgrid(j), T0 )


              write(400) sngl(vid), sngl(slon), sngl(slat) ,      &
                         sngl(plevel),                            &
                         sngl(ave_box( i , j , k , l)),           &
                         sngl(ave_err_box( i , j , k , l)),       &
                         sngl(l)
            end if

          end if
        end do
      end do
    end do
  end do

  !call cpu_time(t2)
  write(6,*)  '####   CPU Time in searching    ' , t2-t1 , 'sec     ####'
  !write(6,*)  'estimation 1 million data time ', 1.0d+6*(t2-t1)/60.d0 , 'min'
  !write(6,*)  'estimation 1 million data time ', 1.0d+2*(t2-t1)/60.d0 , 'min'
  close(300)
  close(400)

  ! deallocate
  call deallocate_fnltemp()
  deallocate( ave_box, ave_err_box, ave_count, ave_err_count   )

end program driver_hssobs

!------------------------------------------------------------------------------

subroutine  get_rgvalue(                  &
                        gl,               &
                        rl,               &
                        gall,             &  !--IN
                        lall,             &  !--IN
                        ctlon,            &  !--IN
                        ctlat,            &  !--IN
                        offline,          &  !--IN
                        opt_bin,          &  !--IN
                        rgnfilename,      &  !--IN
                        inputfiledir,     &  !--IN
                        hss_g,            &  !--OUT <-- int
                        hss_rgn,          &  !--OUT <-- int
                        hss_dist,         &  !--OUT <-- long
                        DoLoopIs,         &  !--OUT
                        tmp_rgnnumber,    &  !--INOUT
                        BuildRgnTableIs   &  !--INOUT
                       )

    !
    !++ module
    use   mod_hssconst

    use  mod_hssll,         only :  &
         setup_latlontab

    use  mod_hssrgn,        only :  &
         RgnLoop

    implicit none
    !
    ! Parameter
    !integer,         parameter   :: kind_real = 8
    !integer,         parameter   :: kind_char = 256
    !integer,         parameter   :: RgnCandMaxID = 9
    !integer,         parameter   :: rgmin_elem = 5
    !integer,         parameter   :: td_dim     = 2
    !integer,         parameter   :: grd_Xdim   = 1
    !integer,         parameter   :: grd_Ydim   = 2
    !
    ! INPUT Variables
    integer,         intent(in)     :: gall
    integer,         intent(in)     :: lall
    real(kind_real), intent(in)     :: ctlon
    real(kind_real), intent(in)     :: ctlat
    !
    !logical ,         intent(in)    :: opt_pass
    logical ,         intent(in)    :: opt_bin
    logical ,         intent(in)    :: offline
    !
    character(kind_char), intent(in)   :: rgnfilename
    character(kind_char), intent(in)   :: inputfiledir
    !
    ! OUTPUT Variables
    !real(kind_real), intent(out)    :: hss_g
    !real(kind_real), intent(out)    :: hss_rgn
    integer,         intent(out)    :: hss_g
    integer,         intent(out)    :: hss_rgn
    real(kind_real), intent(out)    :: hss_dist
    logical,         intent(out)    :: DoLoopIs
    integer,         intent(inout)  :: tmp_rgnnumber
    logical,         intent(inout)  :: BuildRgnTableIs


    ! Local Variables
    ! Int
    integer                      :: gl
    integer                      :: rl
    integer                      :: rgnid_dim
    integer                      :: ierr
    !
    ! Int -Array
    integer                      :: RgnCandArray(RgnCandMaxID)
    integer,         allocatable :: InRgnCandArray(:)
    !
    ! Real
    real(kind_real)              :: rgngrid_mini(rgmin_elem)
    !                               target-G&R-grid info
    !                               rgnmin_elem ==> 5
    !                                 1 : g-grid id
    !                                 2 : longitude
    !                                 3 : latitude
    !                                 4 : distance
    !                                 5 : rgn id
    !
    real(kind_real)              :: tmp_grid(td_dim)
    !
    logical                      :: opt_pass


    !----------------------------------------------------------------------------

    !++ Initialize
    opt_pass            = .true.
    tmp_grid(grd_Xdim)  = ctlon
    tmp_grid(grd_Ydim)  = ctlat

    if  (   BuildRgnTableIs  ) then

        ! module ll
        call    setup_latlontab(              &
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

        write(*, *) '  !!!!!!! check   rgnid_dim ', rgnid_dim
        allocate(  InRgnCandArray( rgnid_dim  ) )

        ! move value
        InRgnCandArray(1:rgnid_dim ) =  RgnCandArray(1:rgnid_dim)

    else  if (    .not. BuildRgnTableIs   ) then
        rgnid_dim       = 1
        allocate(InRgnCandArray(rgnid_dim))
        InRgnCandArray(rgnid_dim) = tmp_rgnnumber
    else !
      !
      write(*, *) '  Stop Error Processing '
      stop
      !
    end if ! table loop


    ! module rgn
    call  RgnLoop(                  &
                  gl,               & ! IN
                  rl,               & ! IN
                  gall,             & ! IN
                  tmp_grid,         & ! IN
                  rgnid_dim,        & ! IN
                  offline,          & ! IN
                  opt_bin,          & ! IN
                  rgnfilename,      & ! IN
                  inputfiledir,     & ! IN
                  InRgnCandArray,   & ! IN
                  rgngrid_mini      & ! OUT
                  )

    ! judge loop table or use table again
    if (      rgngrid_mini(1) < 0  ) then
        !
        !write(*,*) '  Re-Construct Look-Up-Table    '
        !write(*,*) '                                '
        BuildRgnTableIs = .true.
        deallocate( InRgnCandArray  )
        !
        DoLoopIs = .true.
        tmp_rgnnumber = -1
        return
        !
        !++ Just In Case Safety Stopper
        !   ==> if 'return ' did not work, forcely terminate entire program
        !
        write(*,*)  '   !!      Termination Entire Program      !!   '
        stop

    else if (   rgngrid_mini(1) >= 0 ) then
        !
        BuildRgnTableIs = .false.
        DoLoopIs        = .false.
        tmp_rgnnumber   = int(rgngrid_mini(5))
        !rgnid_dim       = 1

        deallocate( InRgnCandArray  )
        !allocate(InRgnCandArray(rgnid_dim))
        !InRgnCandArray(rgnid_dim) = tmp_rgnnumber
    end if

    ! Insert Return Data
    hss_g   = int(rgngrid_mini(1)) ! grid
    hss_dist= rgngrid_mini(4)      ! distance
    hss_rgn = int(rgngrid_mini(5)) ! rgn


  return
end subroutine

!------------------------------------------------------------------------------

subroutine  open_prepbufr_data(                     &
                               input_size,          &  !--IN
                               obsdatadir,          &  !--IN
                               obsinputfile,        &  !--IN
                               ctlon,               &  !--OUT
                               ctlat,               &  !--OUT
                               zlev,                &  !--OUT
                               physval,             &  !--OUT
                               err,                 &  !--OUT
                               pbtype,              &  !--OUT
                               irec                 &  !--INOUT
                              )


  !
  !++ History Log
  !     log-No.    date         comment
  ! -------------------------------------------------------------------
  !       1    2018-06-05   change def(open_data) for NCEP prepbufr
  !       2    2018-06-09   Split as new subroutine: add matrix dim for output
  !       3    2018-06-10   input longitude is standerlized as 0-360 longitude
  !                         system
  !
  !++ module
  use   mod_hssconst

  !use   mod_get_avefid,   only  : &
  !      MISC_get_available_fid

  use   mod_conv_vcoords_fnl,  only  : &
        fnl_temp,                      &
        get_fnllatlon_ID,              &
        get_zval


  implicit none

  ! fixme temporarly-parameter
  integer,  parameter                   :: obs_elem   = 6 ! prepbufr data
  integer,  parameter                   :: grd_OBSdim = 4 ! Dim. of OBS
  integer,  parameter                   :: grd_ERRdim = 5 ! Dim. of Errer
  integer,  parameter                   :: grd_TYPdim = 6 ! Dim. of Type
  !
  integer,              intent(in)      :: input_size
  character(kind_char), intent(in)      :: obsdatadir
  character(kind_char), intent(in)      :: obsinputfile
  real(kind_real),      intent(out)     :: ctlon
  real(kind_real),      intent(out)     :: ctlat
  real(kind_real),      intent(out)     :: physval
  real(kind_real),      intent(out)     :: zlev
  real(kind_real),      intent(out)     :: err
  real(kind_real),      intent(out)     :: pbtype
  integer,              intent(inout)   :: irec
  !
  integer                               :: fid
  integer                               :: lon_id
  integer                               :: lat_id
  real(kind_sngl)                       :: sngl_ico_grid(obs_elem)
  real(kind_real)                       :: plev
  real(kind_real)                       :: T0

  !----------------------------------------------------------------------------

  !fid = MISC_get_available_fid()
  fid = 1000
  if ( input_size  == 4 ) then
    open(fid, file=trim(obsdatadir)//'/'//trim(obsinputfile),   &
              form='unformatted',  access='direct',             &
              recl=kind_sngl*obs_elem, action='read'            &
        )
      read(fid, rec=irec)  sngl_ico_grid(:)
      irec = irec + 1

      ctlon               = dble( sngl_ico_grid(grd_Xdim) )
      if ( ctlon < 0.d0  ) ctlon = ctlon+ 360.d0
      ctlat               = dble( sngl_ico_grid(grd_Ydim) )
      plev                = dble( sngl_ico_grid(grd_Pdim) )
      !write(*,*) ctlon , ctlat , plev

      ! get nearest T value
      call  get_fnllatlon_ID(           &
                              ctlon,    & !--IN
                              ctlat,    & !--IN
                              lon_ID,   & !--OUT
                              lat_ID    & !--OUT
                            )

      T0                  = fnl_temp(lon_id, lat_id)
      zlev                = get_zval(plev, T0)
      physval             = dble( sngl_ico_grid(grd_OBSdim) )
      err                 = dble( sngl_ico_grid(grd_ERRdim) )
      pbtype              = dble( sngl_ico_grid(grd_TYPdim) )
     !
   else 
     !
     write(*,*) '    !!!!      Stop Main Program     !!!!     '
     write(*,*) '    Error : Inappropriate input_size number other than 4 or 8   '
     !
     stop
     !
   end if 
  !close(fid)
  write(*,*) ctlon, ctlat, plev, physval, err, pbtype
  !write(*,*) ' hello stop'
  !if (irec == 3 ) stop
  return

end subroutine

!------------------------------------------------------------------------------

subroutine  open_data(                     &
                      input_size,          &  !--IN
                      obsdatadir,          &  !--IN
                      obsinputfile,        &  !--IN
                      sngl_ico_grid,       &  !--OUT
                      dble_ico_grid,       &  !--OUT
                      ctlon,               &  !--OUT
                      ctlat,               &  !--OUT
                      zlev,                &  !--OUT
                      physval,             &  !--OUT
                      irec                 &  !--INOUT
                     )


  !
  !++ History Log
  !     log-No.    date         comment
  ! -------------------------------------------------------------------
  !       1    2018-06-05   change def(open_data) for NCEP prepbufr
  !
  !++ module
  use   mod_hssconst

  !use   mod_get_avefid,   only  : &
  !      MISC_get_available_fid

  use   mod_conv_vcoords_fnl,  only  : &
        fnl_temp,                      &
        get_fnllatlon_ID,              &
        get_zval


  implicit none

  ! fixme temporarly-parameter
  integer,  parameter                   :: obs_elem = 4 ! prepbufr data
  !
  integer,              intent(in)      :: input_size
  character(kind_char), intent(in)      :: obsdatadir
  character(kind_char), intent(in)      :: obsinputfile
  real(kind_sngl),      intent(out)     :: sngl_ico_grid(obs_elem)
  real(kind_real),      intent(out)     :: dble_ico_grid(obs_elem)
  real(kind_real),      intent(out)     :: ctlon
  real(kind_real),      intent(out)     :: ctlat
  real(kind_real),      intent(out)     :: physval
  real(kind_real),      intent(out)     :: zlev
  integer,              intent(inout)   :: irec
  integer                               :: lon_id
  integer                               :: lat_id
  real(kind_real)                       :: plev
  real(kind_real)                       :: T0

  !----------------------------------------------------------------------------

  if ( input_size  == 4 ) then
    open(1000, file=trim(obsdatadir)//'/'//trim(obsinputfile),   &
               form='unformatted',  access='direct',             &
               recl=kind_sngl*obs_elem, action='read'            &
        )
      read(1000, rec=irec)  sngl_ico_grid(:)
      irec = irec + 1

      !fixme  here
      ctlon               = dble( sngl_ico_grid(grd_Xdim) )
      if ( ctlon < 0.d0  )  ctlon = ctlon + 360.d0
      ctlat               = dble( sngl_ico_grid(grd_Ydim) )
      plev                = dble( sngl_ico_grid(grd_Pdim) )
      !write(*,*) ctlon , ctlat , plev

      ! get nearest T value
      call  get_fnllatlon_ID(           &
                              ctlon,    & !--IN
                              ctlat,    & !--IN
                              lon_ID,   & !--OUT
                              lat_ID    & !--OUT
                            )

      T0                  = fnl_temp(lon_id, lat_id)
      zlev                = get_zval(plev, T0)
      physval             = dble( sngl_ico_grid(obs_elem) )

  else if (  input_size == 8) then
    open(1000, file=trim(obsdatadir)//'/'//trim(obsinputfile) , &
               form='unformatted', access='direct',             &
               recl=kind_real*obs_elem,    action='read'        &
        )
      ! lon
      read(1000, rec=irec)  dble_ico_grid(:)
      irec = irec + 1

      !fixme  here
      ctlon               =  dble_ico_grid(grd_Xdim)
      ctlat               =  dble_ico_grid(grd_Ydim)
      plev                =  dble_ico_grid(grd_Zdim)

      ! get nearest T value
      call  get_fnllatlon_ID(           &
                              ctlon,    & !--IN
                              ctlat,    & !--IN
                              lon_ID,   & !--OUT
                              lat_ID    & !--OUT
                            )

      T0                  =  fnl_temp(lon_id, lat_id)
      zlev                =  get_zval(plev, T0)
      physval             =  dble_ico_grid(4)

   else 
     !
     write(*,*) '    !!!!      Stop Main Program     !!!!     '
     write(*,*) '    Error : Inappropriate input_size number other than 4 or 8   '
     !
     stop
     !
   end if 
  !write(*,*) ctlon, ctlat, plev, T0, zlev, physval, irec
  !write(*,*) ' hello stop'
  !stop
  return
end subroutine  open_data
