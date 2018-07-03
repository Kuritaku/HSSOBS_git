!------------------------------------------------------------------------------
!
! ++ Document 
!     select lat lon info. by delta latlon
!
!     !! Should fix some Develeping points | 0.0
!
! ++ History
!   0.0   2018 03 12``Takuya Kurihana
!
!------------------------------------------------------------------------------
module mod_deltall

  public  :: deltalatlon

  private :: setup_epsval

contains 
  !
  subroutine deltalatlon(   &
              gl, rl,       & ! IN
              gall,         & ! IN
              tgrid,        & ! a temporarly evaluated obs info
              tab_id,       & ! ID of (a) potential rgn grids(grid)
              rgn_tabs,     & ! array of potential rgn grids
              opt_bin,      & ! IN
              dist_mini     & ! OUT
              )
    
    use mod_setupfname, only : &
          get_ifilename

    use mod_qs, only : &
          main_qs

    implicit none 

    integer,         parameter    :: kind_real = 8
    integer,         parameter    :: dist_id   = 9  ! number of return data from quick sort
    integer,         parameter    :: Array_Dim = 4  ! array dimension of sortting data
    !
    !                                Dim 1 : Euclidean distance
    !                                Dim 2 : longitude
    !                                Dim 3 : latitude
    !                                Dim 4 : g grid ID
    !
    integer,         parameter    :: Array_Index = 1  ! sortting target dimension number
    integer,         parameter    :: grd_Xdim  = 1
    integer,         parameter    :: grd_Ydim  = 2
    real(kind_real), parameter    :: undef     = -3215.d0 ! undef value
    real(kind_real), parameter    :: QSundef   = 3215.d0  ! undef value for QS

    integer,         intent(in)   :: gl
    integer,         intent(in)   :: rl
    integer,         intent(in)   :: gall
    integer,         intent(in)   :: tab_id
    integer,         intent(in)   :: rgn_tabs(tab_id)
    real(kind_real), intent(in)   :: tgrid(2)
    logical,         intent(in)   :: opt_bin
    real(kind_real), intent(out)  :: dist_mini(dist_id, Array_Dim, tab_id)

    character(256)                :: basename ! future intent(in)
    character(128)                :: inputdir ! future intent(in)
    character(128)                :: gridnumber
    character(40)                 :: prefix

    integer   :: i, j, k, l
    integer   :: ll
    integer   :: get_num
    integer   :: icount
    integer   :: ignum
    integer   :: idnum
    integer   :: jdnum
    integer   :: itab
    integer   :: log_id
    integer   :: idcount
    integer   :: ArrayID

    real(kind_real)   :: tmp_lon ! deg
    real(kind_real)   :: dlon ! deg
    real(kind_real)   :: dlat ! deg
    real(kind_real)   :: eps_dlon ! deg
    real(kind_real)   :: eps_dlat ! deg
    !real(kind_real)   :: tmp_val(3) future in List

    real(kind_real), allocatable   :: ico_grid(:,:)
    real(kind_real), allocatable   :: tmp_icogrid(:,:)
    real(kind_real), allocatable   :: tmp_near_points(:,:)
    real(kind_real), allocatable   :: near_points(:,:)
    real(kind_real), allocatable   :: tmp_dist_mini(:,:)

    ! fixme Future
    ! Type 
    !type mylist
    !    real(kind_real), pointer  :: array(:,:)
    !end type
    !type(mylist) neighbor_list

    namelist / dll_param  / &
      inputdir , basename

    ! Coefficient
    log_id = 100
    !--------------------------------------------------------------------------

    !open(1, file='dll.cnf')! for test
    open(1, file='hssdriver.cnf') !
      read(1 , nml=dll_param)
    close(1)


    allocate( ico_grid(gall , 2) )

    !------------------------------------------------------
    ! Setup Epsilon latlon
    !------------------------------------------------------
    call setup_epsval(             &
                      gl,          &  ! IN
                      !tlon,       &  ! IN
                      tgrid(2),    &  ! IN
                      eps_dlon,    &  ! OUT
                      eps_dlat     &  ! OUT
                     )

    !------------------------------------------------------
    ! Potential Rgn Outer loop
    !------------------------------------------------------
    do ll = 1, tab_id ! while candidate rgn-ids

      itab = rgn_tabs(ll)
      !write(*,*) itab

      ico_grid(:,:) = 0.d0

      call get_ifilename(   &
                itab,       &
                gridnumber, &
                prefix      &
                )

      ! Open Horizontal grid files
      if( opt_bin ) then
        open(log_id, &
           file=trim(inputdir)//'/'//trim(basename)//'.rgn'//trim(prefix)//trim(gridnumber) , &
           !form='unformatted', access='sequential',  action='read' )
           form='unformatted', access='direct', recl=8*gall*2 ,  action='read' )

            read(log_id, rec=1) ico_grid( : , : )
            !read(log_id) ico_grid( : , grd_Xdim )
            !read(log_id) ico_grid( : , grd_Ydim )

        close(log_id)
        !
      else if ( .not. opt_bin ) then
        !
      ! Open Horizontal grid files
        open(log_id, & 
           file=trim(inputdir)//'/'//trim(basename)//trim(prefix)//trim(gridnumber) , &
           form='formatted', action='read' )

           allocate(tmp_icogrid(gall,2))
           do j = 1, gall
            read(log_id, * ) tmp_icogrid(j, :)

            ico_grid(j,grd_Xdim) = tmp_icogrid(j,grd_Xdim) ! lon
            ico_grid(j,grd_Ydim) = tmp_icogrid(j,grd_Ydim) ! lat
           end do
           deallocate(tmp_icogrid)

        close(log_id)

      else
        write(*,*)  '   warning  '
      end if ! option read binary / text hgrid file


      ! Modify longitude -180 - 180 to 0 360
      do i = 1, gall
          if ( ico_grid(i,1) < 0.d0  ) then
              tmp_lon = ico_grid(i,1)
              ico_grid(i,1) = tmp_lon + 360
          end if
      end do

      ! store delta latlon info into array
      idcount = 1
      allocate( tmp_near_points(gall,Array_Dim) ) ! fixme future Type & pointer
      tmp_near_points(:,:) = undef

      !  Evaluate by loop all g-grid points
      do ignum = 1, gall
          dlon = abs( tgrid(1)  -  ico_grid(ignum,1) )
          dlat = abs( tgrid(2)  -  ico_grid(ignum,2) )
          !write(90,*) dlon, dlat, ico_grid(ignum, 1),ico_grid(ignum, 2)  

        if ( dlon <= eps_dlon .and. dlat <= eps_dlat  ) then
           !tmp_near_points(idcount,1) = dsqrt( dlon*dlon + dlat*dlat ) ! distance on the beta/X-Y plane
           tmp_near_points(idcount,1) = dlon + dlat ! L1 Norm
           tmp_near_points(idcount,2) = ico_grid(ignum,1)
           tmp_near_points(idcount,3) = ico_grid(ignum,2)
           !tmp_near_points(idcount,2) = dlon
           !tmp_near_points(idcount,3) = dlat
           tmp_near_points(idcount,4) = ignum
           write(*,*)' ignum ==     ' , ignum
           idcount = idcount + 1
           !----------------------------------------- 
           ! Now planning below
           !tmp_val(1) = ignum
           !tmp_val(2) = dlon
           !tmp_val(3) = dlat

           ! Linked list will be suitable !
           !Formally stored  Lists 
           !allocate( neighbor_list%array(icount,3) )
           !----------------------------------------- 

        end if
      end do ! ignum

      ! get valid array length
      ArrayID = count( tmp_near_points(:,1) > undef  ) !/  or idcount -1
      write(6,*)  ' number of data dim in ArrayID   ', ArrayID

      ! set array for sort module
      allocate( near_points(ArrayID, Array_Dim) )
      do idnum = 1, Array_Dim
        do jdnum = 1, ArrayID
          near_points(jdnum, idnum) = tmp_near_points(jdnum, idnum) 
        end do
      end do
      deallocate(tmp_near_points)

      ! set for sorting
      allocate ( tmp_dist_mini( ArrayID, Array_Dim ) )
      tmp_dist_mini(:,:) = QSundef ! initialization
      get_num            = ArrayID ! number to get Num. of data
      !if ( ArrayID < dist_id  ) then
      !    get_num            = ArrayID ! number to get Num. of data
      !else if ( ArrayID >= dist_ID  ) then
      !    get_num            = dist_id ! number to get Num. of data
      !else
      !    write(*,*) '  Stop : Error at selecting get-num number '
      !    stop
      !end if


      call main_qs(                 &
                    ArrayID,        & ! IN
                    Array_Dim,      & ! IN
                    Array_Index,    & ! IN
                    get_num,         & ! IN
                    !dist_id,        & ! IN
                    near_points,    & ! IN
                    tmp_dist_mini   & ! OUT
                   )

      ! check sorted array
      write(6,*) '  '
      do i = 1, get_num
         write(6,"(4f12.5)") ( tmp_dist_mini(i,j) , j=1,Array_Dim)
      end do
      !write(6,*) '  '

      icount = 1
      if ( ArrayID < dist_id  ) then
        do jdnum = 1, get_num
            dist_mini(  icount ,: , ll ) = tmp_dist_mini( jdnum , : )
            icount = icount + 1
        end do
        !
        ! filled by undef in unnecessaty parts in dist-Array
        dist_mini( icount:dist_id , :  , ll ) = undef
        !
      else if ( ArrayID >=  dist_id  ) then
        do jdnum = 1, get_num
            dist_mini(  icount ,: , ll ) = tmp_dist_mini( jdnum , : )
            icount = icount + 1
        end do
        !
       end if 

      deallocate(  tmp_dist_mini  )

    end do ! itab

    return
  end subroutine deltalatlon

  !----------------------------------------------------------------------------

  subroutine setup_epsval(            &
                          glevel,     &  ! IN
                          !tlon,       &  ! IN
                          tlat,       &  ! IN
                          eps_lon,    &  ! OUT
                          eps_lat     &  ! OUT
                         )

    !
    ! ++ Document 
    !     Define appropriate epsilon value for delta-latlon
    !
    !     Resolution and Glevel (globally homogenous grid case)
    !         glevel=11, dx~3.5km
    !         glevel=10, dx~7km
    !         glevel= 9, dx~14km
    !         glevel= 8, dx~28km
    !         glevel= 7, dx~56km
    !         glevel= 6, dx~112km
    !         glevel= 5, dx~224km
    !     # 'dx' doesn't depend on rlevel.
    !
    ! ++ history
    !
    !     ver   date            editer
    ! ------------------------------------------------------------
    !     0.0  2018 3 19th    Takuya Kurihana
    !          Assume 1deg of lat-lon is about 100km
    !          * future ====> automatic parameter estimation
    !          ===> bug included
    !
    !     1.0  2018 6 6th    Takuya Kurihana
    !          Bug fix
    !          add auto-eps_lon get function for safe system run

    implicit none

    integer,         parameter    :: kind_real = 8
    integer,        intent(in)    :: glevel
    !real(kind_real),intent(in)    :: tlon     ! target longitude
    real(kind_real),intent(in)    :: tlat     ! target latitude
    real(kind_real),intent(out)   :: eps_lon
    real(kind_real),intent(out)   :: eps_lat

    !--------------------------------------------------------------------------

    eps_lon = get_eps_lon( tlat,glevel )
    if       ( glevel  == 5 ) then
      !eps_lat = 3.d0
      eps_lat = 3.d0
    else  if ( glevel  == 6 ) then
      !eps_lat = 1.5d0
      eps_lat = 2.d0
    else  if ( glevel  == 7 ) then
      !eps_lat = 0.75d0
      eps_lat = 1.d0
    else  if ( glevel  == 8 ) then
      !eps_lat = 0.275d0
      eps_lat = 0.5d0
    else  if ( glevel  == 9 ) then
      !eps_lat = 0.135d0
      eps_lat = 0.25d0
    else  if ( glevel  >= 10 ) then
      !eps_lat = 0.0625d0
      eps_lat = 0.125d0
    else
        ! Exception 
        !
        write(6,*)    '   !  Exception of eps-latlon parameter  !  '
        write(6,*)    '  Error at [ subroutine ] setup_epsval / [module] mod_dlatlon.f90  '
        !
        stop
    end if 

    return
  end subroutine  setup_epsval
  !
  !----------------------------------------------------------------------------
  !
  function  get_eps_lon( lat , gl ) & 
            result( eps_lon )
    !
    !++ module 
    use mod_hssconst,    only : &
        kind_real
    !
    implicit none

    real(kind_real), parameter  :: r        = 6356.d0
    real(kind_real), parameter  :: pi       = 3.14159265358979323846
    real(kind_real), parameter  :: dist_gl5 = 224.d0

    integer                     :: gl       ! IN
    real(kind_real)             :: lat      ! IN
    real(kind_real)             :: eps_lon  ! OUT [km]
    real(kind_real)             :: unit_lon      ![km]
    real(kind_real)             :: gl_dist       ![km]
    real(kind_real)             :: coef_r2d

    !--------------------------------------------------------------------------

    ! delta longitude of 1 deg.
    coef_r2d = 2.d0*pi/360.d0
    unit_lon = coef_r2d*r*cos(lat*coef_r2d)

    ! threshold distance for glevel
    gl_dist = dist_gl5/(gl-5+1)

    ! return epsilon longitude
    eps_lon  =  gl_dist / unit_lon

  end function
  !
end module 
