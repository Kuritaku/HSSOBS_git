!------------------------------------------------------------------------------
!
! ++ Document
!     search nearest g-grid from proir selected LookUpTable (mod_dlatdlon.f90)
!
! ++ History
!     0.0   2018 3 19th  Takuya Kurihana  :  code alpha-version
!     0.1   2018 3 21th  Takuya Kurihana  :  All test case is clear
!
!------------------------------------------------------------------------------
module  mod_gridsearch
  implicit none

  public   :: nearest_Ggrid , & 
              setup_lutelem , & 
              setup_gridtable

  !private  :: setup_gridtable

contains
  !
  subroutine  nearest_Ggrid(         &
                      !gl ,           & ! IN
                      !rl ,          & ! IN
                      lut_dim,       & ! IN
                      lut_elem,      & ! IN
                      sort_id,       & ! IN
                      grid_lutab,    & ! IN
                      dist_mini      & ! OUT
                           )
    !
    use mod_qs , only : &
        main_qs

    implicit none

    integer,        parameter     :: kind_real  = 8
    integer,        parameter     :: get_num    = 1
    real(kind_real),parameter     :: undef      = 9999.d0 ! for DESO sort
    !integer,        intent(in)    :: gl
    !integer,        intent(in)    :: rl
    integer,        intent(in)    :: lut_dim    ! dimension of Look up table(LUT)
    integer,        intent(in)    :: lut_elem   ! number of elements in LUT
    integer,        intent(in)    :: sort_id    ! ID to sort
    real(kind_real),intent(in)    :: grid_lutab( lut_dim, lut_elem )
    real(kind_real),intent(out)   :: dist_mini(lut_elem)

    integer                       :: ildim , ielem
    real(kind_real)               :: tmp
    real(kind_real)               :: self_grid_lutab( lut_dim , lut_elem )

    !--------------------------------------------------------------------------

    ! convert  minus distance undef value to plus
    self_grid_lutab( : , : ) = grid_lutab( : , : )
    ielem = 4 ! distance
      do ildim = 1, lut_dim
        tmp = grid_lutab( ildim , ielem )
        if ( tmp <  0.d0   ) then
            self_grid_lutab( ildim , ielem ) = undef
        end if
      end do

    ! convert all minus undef value to plus
    !self_grid_lutab( : , : ) = grid_lutab( : , : )
    !do ielem = 1, lut_elem
    !  do ildim = 1, lut_dim
    !    tmp = grid_lutab( ildim , ielem )
    !    if ( tmp <  0.d0   ) then
    !        self_grid_lutab( ildim , ielem ) = undef
    !    end if
    !  end do
    !end do


    call  main_qs(            &
            lut_dim  ,        & ! IN
            lut_elem ,        & ! IN
            sort_id  ,        & ! IN
            get_num  ,        & ! IN
            self_grid_lutab , & ! IN
            dist_mini         & ! OUT
                 )

    !
    ! ++ fixme future add flag dist_mini is within acceptable distance
    !
    return

  end subroutine nearest_Ggrid

  !----------------------------------------------------------------------------

  subroutine    setup_lutelem(            &
                          gl,             & ! IN
                          gall,           & ! IN
                          undef,          & ! IN
                          tlon,           & ! IN
                          tlat,           & ! IN
                          lut_dim,        & ! IN
                          lut_elem,       & ! IN
                          lut_g,          & ! IN
                          rgnid,          & ! IN
                          rgnfilename,    & ! IN
                          inputfiledir,   & ! IN
                          offline,        & ! IN
                          opt_bin,        & ! IN
                          grid_num_array  & ! OUT
                             )
    !
    !
    ! ++ Document 
    !      Connect g-grid look up table and these latlon infomation
    !
    ! ++ Hisotry
    !     0.0     2018 3 20th    :  Takuya Kurihana
    !     0.1     2018 3 21th    :  Takuya Kurihana ; offline[true] test clear
    !
    !

    !use  mod_HSSOBS       , only : &
    !     get_icollarray

    use  mod_setupfname , only : &
         get_ifilename

    use  mod_distdeg    , only : &
         get_distance

    implicit none

    integer,        parameter      :: kind_real = 8

    integer,        intent(in)     :: gl
    integer,        intent(in)     :: gall
    real(kind_real),intent(in)     :: undef  ! assume minus value e.g. -9999.d0
    real(kind_real),intent(in)     :: tlon   ! target obs lon. info
    real(kind_real),intent(in)     :: tlat   ! target obs lat. info
    integer,        intent(in)     :: lut_dim
    integer,        intent(in)     :: lut_elem
    integer,        intent(in)     :: lut_g( lut_dim )
    integer,        intent(in)     :: rgnid
    character(256), intent(in)     :: rgnfilename
    character(256), intent(in)     :: inputfiledir
    logical,        intent(in)     :: offline   
    !                                 True   : read ico-latlon rgn file
    !                                 False  : receive ico-latlon array
    !
    logical,        intent(in)     :: opt_bin
    !                                 True   : read ico-latlon rgn binary file
    !                                 False  : read ico-latlon rgn text file
    !
    real(kind_real),intent(out)    :: grid_num_array( lut_dim , lut_elem )


    integer, parameter             :: log_id     = 300
    integer, parameter             :: icoll_dim  = 2
    integer, parameter             :: icolon_id  = 1
    integer, parameter             :: icolat_id  = 2

    real(kind_real)                :: icoll( gall, icoll_dim )
    !
    !                                 icoll_dim = 1 : longitude
    !                                 icoll_dim = 2 : latitude

    logical             :: judge       ! judge for cal distance
    character(128)      :: gridnumber
    character( 40)      :: prefix
    real(kind_real)     :: icolon
    real(kind_real)     :: icolat
    real(kind_real)     :: dist      ! [km]
    real(kind_real)     :: gl_dist   ! [km]
    integer             :: icodim
    integer             :: igdim
    integer             :: ielem
    integer             :: ignum, jgnum


    !--------------------------------------------------------------------------

    !++ set coef
    ! fixme here should be changed more reasonable!
    gl_dist = 224.d0 / real( gl - 5 + 1 )

    !------------------------------------------------------
    !++ Setup ico-latlon data
    !------------------------------------------------------
    if ( offline  ) then

      !----------------------------------------------------
      !++ Open & Read ico-latlon grid info
      !----------------------------------------------------
      call  get_ifilename(        &
                    rgnid,        & ! IN
                    gridnumber,   & ! OUT
                    prefix        & ! OUT
                    )

      if (  opt_bin  ) then
        open(log_id , &
            file=trim(inputfiledir)//'/'//trim(rgnfilename)//'.rgn'//trim(prefix)//trim(gridnumber), &
            form='unformatted' , access='direct', recl=kind_real*gall, action='read' )

             read( log_id, rec=1 ) icoll( : ,  icolon_id )
             read( log_id, rec=2 ) icoll( : ,  icolat_id )

        close(log_id)

      else if (  opt_bin .ne. .true.  ) then
        open(log_id , &
            file=trim(inputfiledir)//'/'//trim(rgnfilename)//'.rgn'//trim(prefix)//trim(gridnumber), &
            action='read' )

            do ignum = 1, gall
             read( log_id, * ) icoll( ignum ,  icolon_id ) , icoll( ignum ,  icolat_id )
            end do

        close(log_id)

      else
        write( 6, * ) '   !  Read Option Error /  Check opt_bin option  ! '
        write( 6, * ) '   Error at mod_gridsearch.f90 [ subroutine | setup_lutelem ] '
        stop
      end if
      !
    else  ! online
      !
      !----------------------------------------------------
      !++ Read ico-latlon array info
      !----------------------------------------------------
      ! fixme
      !  call  get_icollarray()
      !
      !   get icoll data read in setup_ico2ll
      !

    end if


    !------------------------------------------------------
    !++ Combine LUT and  ico-latlon data
    !------------------------------------------------------

    grid_num_array( : , : ) = undef
    do igdim =  1 , lut_dim
      
      ignum  = lut_g( igdim )
      if (  ignum > 0 ) then

          icolon = icoll(ignum,icolon_id )
          icolat = icoll(ignum,icolat_id )


          call  get_distance(         &
                      tlon ,          & ! IN
                      tlat ,          & ! IN
                      icolon ,        & ! IN
                      icolat ,        & ! IN
                      dist ,          & ! OUT
                      gl_dist,        & ! IN
                      judge           & ! OUT 
                            )

          ! store data into return array
          grid_num_array( igdim , 1 ) = dble( ignum )
          grid_num_array( igdim , 2 ) = icolon
          grid_num_array( igdim , 3 ) = icolat
          grid_num_array( igdim , 4 ) = dist
          write(6,*) grid_num_array( igdim , 4 )
      end if ! check lut_g undef

    end do   ! igdim

    return

  end subroutine setup_lutelem

  !----------------------------------------------------------------------------

  subroutine    setup_gridtable(          &
                          gl ,            & ! IN
                          gall ,          & ! IN
                          g_num,          & ! IN
                          rl,             & ! IN
                          lut_g           & ! OUT
                               )

    !
    ! ++ Document 
    !       make Look Up Table for a g-grid number
    !
    !       0.0  2018 3 20th    :  Takuya Kurihana
    !       0.1  2018 3 20th    :  All test case is clear by Takuya
    !

    implicit none

    integer,         parameter    :: kind_real      = 8
    integer,         parameter    :: lut_dim        = 9 ! Inner grid
    integer,         parameter    :: side_lut_dim   = 6 ! Side
    integer,         parameter    :: crnr_lut_dim   = 4 ! 4 corners
    real(kind_real), parameter    :: undef          = -9999.d0

    integer,        intent(in)    :: gl
    integer,        intent(in)    :: gall
    integer,        intent(in)    :: g_num   ! g grid number
    integer,        intent(in)    :: rl      ! rgn level
    !real(kind_real),intent(out)   :: lut_g(lut_dim)
    integer,        intent(out)   :: lut_g(lut_dim)
    !
    !
    !            Look Up  Table
    !        __________________
    !        |     |     |     |
    !        |     |     |     |
    !        |_____|_____|_____|   j ++
    !        |     |     |     |    /\
    !        |     |  g  |     |    ||
    !        |_____|_____|_____|    ||
    !        |     |     |     |    ||
    !        |     |     |     |
    !        |_____|_____|_____|
    !
    !
    !          ======> i ++
    !

    ! local variables
    logical             :: debug = .false.
    integer             :: gall_1d       ! num. of data in an grid / sqrt(gall)
    integer             :: location_type ! 
    integer             :: inc 
    integer             :: i, j, k
    integer             :: ignum , jgnum

    !--------------------------------------------------------------------------

    debug = .true.
    ! set relavant  g-grid values
    !gall_1d = int ( sqrt( gall ) )
    gall_1d = 2**(gl - rl)+2

    ! Initialize ; array filled out by undef
    lut_g(:) = undef

    ! setup ingrid / boundary types
    !
    ! Boundary ==> location_type = 0,1,2,3
    ! Inner    ==> location_type = 4
    !
    !
    if (      mod(g_num , gall_1d ) == 1  .and. &
              g_num /= 1  .and. g_num /= gall-gall_1d+1  )  then  ! Left  side boundary
        !
        location_type = 0
        !
    else if ( mod(g_num , gall_1d ) == 0  .and. &
              g_num /= gall_1d .and. g_num /= gall )  then ! Right side boundary
        !
        location_type = 1
        !
    else if ( g_num >= 1 .and. g_num <= gall_1d  ) then  ! Bottom side boundary
        location_type = 2
        !
    else if ( g_num >= gall - gall_1d + 1  .and. & 
              g_num <= gall                      )  then ! Upper  side boundary
        !
        location_type = 3
        !
    else ! In
        location_type = 4
        !
    end if ! In / Boundary outer loop

    if (debug) write(6,*) '  location_type  ' ,  location_type

    !------------------------------------------------------
    !++ generate Look Up Table
    !------------------------------------------------------
    if (       location_type >= 4 ) then

      inc = -1    ! increment consective grid moves
      do ignum = 1, lut_dim/3
        lut_g( ignum   ) = g_num - gall_1d + inc
        lut_g( ignum+3 ) = g_num + inc
        lut_g( ignum+6 ) = g_num + gall_1d + inc

        ! increment
        inc = inc + 1  
      end do
      !
      do i = 1, lut_dim
          if ( lut_g(i)  >  0  )  then
              write(*,*) '   ######## lut_g type 4 ',i, lut_g(i)
          end if
          if ( lut_g(i)  <= 0  .and. lut_g(i) > -3000 )  then
              stop
          end if
      end do
    else if (  location_type <= 3 ) then

      if (      location_type == 2  ) then
        !
        !  Bottom side boundary
        !

        !  inner part of side
        if ( g_num  /=  1 .and. g_num /= gall_1d  ) then
          inc = -1
          do ignum = 1, side_lut_dim/2
            lut_g( ignum   ) = g_num  + inc
            lut_g( ignum+3 ) = g_num  + gall_1d + inc

            ! increment
            inc = inc + 1
          end do

        !  left under corner
        else if (  g_num == 1)  then
          inc = 0
          do ignum = 1, crnr_lut_dim/2
            lut_g( ignum   ) = g_num  + inc
            lut_g( ignum+2 ) = g_num  + gall_1d + inc

            inc = inc + 1
          end do

        !  right under corner
        else if (  g_num == gall_1d)  then
          inc = -1
          do ignum = 1, crnr_lut_dim/2
            lut_g( ignum   ) = g_num  + inc
            lut_g( ignum+2 ) = g_num  + gall_1d + inc

            inc = inc + 1
          end do

        end if 


      else if ( location_type == 0  .and. & 
                 g_num /= 1 .and. g_num /= gall-gall_1d+1  ) then
        !
        !  Left  side boundary
        !    without left lower and left upper
        !
        inc = 0
        do ignum = 1, side_lut_dim/3
          lut_g( ignum    ) = g_num - gall_1d + inc
          lut_g( ignum+2  ) = g_num + inc
          lut_g( ignum+4  ) = g_num + gall_1d + inc

          ! increment 
          inc = inc + 1
        end do


      else if ( location_type == 1  .and. & 
                 g_num /= gall_1d   .and. g_num /= gall  ) then
        !
        !  Right  side boundary
        !    without right lower and right  upper
        !
        inc = -1
        do ignum = 1, side_lut_dim/3
          lut_g( ignum    ) = g_num - gall_1d + inc
          lut_g( ignum+2  ) = g_num + inc
          lut_g( ignum+4  ) = g_num + gall_1d + inc

          ! increment 
          inc = inc + 1
        end do


      else if ( location_type == 3  ) then
        !
        !  Upper side boundary
        !

        !  inner part of side
        if ( g_num  /=  gall-gall_1d+1  .and. g_num /= gall  ) then
          inc = -1
          do ignum = 1, side_lut_dim/2
            lut_g( ignum   ) = g_num  - gall_1d + inc
            lut_g( ignum+3 ) = g_num  + inc

            ! increment
            inc = inc + 1
          end do

        !  left upper corner
        else if (  g_num == gall-gall_1d+1 )  then
          inc = 0
          do ignum = 1, crnr_lut_dim/2
            lut_g( ignum   ) = g_num  - gall_1d + inc
            lut_g( ignum+2 ) = g_num  + inc

            inc = inc + 1
          end do

        !  right upper corner
        else if (  g_num == gall )  then
          inc = -1
          do ignum = 1, crnr_lut_dim/2
            lut_g( ignum   ) = g_num  - gall_1d + inc
            lut_g( ignum+2 ) = g_num  + inc

            inc = inc + 1
          end do

        end if

      else
        !
        write( 6, * ) '   Error Occurred in the IF loop of boundary LUT' 
        stop
        !
      end if ! boundary if 

    else
      ! Exception Case
      write(6, *) '  !  Classification Error Occurred  !  '
      write(6, *) '  Stop at [ subroutine |  setup_gridtable ]     '
      stop

    end if


    return

  end subroutine setup_gridtable

end module mod_gridsearch
