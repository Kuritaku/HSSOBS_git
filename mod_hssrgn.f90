!------------------------------------------------------------------------------
!
!++  Document 
!
!
!++  Hisotry
!
!       Version    Date       Editer        Description
!    ----------------------------------------------------------------
!        0.0      2018-04-06  T. Kurihana   Building RGN Loop
!
!
!++  API
!
!------------------------------------------------------------------------------
module  mod_hssrgn
  implicit none

  public  ::  RgnLoop

contains
  !
  subroutine  RgnLoop(                &
                      gl,             & ! IN
                      rl,             & ! IN
                      gall,           & ! IN
                      tmp_grid,       & ! IN
                      !lall,         &
                      rgnid_dim,      & ! IN
                      offline,        & ! IN
                      opt_bin,        & ! IN
                      rgnfilename,    & ! IN
                      inputfiledir,   & ! IN
                      RgnCandArray,   & ! INOUT
                      !In_dist_array,     & ! INOUT
                      !In_CandDist_mini,  & ! INOUT
                      rgngrid_mini    & ! OUT
                     )

    !
    ! ++ Used Modules
    !
    use  mod_deltall,      only :  &  ! loughly compute nearest g_number
         deltalatlon

    use  mod_getrgnid,     only :  &  ! set rgnid
         setup_targetrgn

    use  mod_gridsearch,   only :  &  ! grid search from a rgn-grid
         nearest_Ggrid,            &
         setup_lutelem,            &
         setup_gridtable

    use  mod_qs ,          only :  &  ! quick sort
         main_qs


    implicit none
  
    integer,         parameter     :: kind_real  = 8
    integer,         parameter     :: kind_char  = 256
    integer,         parameter     :: lut_dim    = 9
    integer,         parameter     :: lut_elem   = 4
    integer,         parameter     :: rgmin_elem = 5
    integer,         parameter     :: sort_id    = 4
    !                sort_id  dim 1 : grid number
    !                         dim 2 : longitude
    !                         dim 3 : latitude
    !                         dim 4 : distance
    integer,         parameter     :: get_number = 1
    integer,         parameter     :: dist_id    = 9  ! number of return data
    integer,         parameter     :: Array_Dim  = 4
    integer,         parameter     :: tab_id     = 1
    integer,         parameter     :: tmpvar_dim = 2
    integer,         parameter     :: grd_Xdim   = 1
    integer,         parameter     :: grd_Ydim   = 2

    real(kind_real), parameter     :: LUTundef   = -9999.d0

    integer,         intent(in)    :: gl
    integer,         intent(in)    :: rl
    integer,         intent(in)    :: gall
    integer,         intent(in)    :: rgnid_dim

    real(kind_real), intent(in)    :: tmp_grid(tmpvar_dim)

    logical,         intent(in)    :: offline
    logical,         intent(in)    :: opt_bin

    character(kind_char), intent(in)  :: rgnfilename
    character(kind_char), intent(in)  :: inputfiledir

    !integer,         intent(inout) :: RgnCandArray( rgnid_dim  )
    !integer,         intent(in)    :: RgnCandArray( rgnid_dim  )
    integer,       intent(inout)    :: RgnCandArray( rgnid_dim  )
    !real(kind_real), intent(in)    :: In_dist_array(   dist_id   , Array_Dim,tab_id)
    !real(kind_real), intent(in)    :: In_CandDist_mini(rgnid_dim , rgmin_elem)
    real(kind_real), intent(out)   :: rgngrid_mini( rgmin_elem )


    real(kind_real)                :: garray(lut_dim, lut_elem)
    real(kind_real)                :: dist_mini(lut_elem)
    real(kind_real)                :: dist_array(   dist_id   , Array_Dim,tab_id)
    real(kind_real)                :: CandDist_mini(rgnid_dim , rgmin_elem)
    real(kind_real)                :: ctlon
    real(kind_real)                :: ctlat

    integer                        :: lut_g(lut_dim)
    integer                        :: rgn_tabs(tab_id)
    integer                        :: l
    integer                        :: rgnid
    integer                        :: ielem
    integer                        :: i, j, k
    integer                        :: g_num

    !--------------------------------------------------------------------------

    ! move value
    ctlon = tmp_grid(grd_Xdim)
    ctlat = tmp_grid(grd_Ydim)

    !dist_array( : , : , : ) = In_dist_array( :, : , :)
    !CandDist_mini( : , : )  = In_CandDist_mini( :, :)

    do l = 1 , rgnid_dim
       rgnid = RgnCandArray(l)
       rgn_tabs(1) = rgnid
       dist_array( : , : , : ) = 0.d0
       write(*,*) '  Current Processing Rgn ID ==  ', rgnid

      ! mod_dlatlon.f90
      call    deltalatlon(   &
                gl, rl,       & ! IN
                gall,         & ! IN
                tmp_grid,     & ! a temporarly evaluated obs info
                tab_id,       & ! ID? Dim of (a) potential rgn grids(grid)
                rgn_tabs,     & ! array of potential rgn grids
                opt_bin,      & ! IN
                dist_array    & ! OUT
                )

      g_num = int(dist_array( 1 , sort_id , tab_id )) ! get g-grid number
      write(6,*)  '  estimated g_num == ' , g_num
      if ( g_num < 0  ) then
          write(*,*) dist_id, sort_id, tab_id
          write(911, *) dist_array(:,sort_id,tab_id)
          write(912, *) dist_array(:,:,:)
          write(912,*) '-------------------------------'
      end if

      ! First set up garray
      call setup_gridtable(   &
            gl,               &
            gall,             &
            g_num,            &
            rl,               &
            lut_g             &
                              )

      write(6,*) lut_g


      ! check distance
      call  setup_lutelem(            &
                gl,             & ! IN
                gall,           & ! IN
                LUTundef,       & ! IN
                ctlon,          & ! IN
                ctlat,          & ! IN
                lut_dim,        & ! IN
                lut_elem,       & ! IN
                lut_g,          & ! IN
                rgnid,          & ! IN
                rgnfilename,    & ! IN
                inputfiledir,   & ! IN
                offline,        & ! IN
                opt_bin,        & ! IN
                garray          & ! OUT
                     )


    ! view table
    !do i = 1, lut_dim
    !    write(6,'(4f12.6)')  ( garray(i , j), j = 1, lut_elem )
    !end do
    !end do! rgn loop : old

    !dist_mini(:) = 0.d0
    !write(6,*) garray
    !write(6,*) lut_dim , lut_elem , sort_id

      call nearest_Ggrid(      &
              lut_dim,          &
              lut_elem,         &
              sort_id,          &
              garray,           &
              dist_mini         &
                         )

      ! bug check for more than one call mod_qs
      write(*,*)'  '
      write(*,*)'    ========== View Searching Process  =============      '
      write(6,"(5A12)") 'id ', 'g_num'  , 'lon' , 'lat' , 'dist'
      write(6, '(1I10,1A, 4f12.6)') rgnid  ,'  ' , ( dist_mini(j), j=1, lut_elem) 
      write(*,*)'  '
      !write(6,*) ''
      !write(6,*) ''
      !stop

      ! fixme compare 

      do ielem = 1,  rgmin_elem
        if (  ielem  <  rgmin_elem ) then
            CandDist_mini( l , ielem ) = dist_mini(ielem) ! id 1 - 4
        elseif (  ielem == rgmin_elem) then
            CandDist_mini( l , ielem ) = dble(rgnid)  ! id 5
        end if
      end do

    end do! rgn loop


    ! ======> until here search potential grid number
    !----------------------------------------------------------------------------

    ! fixme add sub. to decide dupulicated grid  
    ! quick sort to get number
    call    main_qs(                  &
                     rgnid_dim,       & ! IN  ! Number of data length
                     rgmin_elem,      & ! IN  ! Ar
                     sort_id,         & ! IN  ! Index to sort
                     get_number,      & ! IN  ! Numeber of return vals
                     CandDist_mini,   & ! IN
                     rgngrid_mini     & ! OUT
                    )

     write(6,*) '         Minimun         '
     write(6,"(5A12)") 'g_num'  , 'lon' , 'lat' , 'dist' , ' rgnid'
     write(6, '(5f12.6)') ( rgngrid_mini(j), j=1, rgmin_elem) 
     !write(600, '(5f12.6)') ( rgngrid_mini(j), j=1, rgmin_elem) 
     write(6,*) ''


    ! ======> until here grt  grid  rgn number
    !----------------------------------------------------------------------------

    return
  end subroutine

end module 
