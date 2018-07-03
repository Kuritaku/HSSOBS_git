!------------------------------------------------------------------------------
!
!++ Document
!       Save and Store grid, height data & Average them
!
!++ Hisotry 
!     0.0   2018 03 23th    Takuya Kurihana : Code alpha version
!     1.0   2018 06 09th    Takuya Kurihana : Add Prepbufr Subroutine
!
!------------------------------------------------------------------------------
module mod_compave
  implicit none

  public  :: setup_ave,           & ! setup array 
             get_aveval,          &
             !
             ! Prepbufr Obs data
             setup_obs_ave,       &
             get_obs_aveval

  !public  :: comp_arithmean , &
  !           comp_weightedAve

contains
  !
  subroutine setup_obs_ave(               &
                            g_id,         & !<-- IN
                            rgn_id,       &
                            zlev_id,      &
                            gall,         &
                            lall,         &
                            zlayer,       &
                            varval,       &
                            errval,       &
                            pbtype,       &
                            ave_box,      & !<-- OUT
                            ave_count,    &
                            ave_err_box,  &
                            ave_err_count &
                          )
    !
    !++ module load
    use  mod_hssconst,  only : &
         kind_real, &
         nobtype

    implicit none

    !integer,    parameter          ::  kind_real = 8
    integer,         intent(in)    :: g_id
    integer,         intent(in)    :: rgn_id
    integer,         intent(in)    :: zlev_id
    integer,         intent(in)    :: gall
    integer,         intent(in)    :: lall
    integer,         intent(in)    :: zlayer
    real(kind_real), intent(in)    :: varval
    real(kind_real), intent(in)    :: errval
    real(kind_real), intent(in)    :: pbtype
    !real(kind_real), intent(in)    :: dist_mini(elem_id)
    real(kind_real), intent(inout) :: ave_box(   gall , zlayer , 0:lall-1, nobtype ) ! save
    real(kind_real), intent(inout) :: ave_count( gall , zlayer , 0:lall-1, nobtype ) ! save
    real(kind_real), intent(inout) :: ave_err_box(   gall , zlayer , 0:lall-1, nobtype ) ! save
    real(kind_real), intent(inout) :: ave_err_count( gall , zlayer , 0:lall-1, nobtype ) ! save

    real(kind_real)   :: tmp_val
    real(kind_real)   :: tmp_count
    real(kind_real)   :: tmp_err_val
    real(kind_real)   :: tmp_err_count
    integer           :: obtype
    !--------------------------------------------------------------------------

    ! Increment variable's value
    obtype = int( pbtype )
    tmp_val     = ave_box( g_id , zlev_id , rgn_id, obtype )
    tmp_err_val = ave_err_box( g_id , zlev_id , rgn_id, obtype )
    ave_box( g_id , zlev_id , rgn_id, obtype )     = tmp_val      + varval
    ave_err_box( g_id , zlev_id , rgn_id, obtype ) = tmp_err_val  + errval

    ! count ++
    tmp_count = ave_count( g_id , zlev_id , rgn_id, obtype )
    ave_count( g_id , zlev_id , rgn_id, obtype ) = tmp_count +1.d0
    !
    tmp_err_count = ave_err_count( g_id , zlev_id , rgn_id, obtype )
    ave_err_count( g_id , zlev_id , rgn_id, obtype ) = tmp_err_count +1.d0

    return

  end subroutine setup_obs_ave
  !
  !----------------------------------------------------------------------------
  !
  subroutine get_obs_aveval(               &
                             gall,         &
                             lall,         &
                             zlayer,       &
                             ave_box,      &
                             ave_count,    &
                             ave_err_box,  &
                             ave_err_count &
                           )
    !
    !++ module load
    use  mod_hssconst,  only : &
         kind_real, &
         nobtype


    implicit none

    !integer,  parameter     ::  kind_real = 8
    integer,  parameter     ::  log_id    = 1000

    integer,  intent(in)    :: gall
    integer,  intent(in)    :: lall
    integer,  intent(in)    :: zlayer
    real(kind_real),  intent(inout) :: ave_box(   gall , zlayer , 0:lall-1, nobtype ) ! save
    real(kind_real),  intent(inout) :: ave_count( gall , zlayer , 0:lall-1, nobtype ) ! save
    real(kind_real),  intent(inout) :: ave_err_box(   gall , zlayer , 0:lall-1, nobtype ) ! save
    real(kind_real),  intent(inout) :: ave_err_count( gall , zlayer , 0:lall-1, nobtype ) ! save

    real(kind_real)         :: sum_count
    real(kind_real)         :: sum_val
    real(kind_real)         :: sum_errcount
    real(kind_real)         :: sum_errval

    integer     :: ignum
    integer     :: ilnum
    integer     :: iznum
    integer     :: iobtype

    !--------------------------------------------------------------------------

    do iobtype = 1, nobtype
      do ilnum = 0, lall-1
        do iznum = 1, zlayer
          do ignum = 1, gall

            ! obsdata
            sum_count = ave_count( ignum , iznum , ilnum, iobtype )
            if (   sum_count /= 0.d0 ) then
              sum_val = ave_box(  ignum , iznum , ilnum, iobtype )
              write(*,*) ' sum ', sum_val , sum_count
              ave_box(  ignum , iznum , ilnum, iobtype  ) =  sum_val / sum_count
            end if
            !
            ! obserr
            sum_errcount = ave_err_count( ignum , iznum , ilnum, iobtype )
            if (   sum_errcount /= 0.d0 ) then
              sum_errval = ave_err_box(  ignum , iznum , ilnum, iobtype )
              !write(*,*) ' sum ', sum_val , sum_count
              ave_err_box(  ignum , iznum , ilnum, iobtype  ) =  sum_errval / sum_errcount
            end if
          end do
        end do
      end do
    end do

    return

  end subroutine get_obs_aveval

  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !
  subroutine setup_ave(            &
                      !gl,          &
                      !rl,          &
                      g_id,        &
                      rgn_id,      &
                      zlev_id,     &
                      !elem_id,     &
                      gall,        &
                      lall,        &
                      zlayer,      &
                      varval,      &
                      !dist_mini,   &
                      ave_box,     &
                      ave_count    &
                      )
    implicit none

    integer,    parameter          ::  kind_real = 8
    !integer,    parameter          ::  ignum     = 1  ! Dim. ID Num. of dist_mini 

    !integer,    intent(in)         :: gl
    !integer,    intent(in)         :: rl
    integer,    intent(in)         :: g_id
    integer,    intent(in)         :: rgn_id
    integer,    intent(in)         :: zlev_id
    !integer,    intent(in)         :: elem_id
    integer,    intent(in)         :: gall
    integer,    intent(in)         :: lall
    integer,    intent(in)         :: zlayer
    real(kind_real), intent(in)    :: varval
    !real(kind_real), intent(in)    :: dist_mini(elem_id)
    real(kind_real), intent(inout) :: ave_box(   gall , zlayer , 0:lall-1 ) ! save
    real(kind_real), intent(inout) :: ave_count( gall , zlayer , 0:lall-1 ) ! save


    real(kind_real)   :: tmp_val
    real(kind_real)   :: tmp_count

    !integer           :: gnum

    !--------------------------------------------------------------------------

    !gnum  = int( dist_mini( ignum ) )

    ! Increment variable's value
    tmp_val = ave_box( g_id , zlev_id , rgn_id )
    ave_box( g_id , zlev_id , rgn_id ) = tmp_val + varval

    ! count ++
    tmp_count = ave_count( g_id , zlev_id , rgn_id )
    ave_count( g_id , zlev_id , rgn_id ) = tmp_count +1.d0

    return

  end subroutine setup_ave
  !
  !----------------------------------------------------------------------------
  !
  subroutine get_aveval(            &
                        gall,       &
                        lall,       &
                        zlayer,     &
                        ave_box,    &
                        ave_count   &
                       )
    implicit none

    integer,  parameter     ::  kind_real = 8
    integer,  parameter     ::  log_id    = 1000

    integer,  intent(in)    :: gall
    integer,  intent(in)    :: lall
    integer,  intent(in)    :: zlayer
    real(kind_real),  intent(inout) :: ave_box(   gall , zlayer , 0:lall-1 ) ! save
    real(kind_real),  intent(inout) :: ave_count( gall , zlayer , 0:lall-1 ) ! save

    real(kind_real)         :: sum_count
    real(kind_real)         :: sum_val

    integer     :: ignum
    integer     :: ilnum
    integer     :: iznum

    !--------------------------------------------------------------------------

    do ilnum = 0, lall-1
      do iznum = 1, zlayer
        do ignum = 1, gall
          sum_count = ave_count( ignum , iznum , ilnum )
          if (   sum_count /= 0.d0 ) then
            sum_val = ave_box(  ignum , iznum , ilnum )
            write(6,*) ' sum ', sum_val , sum_count
            ave_box(  ignum , iznum , ilnum  ) =  sum_val / sum_count
          end if
        end do
      end do
    end do

    return

  end subroutine get_aveval

  !----------------------------------------------------------------------------

  !subroutine comp_arithmean(        &
  !                         )

  !end subroutine comp_arithmean

  !----------------------------------------------------------------------------

  !subroutine comp_weightedAve(        &
  !                           )
  !
  !end subroutine  comp_weightedAve

end module mod_compave
