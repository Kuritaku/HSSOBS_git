!------------------------------------------------------------------------------
!
! ++ module class for quick sort with ID return
!
! ++ hisotry
!     0.0   2018, March 19th  Takuya Kurihana
!
!------------------------------------------------------------------------------
module mod_qs
  implicit none

  public  :: main_qs

  private :: sorted_qs

contains
  !
  subroutine main_qs(       &
             nall,          & ! IN  ! Number of data length
             array_dim,     & ! IN  ! Array dimension ( nall , array_dim )
             array_index,   & ! IN  ! Index to sort
             get_number,    & ! IN  ! Numeber of return vals
             array,         & ! IN
             ! Add Me!! DESC_order (logical) descending order
             array_tops     & ! OUT
                    )
    implicit none 

    integer, parameter    :: kind_real = 8
    integer, intent(in)   :: nall
    integer, intent(in)   :: array_dim
    integer, intent(in)   :: array_index
    integer, intent(in)   :: get_number

    real(kind_real), intent(in)   :: array(nall, array_dim) ! assume ( num. of data, elem )
    real(kind_real), intent(out)  :: array_tops(get_number, array_dim) ! assume ( num. of selected data, elem )


    logical           :: flag
    logical           :: debug  = .false.
    !real(kind_real)   :: array_1d(nall)
    real(kind_real)   :: array_2d(nall,array_dim)
    integer           :: istat
    integer           :: iend
    integer           :: restart
    integer           :: icount  ! If IN count 
    integer           :: olcount ! Outer Loop count
    integer           :: ii, jj
    integer           :: i,j,k

    !--------------------------------------------------------------------------

    !--------------------------------
    ! initialization
    !--------------------------------
    istat = 1
    !icount = 0
    iend  = nall
    flag = .true.
    !--------------------------------

    ! Extract an index data to reshape 1D array
    !array_1d(:) = array(:,array_index)

    ! copy all data for computation
    array_2d(:,:) = array(:,:)

    istat   = 1
    iend    = nall
    olcount = -1
    do while ( flag )
      !
      !  first sorting
      !
      !   ----------------------------
      !   |         2D  Array        |
      !   ----------------------------
      !
      !               ||
      !              \  /
      !               \/
      !
      !   ----------------------------
      !   |    Left    |    Right    |
      !   ----------------------------
      !
      call  sorted_qs(        &
                nall,         & ! IN
                array_dim,    & ! IN
                array_index,  & ! IN
                istat, iend,  & ! IN
                ii, jj,       & ! OUT
                array_2d      & ! INOUT
                     )

      ! judge first outer loop or not
      if ( olcount == -1 ) then
          restart  = ii
          olcount  = 0
      end if

      ! Left shufle
      if ( 1 < ii - 1 .and. icount == 0 ) then
        istat = 1
        iend  = ii - 1
        call  sorted_qs(        &
                  nall,         & ! IN 
                  array_dim,    & ! IN
                  array_index,  & ! IN
                  istat, iend,  & ! IN
                  ii, jj,       & ! OUT
                  array_2d      & ! INOUT
                       )
      else
        icount = 1
      end if

      ! Right shufle
      if ( nall > jj + 1 .and. icount == 1) then
        istat = jj + 1
        iend  = nall
        call  sorted_qs(        &
                  nall,         & ! IN 
                  array_dim,    & ! IN
                  array_index,  & ! IN
                  istat, iend,  & ! IN
                  ii, jj,       & ! OUT
                  array_2d      & ! INOUT
                       )
      else
        ! should pass here twice !
        olcount = olcount + 1
        if (debug) write(6,*) ' olcount ++ ', olcount
      end if

      ! Break Requisite
      if ( olcount  <= 1 .and. icount == 1 ) then
          istat   = restart
          iend    = nall
          icount  = 0
      else if  ( olcount >= 2 .and. icount == 1) then
          flag = .false.
      end if

    !
    end do ! flag loop


    !  check all data-set
    !write(6,*)  '   CIO '
    !write(6,*)  ( array_2d(i, array_index) , i=1,nall)
    !write(6,*) 'get_number , array_dim , array_index'
    !write(6,*)  get_number , array_dim , array_index

    !  select number of data from all data-set
    do i = 1, get_number
      array_tops(i, 1:array_dim) = array_2d(i, 1:array_index)
    end do


  end subroutine main_qs

  !----------------------------------------------------------------------------

  subroutine sorted_qs(     &
              n,            & ! IN  number of all data
              ndim,         & ! IN  dimension of array
              id,           & ! IN  sort target dimension
              istat, iend,  & ! IN
              iout, jout,   & ! OUT
              array_2d      & ! INOUT
                      )
    !
    implicit none

    integer, parameter             :: kind_real = 8
    integer, intent(in)            :: n
    integer, intent(in)            :: ndim
    integer, intent(in)            :: id
    integer, intent(in)            :: istat
    integer, intent(in)            :: iend
    integer, intent(out)           :: iout
    integer, intent(out)           :: jout

    real(kind_real), intent(inout) :: array_2d(n,ndim)

    logical           :: flag = .true.  ! flag to break loop
    integer           :: iistat
    integer           :: iiend
    integer           :: ii, jj
    integer           :: i,j,k
    integer           :: pivot
    integer           :: alength
    real(kind_real)   :: pval  ! pivot value
    real(kind_real)   :: tmp(1,ndim)   ! tmp value for swap : 2d

    !--------------------------------------------------------------------------
    alength = iend - istat + 1
    !----------------------------------
    ! ++ set up pivot & pval
    !----------------------------------
    if (  mod( alength ,2) == 0 ) then  ! even
          pivot = int((istat + iend)/2) + 1
    else                       ! odd
          pivot = (istat + iend)/2
    end if

    !write(6,*) 'pivot', pivot
    !----------------------------------
    ! ++ start sorting
    !----------------------------------

    ! Duplicate IN vars for iteration
    iistat    = istat
    iiend     = iend

    ii    = istat
    jj    = iend
    pval  = array_2d(pivot,id)
    do while ( flag )

      ! left hand side
      do while( pval  > array_2d(ii, id))
        ii = ii + 1
      end do

      ! right hand side
      do while( pval  < array_2d(jj, id))
        jj = jj - 1
      end do

      ! break loop
      if ( ii >= jj ) then
        flag = .false.
      end if

      if (  flag  ) then

        ! swap
        tmp(  1     , 1:ndim)  = array_2d(ii, 1:ndim)
        array_2d(ii , 1:ndim)  = array_2d(jj, 1:ndim)
        array_2d(jj , 1:ndim)  = tmp( 1     , 1:ndim)

        ! prep id values
        iistat = iistat + 1
        iiend  = iiend  - 1
        ii = iistat
        jj = iiend
        !write(6,*) ' ii jj in loop  ', ii, jj
        !write(6,*) ' pval in loop   ', pval

        !write(6,*)  (array_1d(k), k=1,n)
        !write(6,*)  (array_2d(k, 1), k=1,n)
        !write(6,*) ' ==== checkio array_2d in  qs loop ===='
        !write(6,*)  array_2d
        !write(6,*) '  '
      end if

    end do ! outer flag

    flag = .true.

    ! output
    iout = ii
    jout = jj
    return
  end subroutine

end module
