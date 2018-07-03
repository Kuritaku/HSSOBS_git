!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
module mod_setupfname

  public :: get_ifilename

contains 
  !
  !----------------------------------------------------------------------------
  !
  subroutine get_ifilename(    &
              igrid,      &
              gridnumber, &
              prefix      &
              )

    implicit none 

    integer,        intent(in)      :: igrid
    character(128), intent(out)     :: gridnumber
    character(40),  intent(out)     :: prefix

    !--------------------------------------------------------------------------

     if (      igrid < 10  ) then
       write(gridnumber, "(i1)") igrid
       prefix = '0000'

     else if ( igrid < 100   .and. igrid >= 10 ) then
       write(gridnumber, "(i2)") igrid
       prefix = '000'

     else if ( igrid < 1000  .and. igrid >= 100 ) then
       write(gridnumber, "(i3)") igrid
       prefix = '00'

     else if ( igrid < 10000 .and. igrid >= 1000 ) then
       write(gridnumber, "(i4)") igrid
       prefix = '0'

     end if 

    return
  end subroutine

end module

