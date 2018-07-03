!------------------------------------------------------------------------------
!
! ++ Document 
!     Module class to calculate distance between 2 points &
!                  to return delta lat-lon degree by a distance
!
!------------------------------------------------------------------------------
module mod_distdeg
  implicit none

  public   ::  get_distance !,         &
               !get_distance_HighRes,  & future add 
               !get_deltadeg

  !private  ::

contains
  !
  subroutine get_distance(        &
                   tlon,          &  ! IN
                   tlat,          &  ! IN
                   ico_lon,       &  ! IN
                   ico_lat,       &  ! IN
                   dist,          &  ! OUT
                   gl_dist,       &  ! IN
                   judge          &  ! OUT
                         )
    implicit none

    integer,          parameter   :: kind_real = 8
    real(kind_real),  intent(in)  :: tlon  ! target data longitude
    real(kind_real),  intent(in)  :: tlat  ! target data latitude
    real(kind_real),  intent(in)  :: ico_lon  ! nicam ico data longitude
    real(kind_real),  intent(in)  :: ico_lat  ! nicam ico data latitude
    real(kind_real),  intent(out) :: dist     ! distance betw. A and B
    real(kind_real),  intent(in)  :: gl_dist ! glevel rsolution distance[km]
    logical,          intent(out) :: judge    ! eps_dist > dist => .true.

    
    real(kind_real),parameter     :: r = 6378.137000d0 ! [km]
    real(kind_real)   :: pi
    real(kind_real)   :: thres_dist
    real(kind_real)   :: deg_coef

    real(kind_real)   :: dlon
    real(kind_real)   :: dlat
    real(kind_real)   :: iphi , jphi
    real(kind_real)   :: dFunc

    ! coefficient
    pi          = 4.d0 * datan(1.d0)
    thres_dist  = gl_dist*0.5*dsqrt(2.d0) ! limit length of distance
    !
    !        __________ _________
    !       |         |        *|\      *** = diagonal gl_distance
    !       |         |      *  | \
    !       |         |    *    |  \         = 1/2 * gl_dist * sqrt(2)
    !       |         |  *      |   \
    !       |_________|*________| __ }  Half gl_diatance
    !       {     gl_distance   }
    !
    !
    deg_coef = 1.d0 / 180.d0

    !--------------------------------------------------------------------------

    dlon = ( tlon - ico_lon ) * pi *  deg_coef

    ! coef
    iphi = tlat    * pi *  deg_coef
    jphi = ico_lat * pi *  deg_coef

    dFunc = dsin(iphi)*dsin(jphi) + dcos(iphi)*dcos(jphi)*dcos(dlon)

    ! avoid NaN in arccos function
    if (  dFunc >= 1.d0  ) then
      dFunc = dFunc - 1.0d-12
    else if (  dFunc <= -1.d0  ) then
      dFunc = dFunc + 1.0d-12
    end if

    dist  =  r * dacos(dFunc)

    if ( dist <= thres_dist ) then
        judge = .true.
    else
        judge = .false.
    end if

    return

  end subroutine get_distance

  !----------------------------------------------------------------------------

  ! fixme future add nonlinear optimization to decide epaval parameters
  !subroutine get_deltadeg(       &
  !                       )
  !  implicit none

    !--------------------------------------------------------------------------

  !  return

  !end subroutine get_deltadeg


end module mod_distdeg
