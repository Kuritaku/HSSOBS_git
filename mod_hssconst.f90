!------------------------------------------------------------------------------
!
!++  Document 
!    Main Program of High Speed Super-OBservation System
!
!++  Hisotry
!
!       Version     Date       Editer        Description
!    ----------------------------------------------------------------
!        0.0     2018-05-12  T. Kurihana   Gather HSSOBS parameters
!        1.0     2018-05-12  T. Kurihana   Add Prepbufr  parameters
!
!------------------------------------------------------------------------------
!
module mod_hssconst
  implicit none


  integer,         parameter ,   public  :: kind_real  = 8
  integer,         parameter ,   public  :: kind_sngl  = 4
  integer,         parameter ,   public  :: kind_char  = 256
  !integer,         parameter   :: edgeall    = 5
  !integer,         parameter   :: tgdim      = 2 ! dim. of target grid
  integer,         parameter ,   public  :: lut_dim    = 9
  integer,         parameter ,   public  :: lut_elem   = 4
  integer,         parameter ,   public  :: rgmin_elem = 5
  integer,         parameter ,   public  :: sort_id    = 4
  integer,         parameter ,   public  :: get_number = 1
  !                sort_id  dim 1 : grid number
  !                         dim 2 : longitude
  !                         dim 3 : latitude
  !                         dim 4 : distance
  integer,         parameter ,   public  :: dist_id    = 9  ! number of return data
  integer,         parameter ,   public  :: Array_Dim  = 4
  integer,         parameter ,   public  :: tab_id     = 1
  integer,         parameter ,   public  :: tmp_rgndim = 3
  integer,         parameter ,   public  :: RgnCandMaxID = 9
  integer,         parameter ,   public  :: grd_Xdim   = 1
  integer,         parameter ,   public  :: grd_Ydim   = 2
  integer,         parameter ,   public  :: grd_Zdim   = 3
  integer,         parameter ,   public  :: grd_Pdim   = 3
  integer,         parameter ,   public  :: nobtype    = 20
  integer,         parameter ,   public  :: td_dim     = 2
  !
  ! Pgrid
  integer,         parameter ,   public  :: pall       = 19


end module
