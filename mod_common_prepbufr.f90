module mod_common_prepbufr
  implicit none

 ! ID for observations
 integer, parameter, public    :: id_u_obs=2819
 integer, parameter, public    :: id_v_obs=2820
 integer, parameter, public    :: id_t_obs=3073
 integer, parameter, public    :: id_q_obs=3330
 integer, parameter, public    :: id_rh_obs=3331
 integer, parameter, public    :: id_ps_obs=14593

 ! # Retuen function for prepbufr unique variable numbers
 public   :: get_prepbufr_varID 

contains
  !
  function  get_prepbufr_varID( pb_id )&
      result( vid )
      !
      implicit none
      integer :: pb_id  !--IN
      !          id
      !           1 : u
      !           2 : v
      !           3 : temp
      !           4 : q
      !           5 : rh
      !           6 : ps
      integer :: vid    !--OUT

      if (      pb_id == 1  ) then
              vid = id_u_obs
      else if(  pb_id == 2  ) then
              vid = id_v_obs
      else if(  pb_id == 3  ) then
              vid = id_t_obs
      else if(  pb_id == 4  ) then
              vid = id_q_obs
      else if(  pb_id == 5  ) then
              vid = id_rh_obs
      else if(  pb_id == 6  ) then
              vid = id_ps_obs
      end if

  end function
  !
end module mod_common_prepbufr
