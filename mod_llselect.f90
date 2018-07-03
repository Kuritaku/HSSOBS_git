!------------------------------------------------------------------------------
!
!++ Document
!
!------------------------------------------------------------------------------
module mod_llselection
  implicit none 

  public  :: latlon_selection

  private :: get_latlon, int2char !, setup_TargetRgn

contains 
  !
  subroutine latlon_selection(   &
              gl, rl,            & ! in
              gall, lall,        & ! in
              GrdEDge,           & ! out
              gwflag             & ! inout
              )

    implicit none 

    integer, parameter    :: kind_real = 8
    integer, parameter    :: edgeall   = 5

    integer,           intent(in)     :: gl, rl
    integer,           intent(in)     :: gall, lall
    real(kind_real),   intent(out)    :: GrdEdge(:,:) ! [rgn][latLeft,lonL,latRight,lonR,lonMid]
    integer,           intent(inout)  :: gwflag(lall) 
                                        !  0  : not over Greenwitch line
                                        !  1  :     over Greenwitch line

    !integer               :: gall, lall
    !real(kind_real)       :: GrdEdge(5) ! [latMin,lonMin,latMax,lonMax,lonMid]
    !real(kind_real), allocatable  :: GrdEdge(:,:) ! [rgn][latLeft,lonL,latRight,lonR,lonMid]

    integer :: i,j


    !--------------------------------------------------------------------------
    !gall = (  2**(gl - rl )+2 )**2
    !lall = (  2**rl )**2 * 10

    !allocate(GrdEdge(lall,edgeall))

    ! get max-min latlon
    call get_latlon(      & 
              gl, rl,     &
              gall, lall, &
              edgeall,    &
              GrdEdge,    &
              gwflag      & ! inout
              )


    !do i = 1, lall
    !  do j = 1, edgeall
    !    write(6,*) GrdEdge(i,j)
    !  end do
    !end do

    ! setup latlon boundary
    !call latlon_setup

    ! Compute target rgn grid number
    !call setup_TargetRgn(   &
    !      gl, rl,           &
    !      gall, lall,       &
    !      GrdEdge,          &
    !      RgnNumber         &
    !      )

    !--> get new_lnum

  end subroutine latlon_selection

!------------------------------------------------------------------------------

  subroutine get_latlon(   &
              gl, rl,      &
              gall, lall,  &
              edgeall,     &
              GrdEdge,     &
              gwflag       &
              )
    implicit none 
    
    integer, parameter    :: kind_real = 8
    
    integer, intent(in)           :: gl, rl 
    integer, intent(in)           :: gall, lall 
    integer, intent(in)           :: edgeall 
    real(kind_real), intent(out)  :: GrdEdge(lall,edgeall)
    integer, intent(inout)        :: gwflag(lall)

    integer               :: filenum
    character(40)         :: prefix
    character(128)        :: gridnumber
    character(256)        :: inputinfofile
    character(256)        :: inputgridfile
    
    real(kind_real), allocatable       :: lon(:)
    real(kind_real), allocatable       :: lat(:)
    real(kind_real), allocatable       :: w1(:)
    real(kind_real), allocatable       :: w2(:)
    real(kind_real), allocatable       :: w3(:)
    real(kind_real), allocatable       :: lldata(:,:)
    real(kind_real), allocatable       :: ico_data(:,:)

    real(kind_real), allocatable       :: lonArray(:) ! tmp store lon data
    real(kind_real), allocatable       :: latArray(:) ! tmp store lat data

    integer, allocatable  :: lon_index(:)
    integer, allocatable  :: lat_index(:)
    integer, allocatable  :: n1_index(:)
    integer, allocatable  :: n2_index(:)
    integer, allocatable  :: n3_index(:)
    integer, allocatable  :: maxnum_latlon(:)
    integer, allocatable  :: nstat(:)
    integer, allocatable  :: nend(:)

    integer               :: igrid
    integer               :: imax
    integer               :: jmax
    integer               :: ArrayLength
    integer               :: i, j ,k, l, ij
    integer               :: ierr

    real(kind_real)       :: pi
    real(kind_real)       :: edge(edgeall)
    real(kind_real)       :: tmp_lon
    real(kind_real)       :: width_edge
    real(kind_real)       :: edge_length

    logical               :: debug = .false. ! write-option

    namelist / latlon_param / &
      inputinfofile,  &
      inputgridfile

    ! default settings
    filenum = 10
    inputinfofile = './llmap/llmap.info'
    inputgridfile = './llmap/llmap.rgn'

    pi = 4.0 * atan(1.0)

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------
    !++  namelist
    !--------------------------------------------------------------
    open(1, file='hssdriver.cnf' )
      read(1, nml=latlon_param, iostat=ierr)
    close(1)

    !--------------------------------------------------------------
    !++  read llmap info 
    !--------------------------------------------------------------
    open(filenum, file=trim(inputinfofile), form='unformatted' )

      ! read longitude 
      read(filenum) imax
        write(6,*) '   number of longitude  == ', imax
        allocate( lon(imax))
      read(filenum) lon(:)

      ! read latitude
      read(filenum) jmax
        write(6,*) '   number of latitude   == ', jmax
        allocate( lat(jmax))
      read(filenum) lat(:)

    close(filenum)

    !--------------------------------------------------------------
    !++  allocate relevant arraies
    !--------------------------------------------------------------
    ! llamp grid info
    allocate (lon_index(imax*jmax))
    allocate (lat_index(imax*jmax))
    allocate (n1_index(imax*jmax))
    allocate (n2_index(imax*jmax))
    allocate (n3_index(imax*jmax))
    allocate (w1(imax*jmax))
    allocate (w2(imax*jmax))
    allocate (w3(imax*jmax))
    
    ! iteration number
    allocate (maxnum_latlon(lall))
    allocate (nstat(lall))
    allocate (nend(lall))
    
    !--------------------------------------------------------------
    !++  selection process while each grid
    !--------------------------------------------------------------
    ! outer loop for rgn grid number 
    do igrid = 0, lall-1

        call int2char(    &
              igrid,      &
              gridnumber, &
              prefix      &
              )

        l = igrid + 1 ! conver for fortran array reference
        open(filenum, file=trim(inputgridfile)//trim(prefix)//trim(gridnumber), &
                      form='unformatted' )
            read(filenum) maxnum_latlon(l)
                nend(l)  = sum(maxnum_latlon(1:l))
                nstat(l) = nend(l) - maxnum_latlon(l) + 1

                if ( debug )  write(6,*)  nstat(l) , nend(l) , nend(l) - nstat(l) + 1
                ! read indexes from grid file
                if ( maxnum_latlon(l)  /= 0 ) then
                  read(filenum) lon_index(nstat(l):nend(l))
                  read(filenum) lat_index(nstat(l):nend(l))
                  read(filenum) n1_index(nstat(l):nend(l))
                  read(filenum) n2_index(nstat(l):nend(l))
                  read(filenum) n3_index(nstat(l):nend(l))
                  read(filenum) w1(nstat(l):nend(l))
                  read(filenum) w2(nstat(l):nend(l))
                  read(filenum) w3(nstat(l):nend(l))
                else
                  write(6,*) '   !! Inappropriate latlon length  !!  '
                  write(6,*) '   Stop Program for detecting 0 in Array Length '
                  stop
                end if

                ! tmp lat-lon grid
                ArrayLength = nend(l) - nstat(l) + 1
                allocate (lonArray(ArrayLength))
                allocate (latArray(ArrayLength))

                lonArray(:) = 0.00
                latArray(:) = 0.00
                do i = 1, ArrayLength
                  ij = nstat(l) + i -1
                  tmp_lon = lon(lon_index(ij)) * 180.0/pi

                  ! -180 to 180 ==> 0 to 360 
                  if ( tmp_lon < 0.d0 ) then
                       lonArray(i) = 360.d0 + lon(lon_index(ij)) * 180.0/pi
                  else 
                       lonArray(i) = lon(lon_index(ij)) * 180.0/pi
                  end if

                  latArray(i) = lat(lat_index(ij)) * 180.0/pi
                end do

                edge(1) = minval(latArray)
                edge(2) = minval(lonArray)
                edge(3) = maxval(latArray)
                edge(4) = maxval(lonArray)


                !------------------------------------------------------------ 
                ! ++  Check width if Greenwich line lies in mid of rgn
                !------------------------------------------------------------ 
                !
                width_edge = edge(4) - edge(2)
                edge_length =  360.d0 / real(lall) * 2.d0**(rl+1)

                ! No overline case
                if ( width_edge <  edge_length  ) then
                    edge(5) = ( maxval(lonArray)+ minval(lonArray) ) / 2.0
                    gwflag(l)  = 0

                ! Overline case
                else
                    !  fixme is it ok ??
                    do i = 1, ArrayLength
                       ij = nstat(l) + i -1
                       lonArray(i) = lon(lon_index(ij)) * 180.0/pi
                    end do

                    edge(2) = 360.d0 + minval(lonArray)
                    edge(4) = maxval(lonArray)
                    edge(5) = 0.00
                    gwflag(l)  = 1
                    write(6,*) '  Modify Edge Value due to the Greenwitch line '
                end if
                !------------------------------------------------------------ 

                ! store Edge info
                !write(6,*) '  Now Store ',l
                do j = 1, edgeall
                    GrdEdge(l,j) = edge(j)
                    !write(6,*) GrdEdge(l,j)
                end do

        close(filenum)
        filenum = filenum + 1


        ! tmp lat-lon grid
        deallocate (lonArray)
        deallocate (latArray)

    end do !

    !--------------------------------------------------------------
    !++  deallocate relevant arraies
    !--------------------------------------------------------------
    ! llamp grid info
    deallocate (lon_index)
    deallocate (lat_index)
    deallocate (n1_index)
    deallocate (n2_index)
    deallocate (n3_index)
    deallocate (w1)
    deallocate (w2)
    deallocate (w3)

    ! iteration number
    deallocate (maxnum_latlon)
    deallocate (nstat)
    deallocate (nend)


    return
  end subroutine

!------------------------------------------------------------------------------

  subroutine int2char(    &
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
