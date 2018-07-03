#!/bin/bash


###############################################################################
#
#   =============               User Settings             =================
#
###############################################################################

# nicam
gl=5
rl=0
zall=38

#decoded prepbudfr  date
export yyyy=2011
export mm=11
export dd=01
export hh=00
ymd=${yyyy}${mm}${dd}

# prepbufr varid from 1-6
nobs_prepbufr=11380
var_id=6
#varname='u' # 1
#varname='v' # 2
#varname='temp' # 3
#varname='q' # 4
varname='ps' # 6

num_of_date=7
ODATE1=` date -d "${ymd} ${hh}:00 3 hour ago " +%Y%m%d%H `
ODATE2=` date -d "${ymd} ${hh}:00 2 hour ago " +%Y%m%d%H `
ODATE3=` date -d "${ymd} ${hh}:00 1 hour ago " +%Y%m%d%H `
ODATE4=${ymd}${hh}
ODATE5=` date -d "${ymd} ${hh}:00 1 hour  " +%Y%m%d%H `
ODATE6=` date -d "${ymd} ${hh}:00 2 hour  " +%Y%m%d%H `
ODATE7=` date -d "${ymd} ${hh}:00 3 hour  " +%Y%m%d%H `

#ODATE_LIST=( ${ODATE1} ${ODATE2} ${ODATE3} ${ODATE4}
#             ${ODATE5} ${ODATE6} ${ODATE7} )


#var_id=1
#VAR_LIST=( 'u' 'v' 'temp' 'q' 'rh' 'ps' )
#NOBS_LIST=( 12474 12474 745 332 0 7206 )
#for varname , num_of_obs  in ${VAR_LIST[@]} , ${NOBS_LIST[@]} ; do

#for ODATE in ${ODATE_LIST[@]} ; do
#  var_id=1
#  echo ${ODATE}
#  for varname in ${VAR_LIST[@]} ; do
#  obsdir="/work3/kurihana/script_etc/koji/workspace/work_readprepbufr/bin/${ODATE}"
#  obsinfile="bin_${ODATE}_${varname}.dat"
#  ifile=${obsdir}/$obsinfile
#  fbyte=` wc -c < ${ifile}`
#  if [ $fbyte -gt 0  ]; then echo $var_id ; fi
#  var_id=`expr  ${var_id} + 1 `
#  done
#done
#exit


# nicam_database dir
#nicam_databasedir=/work3/kurihana/NICAM14.4_ATM_version/NICAM/NICAM_SAVEDATA
nicam_databasedir=/work3/kurihana/NICAM16.3/NICAM/NICAM_SAVEDATA

# llmap dir
llmapinfo=${nicam_databasedir}/llmap/gl0${gl}/rl0${rl}/llmap.info
llmapfile=${nicam_databasedir}/llmap/gl0${gl}/rl0${rl}/llmap.rgn

# observation data dir
obs_binary_size=4
obsdir="/work3/kurihana/script_etc/koji/workspace/work_readprepbufr/bin/${yyyy}${mm}${dd}${hh}"
obsinfile="bin_${yyyy}${mm}${dd}${hh}_${varname}.dat"

# make binary icogrid-latlon-database
inputhgridname="/work3/kurihana/NICAM16.3/NICAM/NICAM_SAVEDATA/hgrid/gl0${gl}/rl0${rl}"
outputbasename=hgrid

# self-made icollgird dir
ico2lldir='icohgrid'
ico2ll_hgriddir='./icohgrid'
ico2ll_hgridname=bin_hgrid

# vertical grid info
zgridinfo='./ico_38zlayer.txt'

# nicam horizontal grid info
#nicam_hgrid="${nicam_databasedir}/hgrid_intel/gl0${gl}/rl0${rl}"
nicam_hgrid="${nicam_databasedir}/hgrid/gl0${gl}/rl0${rl}"
nicam_hgridbasename='grid.rgn'

# fnl reanalysis data info
fnl_nlon=360
fnl_nlat=181
fnl_tempfiledir='./'
fnl_tempfile='fnl_tmp_201111010000.dat'

# log file name
logfname="${yyyy}${mm}${dd}${hh}_${varname}"

###############################################################################

mkdir -p ${ico2lldir}

rm -f *.cnf
cat >>hssdriver.cnf<<EOF
  &hssdriver_param
    gl     = ${gl},
    rl     = ${rl},
    zall   = ${zall},
    obsall = ${nobs_prepbufr},
    offline           = .true.  ,
    opt_bin           = .true.  ,
    opt_outlon_0to360 = .false. ,
    opt_zcoords       = .false. ,
    opt_pcoords       = .true.  ,
    obsdatadir    ="${obsdir}"
    obsinputfile  ="${obsinfile}"
    inputdir      ="${inputhgridname}"
    input_size    = ${obs_binary_size}, 
    outputdir     ="${ico2lldir}"
    outputname    ="${outputbasename}"
    rgnfilename   ="${ico2ll_hgridname}",
    inputfiledir  ="${ico2ll_hgriddir}",
    zgridinfo     ="${zgridinfo}",
    logfname      ="${logfname}",
  /

  &latlon_param
    inputinfofile = "${llmapinfo}" ,
    inputgridfile = "${llmapfile}" ,
  /

  &dll_param
    inputdir  ="${ico2ll_hgriddir}",
    basename  ="${ico2ll_hgridname}",
  /

  &icoid_param
    opt_bin        = .true. ,
    use_icohgrid   = .true. ,
    input_size     = 8,
    icodatadir     = "${nicam_hgrid}",
    icobasename    = "${nicam_hgridbasename}"
  /

  &p2z_param
    nlon             = ${fnl_nlon},
    nlat             = ${fnl_nlat},
    fnl_tempfiledir  = "${fnl_tempfiledir}",
    fnl_tempfile     = "${fnl_tempfile}",
  /

  &prepbufr_param
    pb_id     = ${var_id},
  /

EOF

echo ${yyyy}${mm}${dd}${hh} ${var_id} ${varname}

#
# L1 compiler option
# L2 intger rgn name to character
# L3 ico2ll
# L4 rgn-grid selection modules
# L5 g-grid number selection & dist. comp.  modules
#
#rm -f testhss
#ifort -assume byterecl -traceback -g -debug all -warn all  -CB   \
#ifort -assume byterecl -traceback -CB   \
#     mod_int2fname.f90  \
#     mod_ico2latlon.f90 \
#     mod_qs.f90  \
#     mod_llselect.f90 mod_getrgn.f90 \
#     mod_dlatlon.f90 mod_distdeg.f90  mod_gridsearch.f90 \
#     mod_hssll.f90 \
#     mod_hssrgn.f90 \
#     mod_compave.f90 \
#     test_hssobs_driver.f90 \
#     -o testhss

#./testhss
exit
