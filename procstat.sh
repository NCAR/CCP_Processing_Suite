#!/bin/sh
#
# ASB - procstat script to monitor CCP_Processing_Suite as it runs
# from the script generated by Process_Setup or Process_CMIP5
# usage: procstat.sh [error | start | complete] [procedure name (for error only)]
# example: ./procstat.sh error mss_read
#
# ASB - 3/5/2012 - updated to comment out calls to log df -h, ps auwx, and hsi as the time delay
# makes this information not particularly useful
#
_DEBUG="on"
function DEBUG()
{
 [ "$_DEBUG" == "on" ] &&  $@
}
#
# for CMIP5 processing TPER not set but TABL is so use that for TPER
#
if ! [ $TPER ] ; then
   TPER=$TABL
fi
#
# generate a base filename
#
FILEBASE=${CASE}.${HIST}.${TPER} ; export FILEBASE
#
# generate the filenames to be used
#
FILELOG=${FILEBASE}.log ; export FILELOG
#
# get current user name and pp hostname
#
USER=`whoami`
HOSTNAME=`hostname`
#
# need to follow a naming convention for the logfiles
# execute and save to files
#
case "$HOSTNAME" in 
  silver* | tramhill* | hurricane* | mirage* | euclid* )  # NCAR and NERSC machines
    LOGFILE=$HOME/CCP_Processing_Suite/log.${CASE}_${HIST}_${TPER}_process.sh
    if [ -f $LOGFILE ] ; then
      date > ${FILELOG}
      tail -n 40 ${LOGFILE} >> ${FILELOG}
    else
      touch $LOGFILE
      date > ${FILELOG}
      tail -n 40 ${LOGFILE} >> ${FILELOG}
    fi ;;
  lens* )                                  # lens @ ORNL
    LOGFILE=`ls ~/CCP_Processing_Suite/*.ER`
    date > ${FILELOG}
    tail -n 40 ${LOGFILE} >> ${FILELOG} ;;
  * )
    date > ${FILELOG}
    echo "procstat : LOGFILE does not exist for "${LOCATION}" - "${FILEBASE} >> ${FILELOG}
    ;;
esac
#
# create the mail messages
#
MAILTO=procstat@cgd.ucar.edu
mail -s "procstat ${FILELOG}" ${MAILTO} < ${FILELOG}
#
# check if procstat.sh was called with an argument either start, error or complete
# if so, then send an email to trigger an automatic update of the
# database status for this run
#
if [ $# -ge 1 ] ; then
#
# get the user name and corresponding email address to notify user
#
  case "$USER" in
  'aliceb' | 'ilana' | 'strandwg' | 'hteng' | 'jma' | 'asphilli' | 'tilmes' | 'higginsm' )
       MAILTOUSER=${USER}@ucar.edu ;;
  'abertini')
       MAILTOUSER=aliceb@ucar.edu  ;;
  'jmarb')
       MAILTOUSER=jma@ucar.edu     ;;
  'wgstrand')
       MAILTOUSER=strandwg@ucar.edu ;;
  esac
#
  rm -f tossme
  if test $1 = "start" then ;  
    FILESTART=start.${FILEBASE} ; export FILESTART
    echo "STATUS = START" >> tossme
    echo "CASE = ${CASE}" >> tossme
    echo "HIST = ${HIST}" >> tossme
    echo "TPER = ${TPER}" >> tossme
    echo "LOCATION = ${LOCATION}" >> tossme
    echo "HTYP = ${HTYP}" >> tossme
    echo "PP_HOST = ${HOSTNAME}" >> tossme
    echo "USER = ${USER}" >> tossme
    echo "NOTE = This email is being sent to you automatically from the procstat.sh program. For more details, go to the procstat web at http://webint.cgd.ucar.edu/project/ccr/procstat/cgi-bin/index.cgi." >> tossme
    mail -s "procstat ${FILESTART}" ${MAILTO},${MAILTOUSER} < tossme
    rm -f tossme
  fi
  if test $1 = "error" then ; 
    FILEERROR=error.${FILEBASE} ; export FILEERROR
    echo "STATUS = ERROR" >> tossme
    echo "CASE = ${CASE}" >> tossme
    echo "HIST = ${HIST}" >> tossme
    echo "TPER = ${TPER}" >> tossme
    echo "LOCATION = ${LOCATION}" >> tossme
    echo "HTYP = ${HTYP}" >> tossme
    echo "PP_HOST = ${HOSTNAME}" >> tossme
    echo "USER = ${USER}" >> tossme
    echo "CALLING_PROGRAM = $2" >> tossme
    echo "NOTE = This email is being sent to you automatically from the procstat.sh program. For more details, go to the procstat web at http://webint.cgd.ucar.edu/project/ccr/procstat/cgi-bin/index.cgi." >> tossme
    mail -s "procstat ${FILEERROR}" ${MAILTO},${MAILTOUSER} < tossme
    rm -f tossme
  fi
  if test $1 = "complete" then ;
    FILECOMPLETE=complete.${FILEBASE} ; export FILECOMPLETE
    echo "STATUS = COMPLETE" >> tossme
    echo "CASE = ${CASE}" >> tossme
    echo "HIST = ${HIST}" >> tossme
    echo "TPER = ${TPER}" >> tossme
    echo "LOCATION = ${LOCATION}" >> tossme
    echo "HTYP = ${HTYP}" >> tossme
    echo "PP_HOST = ${HOSTNAME}" >> tossme
    echo "USER = ${USER}" >> tossme
    echo "NOTE = This email is being sent to you automatically from the procstat.sh program. For more details, go to the procstat web at http://webint.cgd.ucar.edu/project/ccr/procstat/cgi-bin/index.cgi. Please double check log file and mass storage listing for any errors and then manually change the status to COMPLETE." >> tossme
    mail -s "procstat ${FILECOMPLETE}" ${MAILTO},${MAILTOUSER} < tossme
    rm -f tossme
  fi
fi
_DEBUG="off"
#echo "procstat = finish"
