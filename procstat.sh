#!/bin/sh
#
# ASB - procstat script to monitor CCP_Processing_Suite as it runs
# from the script generated by Process_Setup or Process_CMIP5
# usage: procstat.sh [error | start | complete] [procedure name (for error only)]
# example: ./procstat.sh error mss_read
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
FILEPS=${FILEBASE}.ps ; export FILEPS
FILELOG=${FILEBASE}.log ; export FILELOG
FILEDF=${FILEBASE}.df ; export FILEDF
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
    date > ${FILEPS}
    ps auwx | grep ${USER} >> ${FILEPS}
    date > ${FILELOG}
    ;;
  lens* )                                  # lens @ ORNL
    date > ${FILEPS}
    qstat >> ${FILEPS}
    date > ${FILELOG}
    ;;
  * )
    date > ${FILELOG}
    echo "procstat: "${HOSTNAME}" unknown - "${FILEBASE} >> ${FILELOG}
    exit 1
    ;;
esac
#
date > ${FILEDF}
df -h >> ${FILEDF}
#
# create the mail messages
#
MAILTO=procstat@cgd.ucar.edu
mail -s "procstat ${FILEPS}" ${MAILTO} < ${FILEPS}
mail -s "procstat ${FILELOG}" ${MAILTO} < ${FILELOG}
mail -s "procstat ${FILEDF}" ${MAILTO} < ${FILEDF}
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
    'aliceb' | 'ilana' | 'strandwg' | 'hteng' | 'jma' | 'asphilli' )
         MAILTOUSER=${USER}@ucar.edu ;;
    'abertini')
         MAILTOUSER=aliceb@ucar.edu  ;;
    'jmarb')
         MAILTOUSER=jma@ucar.edu     ;;
    'wgstrand')
         MAILTOUSER=strandwg@ucar.edu ;;
    esac

    rm -f tossme
    
    if test $1 = "start" 
	then
	    FILESTART=start.${FILEBASE} ; export FILESTART
	    echo "STATUS = START" > tossme
	    echo "CASE = ${CASE}" >> tossme
	    echo "HIST = ${HIST}" >> tossme
	    echo "TPER = ${TPER}" >> tossme
	    echo "LOCATION = ${LOCATION}" >> tossme
	    echo "HTYP = ${HTYP}" >> tossme
	    echo "PP_HOST = ${HOSTNAME}" >> tossme
	    echo "USER = ${USER}" >> tossme
	    echo "NOTE = This email is being sent to you automatically from the procstat.sh program. For more details, go to the procstat web at http://webint.cgd.ucar.edu/project/ccr/procstat/cgi-bin/index.cgi." >> tossme
	    mail -s "procstat ${FILESTART}" ${MAILTO} -c ${MAILTOUSER} < tossme
	    rm -f tossme
    fi
    if test $1 = "error" 
        then
	    FILEERROR=error.${FILEBASE} ; export FILEERROR
	    echo "STATUS = ERROR" > tossme
	    echo "CASE = ${CASE}" >> tossme
	    echo "HIST = ${HIST}" >> tossme
	    echo "TPER = ${TPER}" >> tossme
	    echo "LOCATION = ${LOCATION}" >> tossme
	    echo "HTYP = ${HTYP}" >> tossme
	    echo "PP_HOST = ${HOSTNAME}" >> tossme
	    echo "USER = ${USER}" >> tossme
	    echo "CALLING_PROGRAM = $2" >> tossme
	    echo "NOTE = This email is being sent to you automatically from the procstat.sh program. For more details, go to the procstat web at http://webint.cgd.ucar.edu/project/ccr/procstat/cgi-bin/index.cgi." >> tossme
	    mail -s "procstat ${FILEERROR}" ${MAILTO} -c ${MAILTOUSER} < tossme
	    rm -f tossme
    fi
    if test $1 = "complete" 
        then
	    FILECOMPLETE=complete.${FILEBASE} ; export FILECOMPLETE
	    echo "STATUS = COMPLETE" > tossme
	    echo "CASE = ${CASE}" >> tossme
	    echo "HIST = ${HIST}" >> tossme
	    echo "TPER = ${TPER}" >> tossme
	    echo "LOCATION = ${LOCATION}" >> tossme
	    echo "HTYP = ${HTYP}" >> tossme
	    echo "PP_HOST = ${HOSTNAME}" >> tossme
	    echo "USER = ${USER}" >> tossme
	    echo "NOTE = This email is being sent to you automatically from the procstat.sh program. For more details, go to the procstat web at http://webint.cgd.ucar.edu/project/ccr/procstat/cgi-bin/index.cgi. Please double check log file and mass storage listing for any errors and then manually change the status to COMPLETE." >> tossme
	    mail -s "procstat ${FILECOMPLETE}" ${MAILTO} -c ${MAILTOUSER} < tossme
	    rm -f tossme
    fi
fi

_DEBUG="off"

#echo "procstat = finish"
