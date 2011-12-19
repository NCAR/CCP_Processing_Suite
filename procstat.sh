#!/bin/sh
#
# ASB - procstat script to monitor CCP_Processing_Suite as it runs
# from the script generated by Process_Setup
# usage: procstat.sh [error | start | complete] [procedure name (for error only)]
# examplet: ./procstat.sh error mss_read
#
_DEBUG="on"
function DEBUG()
{
 [ "$_DEBUG" == "on" ] &&  $@
}

# for debugging, hardcode - will get these from env. vars.
#
#CASE=b40.1850.track1.1deg.006a
#HIST=cam2.h0
#TPER=mon
#HTYP=atm
#LOCATION=NC
#MSSPROC=/CCSM/csm/b40.1850.track1.1deg.006a/atm/proc/tseries/monthly

#
# generate a base filename
#
#FILEBASE=${CASE}_${HIST}_${TPER}_${DO_PROC}_${DO_CMIP} ; export FILEBASE
FILEBASE=${CASE}.${HIST}.${TPER} ; export FILEBASE


#
# generate the filenames to be used
#
FILEPS=${FILEBASE}.ps ; export FILEPS
FILELOG=${FILEBASE}.log ; export FILELOG
FILEDF=${FILEBASE}.df ; export FILEDF
FILEHSI=${FILEBASE}.hsi; export FILEHSI

#
# get current user name and pp hostname
#
USER=`whoami`
HOSTNAME=`hostname`
#
# need to follow a naming convention for the logfiles
# execute and save to files
#
case "${LOCATION}" in
  NC )                                        # NCAR experiments
    LOGFILE=`ls ~/CCP_Processing_Suite/log*${CASE}*${HIST}*`
    date > ${FILEPS}
    ps auwx | grep ${USER} >> ${FILEPS}
    date > ${FILELOG}
    tail -n 40 ${LOGFILE} >> ${FILELOG} ;;
  NE )                                        # NERSC experiments
    LOGFILE=`ls ~/CCP_Processing_Suite/log*${CASE}*${HIST}*` 
    date > ${FILEPS}
    ps auwx | grep ${USER} >> ${FILEPS}
    date > ${FILELOG}
    tail -n 40 ${LOGFILE} >> ${FILELOG} ;;
  OR )                                        # ORNL experiments
    LOGFILE=`ls ~/CCP_Processing_Suite/*.ER`
    date > ${FILEPS}
    qstat >> ${FILEPS}
    date > ${FILELOG}
    tail -n 40 ${LOGFILE} >> ${FILELOG} ;;
  * )
    date > ${FILELOG}
    echo "procstat : LOGFILE does not exist for "${LOCATION}" - "${FILEBASE} >> ${FILELOG}
    ;;
esac

date > ${FILEDF}
df -h >> ${FILEDF}

#
# get the mass store listing to see where we are in the processing using $MSSPROC
#
#
# If NCAR MSS msrcp command exists, use it.
date > ${FILEHSI}

#
# all sites should be using hsi
TEST4HSI=`which hsi 2<&1`
if [ $? -eq 0 ] ; then
    rm -f tossme
    hsi -q "cd $MSSPROC; ls -Fl" >& tossme
    if [ $? -eq 0 ] ; then
      egrep "\.${HIST}\." tossme | cut -c60- >> $FILEHSI
      rm -f tossme
    else
      echo "procstat : (files may not have been written yet) Error on hsi from "$MSSPROC >> $FILEHSI
      rm -f tossme
    fi
else
    echo "procstat.sh: cannot find hsi command." >> $FILEHSI
fi

#
# create the mail messages
#
MAILTO=procstat@cgd.ucar.edu

mail -s "procstat ${FILEPS}" ${MAILTO} < ${FILEPS}
mail -s "procstat ${FILELOG}" ${MAILTO} < ${FILELOG}
mail -s "procstat ${FILEDF}" ${MAILTO} < ${FILEDF}
mail -s "procstat ${FILEHSI}" ${MAILTO} < ${FILEHSI}

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
    'aliceb' | 'ilana' | 'strandwg' | 'hteng' | 'jma' )
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
