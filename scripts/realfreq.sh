#!/bin/bash
#================================================================
# HEADER
#================================================================
#% SYNOPSIS
#+    ${SCRIPT_NAME} -m [directory] -f [reference] [options ...]
#%
#% DESCRIPTION
#%    Realtime methylation frequency computation of a given sequencing directory.
#%
#% OPTIONS
#%
#%    -h, --help                                    Print help message
#%    -i, --info                                    Print script information
#%    -m [directory]                                The sequencing experiment directory to be monitored
#%    -f [reference]                                Reference genome for alignment
#%    -o [output]                                   Output file for modification frequency [default: freq.tsv]
#%    -r                                            Resumes a previous live conversion
#%    -c [port]                                     Server port for realfreq
#%    -t [time]                                     Timeout in seconds [default: 21600]
#%    -p [processes]                                Maximum number of parallel conversion processes [default: 1]
#%    -a [extension]                                Watch for files with extension [default: pod5]
#%    -b                                            Output bedmethyl format
#%
#% ADVANCED/DEBUGGING OPTIONS
#%
#%    -n                                            Specify non-realtime analysis
#%    -d [filename]                                 Specify custom location for the list of attempted files [default: monitor_dir/realfreq_attempted_list.log]
#%    -l [filename]                                 Specify custom log filename [default: monitor_dir/realfreq.log]
#%    -f [file]                                     Specify location for the list of files that failed to convert [default: monitor_dir/realfreq_failed_list.log]
#%    -s [file]                                     Specify custom script for handling conversion [default: script_location/pipeline.sh]
#%    -y, --yes                                     Say yes to 'Are you sure?' message in advance for overwriting
#%
#================================================================
#- IMPLEMENTATION
#-    authors         Hasindu GAMAARACHCHI (hasindu@unsw.edu.au),
#-                    Sasha JENNER (jenner.sasha@gmail.com),
#-                    Suneth SAMARASINGHE (imsuneth@gmai.com)
#-    license         MIT
#-
#-    Copyright (c) 2019 Hasindu Gamaarachchi, 2020 Sasha Jenner, 2024 Suneth Samarasinghe
#-
#-    Permission is hereby granted, free of charge, to any person obtaining a copy
#-    of this software and associated documentation files (the "Software"), to deal
#-    in the Software without restriction, including without limitation the rights
#-    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#-    copies of the Software, and to permit persons to whom the Software is
#-    furnished to do so, subject to the following conditions:
#-
#-    The above copyright notice and this permission notice shall be included in all
#-    copies or substantial portions of the Software.
#-
#-    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#-    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#-    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#-    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#-    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#-    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#-    SOFTWARE.
#-
#================================================================
# END_OF_HEADER
#================================================================

    #== Necessary variables ==#
SCRIPT_HEADSIZE=$(head -200 ${0} | grep -n "^# END_OF_HEADER" | cut -f1 -d:)
SCRIPT_NAME="$(basename ${0})"
SCRIPT_PATH="$( cd "$(dirname "$0")" ; pwd -P )" # Scripts current path

    #== Usage functions ==#
usage() { printf "Usage: "; head -${SCRIPT_HEADSIZE:-99} ${0} | grep -e "^#+" | sed -e "s/^#+[ ]*//g" -e "s/\${SCRIPT_NAME}/${SCRIPT_NAME}/g"; }
usagefull() { head -${SCRIPT_HEADSIZE:-99} ${0} | grep -e "^#[%+]" | sed -e "s/^#[%+-]//g" -e "s/\${SCRIPT_NAME}/${SCRIPT_NAME}/g"; }
scriptinfo() { head -${SCRIPT_HEADSIZE:-99} ${0} | grep -e "^#-" | sed -e "s/^#-//g" -e "s/\${SCRIPT_NAME}/${SCRIPT_NAME}/g"; }

    #== Default variables ==#

# Default script to be copied and run on the worker nodes
PIPELINE_SCRIPT="$SCRIPT_PATH/pipeline.sh"

# Set options by default
resuming=false
realtime=true
say_yes=false

# Default timeout of 6 hours
TIME_INACTIVE=21600

# maximum number of parallel coversion processes
MAX_PROC=1

# Assume necessary options not set
monitor_dir_specified=false
MONITOR_PARENT_DIR=
OUTPUT_FILE=
server_port=""
MONITOR_EXTENSION="pod5"
bedmethyl_output=false

## Handle flags
while getopts "m:f:o:l:hinryd:t:s:f:p:c:a:b" o; do
    case "${o}" in
        m)
            MONITOR_PARENT_DIR=${OPTARG}
            monitor_dir_specified=true
            ;;
        f)
            REF=${OPTARG}
            ;;
        o)
            OUTPUT_FILE=${OPTARG}
            ;;
		l)
            LOG=${OPTARG}
			;;
        h)
            usagefull
            exit 0
            ;;
        i)
            scriptinfo
            exit 0
            ;;
        n)
            realtime=false
            ;;
        r)
            resuming=true
            ;;
        y)
            say_yes=true
            ;;
        d)
            TMP_FILE_PATH="${OPTARG}"
            ;;
        t)
            TIME_INACTIVE="${OPTARG}"
            ;;
        s)
            PIPELINE_SCRIPT=${OPTARG}
            ;;
        f)
            FAILED_LIST=${OPTARG}
            ;;
        p)
            MAX_PROC=${OPTARG}
            ;;
        c)
            server_port=${OPTARG}
            ;;
        a)
            MONITOR_EXTENSION=${OPTARG}
            ;;
        b)
            bedmethyl_output=true
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

# Colour codes for printing
YELLOW="\e[33m"
RED="\e[31m"
NORMAL="\033[0;39m"

if [ ! -e $PIPELINE_SCRIPT ]; then
    echo -e $RED"[$SCRIPT_NAME] Pipeline script $PIPELINE_SCRIPT not found. Exiting."$NORMAL
    exit 1
fi

# If either format or monitor option not set
if ! ($monitor_dir_specified); then
    if ! $monitor_dir_specified; then echo "[$SCRIPT_NAME] No monitor directory specified!"; fi
	usage
	exit 1
fi

# If output file not set
if [ -z ${OUTPUT_FILE} ]; then
    OUTPUT_FILE="$MONITOR_PARENT_DIR/freq.tsv"
    echo -e "[$SCRIPT_NAME] Output file not set. Using default $OUTPUT_FILE"
fi

DUMP_FILE="$MONITOR_PARENT_DIR/realfreq_prog_dump.tmp"

server_port_flag=
if [ ! $server_port == "" ]; then
    server_port_flag="-s $server_port"
fi

bedmethyl_output_flag=
if [ $bedmethyl_output == true ]; then
    bedmethyl_output_flag="-b"
fi

resume_flag=
if [ $resuming == true ]; then
    resume_flag="-r"
fi

if [ -z ${REALFREQ_THREADS} ]; then
    REALFREQ_THREADS=1
fi

# set the temporary file and log file
[ -z ${TMP_FILE_PATH} ] && TMP_FILE_PATH=${MONITOR_PARENT_DIR}/realfreq_success_list.log
[ -z ${FAILED_LIST} ] && FAILED_LIST=${MONITOR_PARENT_DIR}/realfreq_fail_list.log
[ -z ${LOG} ] && LOG=${MONITOR_PARENT_DIR}/realfreq_script.log
MONITOR_TRACE=${MONITOR_PARENT_DIR}/realfreq_monitor_trace.log              #trace of the monitor for debugging
MONITOR_TEMP=${MONITOR_PARENT_DIR}/realfreq_monitor_temp                    #used internally to communicate with the monitor
PIPELINE_TIME_LOG=${MONITOR_PARENT_DIR}/realfreq_pipeline_time.log
REALFREQ_PROG_LOG=${MONITOR_PARENT_DIR}/realfreq_prog.log

which inotifywait &> /dev/null || { echo -e $RED"[$SCRIPT_NAME] inotifywait not found! On ubuntu: sudo apt install inotify-tools"$NORMAL; exit 1; }
[ -z ${REALFREQ} ] && REALFREQ=realfreq
${REALFREQ} -V &> /dev/null || { echo -e $RED"[$SCRIPT_NAME] realfreq not found! Add realfreq to PATH or export REALFREQ=/path/to/realfreq"$NORMAL; exit 1;}

# Perform pipeline tool check
"${PIPELINE_SCRIPT}" -c || { echo -e $RED"[$SCRIPT_NAME] Pipeline script tool check failed. Exiting."$NORMAL; exit 1; }

#== Echo the options ==#
echo "[$SCRIPT_NAME] Current options:"
echo -e "\tMonitor directory:\t $MONITOR_PARENT_DIR"
echo -e "\tOutput file:\t\t $OUTPUT_FILE"
echo -e "\tDump file:\t\t $DUMP_FILE"
echo -e "\tIs resuming:\t\t $resuming"
echo -e "\tIs realtime:\t\t $realtime"
echo -e "\tServer port:\t\t $server_port"
echo -e "\tProcessed list:\t\t $TMP_FILE_PATH"
echo -e "\tFailed list:\t\t $FAILED_LIST"
echo -e "\tIdle time:\t\t $TIME_INACTIVE"
echo -e "\tMax pipeline processes:\t $MAX_PROC"
echo -e "\tMonitor watch for:\t $MONITOR_EXTENSION"
echo -e "\tBedmethyl output:\t $bedmethyl_output"
echo -e "\tREALFREQ_THREADS:\t $REALFREQ_THREADS"
echo -e "\tREALFREQ_AUTO:\t\t $REALFREQ_AUTO"
echo -e "\tPipeline script:\t $PIPELINE_SCRIPT"
echo -e "\tMonitor trace:\t\t $MONITOR_TRACE"
echo -e "\tLog file:\t\t $LOG"
echo -e "\tPipeline time log:\t $PIPELINE_TIME_LOG"
echo -e "\tRealfreq log:\t\t $REALFREQ_PROG_LOG"
echo -e "\tSay yes:\t\t $say_yes"

# wait till the monitor directory is available
if [ ! -d $MONITOR_PARENT_DIR ]; then
    echo "[$SCRIPT_NAME] Monitor directory $MONITOR_PARENT_DIR not found. Waiting for it to be created."
    while [ ! -d $MONITOR_PARENT_DIR ]; do
        sleep 1
    done
fi

# Warn before cleaning logs
if ! $resuming && ! $say_yes; then # If not resuming

    if [ -e ${LOG} ]
        then
            while true; do
                read -p "[$SCRIPT_NAME] A previous log file exist at ${LOG}! Are you sure you want to remove them and start over (y/n)" response
                case $response in
                    [Yy]* )
                        test -e $LOG && rm $LOG # Empty log file
                        test -e $MONITOR_TRACE && rm $MONITOR_TRACE # Empty log file
                        test -e $FAILED_LIST && rm $FAILED_LIST # Empty log file
                        test -e $TMP_FILE_PATH && rm $TMP_FILE_PATH # Empty log file
                        test -e $PIPELINE_TIME_LOG && rm $PIPELINE_TIME_LOG # Empty log file
                        test -e $REALFREQ_PROG_LOG && rm $REALFREQ_PROG_LOG # Empty log file
                        break
                        ;;

                    [Nn]* )
                        exit 0
                        ;;

                    * )
                        echo "[$SCRIPT_NAME] Please answer yes or no."
                        ;;
                esac
        done
    fi

fi

# Function to catch the bam file output from various tools in pipeline
catch_bam() {
    while read line; do
        if [[ $line == *"realfreq-pipeline-output"* ]]; then
            modbam_file=$(echo $line | grep -oP '(?<=realfreq-pipeline-output:).*')
            CURRENT_BAM_FILEPATH=$modbam_file
            echo "$modbam_file"
        fi
    done
}

# start realfreq
echo "[$SCRIPT_NAME] Starting realfreq" | tee -a $LOG
PIPE="pipeline_output.pipe"
if [[ -p $PIPE ]]; then
    rm $PIPE
fi
mkfifo $PIPE
${REALFREQ} -t 1 ${bedmethyl_output_flag} ${server_port_flag} ${resume_flag} -d $DUMP_FILE -o $OUTPUT_FILE $REF -l $TMP_FILE_PATH < $PIPE 2>$REALFREQ_PROG_LOG &

if ! $realtime; then # If non-realtime option set
    echo "[$SCRIPT_NAME] Non realtime conversion of all files in $MONITOR_PARENT_DIR" | tee -a $LOG
    test -e $TMP_FILE_PATH && rm $TMP_FILE_PATH

    find $MONITOR_PARENT_DIR/ -name "*.${MONITOR_EXTENSION}" | "$PIPELINE_SCRIPT" -l $PIPELINE_TIME_LOG -f $FAILED_LIST -p $MAX_PROC | catch_bam > $PIPE

else # Else assume realtime analysis is desired

    if $resuming; then # If resuming option set
        echo "[$SCRIPT_NAME] resuming" | tee -a $LOG
        if [ ! -e $DUMP_FILE ]; then
            echo "[$SCRIPT_NAME] Dump file $DUMP_FILE not found. Exiting." | tee -a $LOG
            exit 1
        fi

        "$SCRIPT_PATH"/monitor/monitor.sh -x ${MONITOR_EXTENSION} -t $TIME_INACTIVE -f -d ${MONITOR_TEMP} $MONITOR_PARENT_DIR/  |
        "$SCRIPT_PATH"/monitor/ensure.sh -x ${MONITOR_EXTENSION} -r -d $TMP_FILE_PATH -l ${MONITOR_TRACE}  |
        "$PIPELINE_SCRIPT" -l $PIPELINE_TIME_LOG -f $FAILED_LIST -p $MAX_PROC | catch_bam > $PIPE
        
    else
        echo "[$SCRIPT_NAME] running" | tee -a $LOG
        test -e $TMP_FILE_PATH && rm $TMP_FILE_PATH

        "$SCRIPT_PATH"/monitor/monitor.sh -x ${MONITOR_EXTENSION} -t $TIME_INACTIVE -f -d ${MONITOR_TEMP} $MONITOR_PARENT_DIR/  |
        "$SCRIPT_PATH"/monitor/ensure.sh -x ${MONITOR_EXTENSION} -d $TMP_FILE_PATH -l ${MONITOR_TRACE}  |
        "$PIPELINE_SCRIPT" -l $PIPELINE_TIME_LOG  -f $FAILED_LIST -p $MAX_PROC  | catch_bam > $PIPE
        
    fi
    
    echo "[$SCRIPT_NAME] No new ${MONITOR_EXTENSION} files found in last ${TIME_INACTIVE} seconds." | tee -a $LOG
    echo "[$SCRIPT_NAME] converting left overs" | tee -a $LO

    find $MONITOR_PARENT_DIR/ -name "*.${MONITOR_EXTENSION}"   |
    "$SCRIPT_PATH"/monitor/ensure.sh -x ${MONITOR_EXTENSION} -r -d $TMP_FILE_PATH -l ${MONITOR_TRACE}  |
    "$PIPELINE_SCRIPT" -l $PIPELINE_TIME_LOG  -f $FAILED_LIST -p $MAX_PROC | catch_bam > $PIPE
    
fi

test -e $FAILED_LIST && echo -e $RED"[$SCRIPT_NAME] $(wc -l $FAILED_LIST) files failed the pipeline. See $FAILED_LIST for the list"$NORMAL | tee -a $LOG
NUMPOD5=$(find $MONITOR_PARENT_DIR/ -name "*.${MONITOR_EXTENSION}" | wc -l)
NUMBAM=$(find $MONITOR_PARENT_DIR/ -name "*.${MONITOR_EXTENSION}" | wc -l)
if [ ${NUMPOD5} -ne ${NUMBAM} ] ; then
    echo -e $RED"[$SCRIPT_NAME] In $MONITOR_PARENT_DIR, $NUMPOD5 pod5 files, but only $NUMBAM bam files. Check the logs for any failures."$NORMAL | tee -a $LOG
else
    echo "[$SCRIPT_NAME] In $MONITOR_PARENT_DIR, $NUMPOD5 pod5 files, $NUMBAM blow5 files." | tee -a $LOG
fi

echo "Scanning for errors in log files" | tee -a $LOG
find $MONITOR_PARENT_DIR/ -name '*.log' -exec cat {} \; | grep -i "ERROR" | tee -a $LOG
echo "Scanning for warnings in log files" | tee -a $LOG
find $MONITOR_PARENT_DIR/ -name '*.log' -exec cat {} \; | grep -i "WARNING" | tee -a $LOG
echo "[$SCRIPT_NAME] exiting" | tee -a $LOG

# Cleanup
wait
rm $PIPE