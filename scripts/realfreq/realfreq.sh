#!/bin/bash
#================================================================
# HEADER
#================================================================
#% SYNOPSIS
#+    ${SCRIPT_NAME} -m [directory] -g [guppy_bin] -f [reference] -x [reference_index] -e [model] [options ...]
#%
#% DESCRIPTION
#%    Runs realtime POD5 to SLOW5 conversion of for a given sequencing directory.
#%
#% OPTIONS
#%
#%    -h, --help                                    Print help message
#%    -i, --info                                    Print script information
#%    -m [directory]                                The sequencing experiment directory to be monitored
#%    -g [guppy_bin]                                Path to guppy binary
#%    -f [reference]                                Reference genome for alignment
#%    -x [reference_index]                          Reference genome index for alignment
#%    -e [model]                                    Model for guppy basecalling
#%    -o [output]                                   Output file for modification frequency [default: freq.tsv]
#%    -r                                            Resumes a previous live conversion
#%    -t [time]                                     Timeout in seconds [default: 21600]
#%    -p [processes]                                Maximum number of parallel conversion processes [default: 1]
#%
#% ADVANCED/DEBUGGING OPTIONS
#%
#%    -n                                            Specify non-realtime analysis
#%    -d [filename]                                 Specify custom location for the list of attempted files [default: monitor_dir/realtime_p2s_attempted_list.log]
#%    -l [filename]                                 Specify custom log filename [default: monitor_dir/realtime_p2s.log]
#%    -f [file]                                     Specify location for the list of files that failed to convert [default: monitor_dir/realtime_p2s_failed_list.log]
#%    -s [file]                                     Specify custom script for handling conversion [default: script_location/pipeline.sh]
#%    -y, --yes                                     Say yes to 'Are you sure?' message in advance for overwriting
#%
#% EXAMPLES
#%    convert
#%        ${SCRIPT_NAME} -m [directory]
#%    resume convert
#%        ${SCRIPT_NAME} -m [directory] -r
#%    one hour timeout
#%        ${SCRIPT_NAME} -m [directory] -t 3600
#%
#================================================================
#- IMPLEMENTATION
#-    authors         Hasindu GAMAARACHCHI (hasindu@unsw.edu.au),
#-                    Sasha JENNER (jenner.sasha@gmail.com),
#-                    Suneth SAMARASINGHE (suneth@unsw.edu.au)
#-    license         MIT
#-
#-    Copyright (c) 2019 Hasindu Gamaarachchi, 2020 Sasha Jenner
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
GUPPY_BIN=
REF=
REFIDX=
MODEL=
OUTPUT_FILE=

## Handle flags
while getopts "ihnyrm:g:f:x:e:o:l:t:d:f:s:p:" o; do
    case "${o}" in
        m)
            MONITOR_PARENT_DIR=${OPTARG}
            monitor_dir_specified=true
            ;;
        g)
            GUPPY_BIN=${OPTARG}
            ;;
        f)
            REF=${OPTARG}
            ;;
        x)
            REFIDX=${OPTARG}
            ;;
        e)
            MODEL=${OPTARG}
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

# If either format or monitor option not set
if ! ($monitor_dir_specified); then
    if ! $monitor_dir_specified; then echo "[realp2s.sh] No monitor directory specified!"; fi
	usage
	exit 1
fi

# If GUPPY_BIN not set
if [ -z ${GUPPY_BIN} ]; then
    echo -e $RED"[realp2s.sh] Guppy binary not set! Set with -g option"$NORMAL
    exit 1
fi

# If reference genome not set
if [ -z ${REF} ]; then
    echo -e $RED"[realp2s.sh] Reference genome not set! Set with -f option"$NORMAL
    exit 1
fi

# If reference genome index not set
if [ -z ${REFIDX} ]; then
    echo -e $RED"[realp2s.sh] Reference genome index not set! Set with -x option"$NORMAL
    exit 1
fi

# If model not set
if [ -z ${MODEL} ]; then
    echo -e $RED"[realp2s.sh] Model not set! Set with -e option"$NORMAL
    exit 1
fi

# If output file not set
if [ -z ${OUTPUT_FILE} ]; then
    OUTPUT_FILE="$MONITOR_PARENT_DIR/freq.tsv"
    echo -e "[realp2s.sh] Output file not set. Using default $MONITOR_PARENT_DIR"
fi

# wait till the monitor directory is available
if [ ! -d $MONITOR_PARENT_DIR ]; then
    echo "[realp2s.sh] $MONITOR_PARENT_DIR not found. Waiting for it to be created."
    while [ ! -d $MONITOR_PARENT_DIR ]; do
        sleep 1
    done
fi

DUMP_FILE="$MONITOR_PARENT_DIR/dump.tmp"

# set the temporary file and log file
[ -z ${TMP_FILE_PATH} ] && TMP_FILE_PATH=${MONITOR_PARENT_DIR}/realtime_p2s_attempted_list.log
echo "[realp2s.sh] Temporary file location that saves the state $TMP_FILE_PATH"
[ -z ${FAILED_LIST} ] && FAILED_LIST=${MONITOR_PARENT_DIR}/realtime_p2s_failed_list.log
echo "[realp2s.sh] Any pod5 files that failed conversions will be written to $FAILED_LIST"
[ -z ${LOG} ] && LOG=${MONITOR_PARENT_DIR}/realtime_p2s.log
echo "[realp2s.sh] Master log file location ${LOG}"
MONITOR_TRACE=${MONITOR_PARENT_DIR}/realtime_p2s_monitor_trace.log              #trace of the monitor for debugging
echo "[realp2s.sh] Monitor trace log ${MONITOR_TRACE}"
START_END_TRACE=${MONITOR_PARENT_DIR}/realtime_p2s_start_end_trace.log          #trace for debugging
echo "[realp2s.sh] Start end trace log ${START_END_TRACE}"
MONITOR_TEMP=${MONITOR_PARENT_DIR}/realtime_p2s_monitor_temp                    #used internally to communicate with the monitor
echo "[realp2s.sh] Idle time with no pod5 files to end the program ${TIME_INACTIVE} seconds"
test -d ${MONITOR_PARENT_DIR} || { echo "[realp2s.sh] Monitor directory does not exist!"; exit 1; }
[ -z ${SLOW5TOOLS} ] && export SLOW5TOOLS=slow5tools
[ -z ${BLUECRAB} ] && export BLUECRAB=blue-crab
[ -z ${BUTTERY_EEL} ] && export BUTTERY_EEL=buttery-eel
[ -z ${MINIMAP2} ] && export MINIMAP2=minimap2
[ -z ${SAMTOOLS} ] && export SAMTOOLS=samtools


${SLOW5TOOLS} --version &> /dev/null || { echo -e $RED"[realp2s.sh] slow5tools not found! Either put slow5tools under path or set SLOW5TOOLS variable, e.g.,export SLOW5TOOLS=/path/to/slow5tools"$NORMAL; exit 1;}
${BLUECRAB} --version &> /dev/null || { echo -e $RED"[realp2s.sh] blue-crab not found! Either put blue-crab under path or set BLUECRAB variable, e.g.,export BLUECRAB=/path/to/blue-crab"$NORMAL; exit 1;}
command -v ${BUTTERY_EEL} &> /dev/null || { echo -e $RED"[realp2s.sh] buttery-eel not found! Either put buttery-eel under path or set BUTTERY_EEL variable, e.g.,export BUTTERY_EEL=/path/to/buttery-eel"$NORMAL; exit 1;}
${MINIMAP2} --version &> /dev/null || { echo -e $RED"[realp2s.sh] minimap2 not found! Either put minimap2 under path or set MINIMAP2 variable, e.g.,export MINIMAP2=/path/to/minimap2"$NORMAL; exit 1;}
${SAMTOOLS} --version &> /dev/null || { echo -e $RED"[realp2s.sh] samtools not found! Either put samtools under path or set SAMTOOLS variable, e.g.,export SAMTOOLS=/path/to/samtools"$NORMAL; exit 1;}
which inotifywait &> /dev/null || { echo -e $RED"[realp2s.sh] inotifywait not found! On ubuntu: sudo apt install inotify-tools"$NORMAL; exit 1; }

#== Begin Run ==#

# Warn before cleaning logs
if ! $resuming && ! $say_yes; then # If not resuming

    if [ -e ${LOG} ]
        then
            while true; do
                read -p "[realp2s.sh] A previous log file exist at ${LOG}! Are you sure you want to remove them and start over (y/n)" response
                case $response in
                    [Yy]* )
                        test -e $LOG && rm $LOG # Empty log file
                        test -e $MONITOR_TRACE && rm $MONITOR_TRACE # Empty log file
                        test -e $START_END_TRACE && rm $START_END_TRACE # Empty log file
                        test -e $FAILED_LIST && rm $FAILED_LIST # Empty log file
                        break
                        ;;

                    [Nn]* )
                        exit 0
                        ;;

                    * )
                        echo "[realp2s.sh] Please answer yes or no."
                        ;;
                esac
        done
    fi

fi

# Create folders to copy the results (slow5 files logs)
# test -d $MONITOR_PARENT_DIR/slow5         || mkdir $MONITOR_PARENT_DIR/slow5            || exit 1
# echo "[realp2s.sh] SLOW5 files will be written to $MONITOR_PARENT_DIR/slow5"
# test -d $MONITOR_PARENT_DIR/slow5_logs    || mkdir $MONITOR_PARENT_DIR/slow5_logs       || exit 1
# echo "[realp2s.sh] SLOW5 p2s individual logs will be written to $MONITOR_PARENT_DIR/slow5_logs"


# # Start realfreq
# echo "[realp2s.sh] Starting realfreq" | tee $LOG
# (
# while read line; do
#     if [[ $line == "Finished pipeline"* ]]; then
#         pod5_file=$(echo $line | grep -oP '(?<=for ).*')
#         file_name=$(basename $pod5_file)
#         modbam_file="$pod5_file".remora.bam
#         echo "$modbam_file" | realfreq -r $REF -o freq.tsv |&
#         tee -a $LOG
#     fi
# done < $LOG
# ) &

catch_bam() {
    while read line; do
        if [[ $line == *"Finished pipeline"* ]]; then
            modbam_file=$(echo $line | grep -oP '(?<=modbam: ).*')
            echo "$modbam_file"
        fi
    done
}


if ! $realtime; then # If non-realtime option set
    echo "[realp2s.sh] Non realtime conversion of all files in $MONITOR_PARENT_DIR" | tee $LOG
    test -e $TMP_FILE_PATH && rm $TMP_FILE_PATH
    find $MONITOR_PARENT_DIR/ -name "*.pod5" | "$PIPELINE_SCRIPT" -g $GUPPY_BIN -r $REF -i $REFIDX -m $MODEL | tee $LOG | catch_bam | realfreq -d $DUMP_FILE -r $REF -o $OUTPUT_FILE |&
    tee $LOG

else # Else assume realtime analysis is desired

    # Monitor the new file creation in pod5 folder and execute realtime f5-pipeline script
    # Close after timeout met
    if $resuming; then # If resuming option set
        echo "[realp2s.sh] resuming" | tee -a $LOG
        "$SCRIPT_PATH"/monitor/monitor.sh -t $TIME_INACTIVE -f -d ${MONITOR_TEMP} $MONITOR_PARENT_DIR/  |
        "$SCRIPT_PATH"/monitor/ensure.sh -r -d $TMP_FILE_PATH -l ${MONITOR_TRACE}  |
        "$PIPELINE_SCRIPT" -d $TMP_FILE_PATH -l $START_END_TRACE -p $MAX_PROC -g $GUPPY_BIN -r $REF -i $REFIDX -m $MODEL | tee $LOG | catch_bam | realfreq -s -d $DUMP_FILE -r $REF -o $OUTPUT_FILE |&
        tee -a $LOG | { read file; echo "[realfreq.sh] $file"; } >> $TMP_FILE_PATH
    else
        echo "[realp2s.sh] running" | tee $LOG
        test -e $TMP_FILE_PATH && rm $TMP_FILE_PATH
        "$SCRIPT_PATH"/monitor/monitor.sh -t $TIME_INACTIVE -f -d ${MONITOR_TEMP} $MONITOR_PARENT_DIR/  |
        "$SCRIPT_PATH"/monitor/ensure.sh -d $TMP_FILE_PATH -l ${MONITOR_TRACE}  |
        "$PIPELINE_SCRIPT" -d $TMP_FILE_PATH -l $START_END_TRACE -f $FAILED_LIST -p $MAX_PROC -g $GUPPY_BIN -r $REF -i $REFIDX -m $MODEL | tee $LOG | catch_bam | realfreq -d $DUMP_FILE -r $REF -o $OUTPUT_FILE |&
        tee -a $LOG | { read file; echo "[realfreq.sh] $file"; } >> $TMP_FILE_PATH
    fi
    echo "[realp2s.sh] No new pod5 files found in last ${TIME_INACTIVE} seconds." | tee -a $LOG
    echo "[realp2s.sh] converting left overs" | tee -a $LOG
    find $MONITOR_PARENT_DIR/ -name "*.pod5"   |
    "$SCRIPT_PATH"/monitor/ensure.sh -r -d $TMP_FILE_PATH -l ${MONITOR_TRACE}  |
    "$PIPELINE_SCRIPT" -d $TMP_FILE_PATH -l $START_END_TRACE -f $FAILED_LIST -p $MAX_PROC -g $GUPPY_BIN -r $REF -i $REFIDX -m $MODEL | tee $LOG | catch_bam | realfreq -s -d $DUMP_FILE -r $REF -o $OUTPUT_FILE |&
    tee -a $LOG | { read file; echo "[realfreq.sh] $file"; } >> $TMP_FILE_PATH

fi

test -e $FAILED_LIST && echo -e $RED"[realp2s.sh] $(wc -l $FAILED_LIST) pod5 files failed to convert. See $FAILED_LIST for the list"$NORMAL | tee -a $LOG
NUMPOD5=$(find $MONITOR_PARENT_DIR/ -name '*.pod5' | wc -l)
NUMBLOW5=$(find $MONITOR_PARENT_DIR/ -name '*.blow5' | wc -l)
if [ ${NUMPOD5} -ne ${NUMBLOW5} ] ; then
    echo -e $RED"[realp2s.sh] In $MONITOR_PARENT_DIR, $NUMPOD5 pod5 files, but only $NUMBLOW5 blow5 files. Check the logs for any failures."$NORMAL | tee -a $LOG
else
    echo "[realp2s.sh] In $MONITOR_PARENT_DIR, $NUMPOD5 pod5 files, $NUMBLOW5 blow5 files." | tee -a $LOG
fi
POD5_SIZE=$(find $MONITOR_PARENT_DIR/ -name '*.pod5' -printf "%s\t%p\n" | awk 'BEGIN{sum=0}{sum=sum+$1}END{print sum/(1024*1024*1024)}')
BLOW5_SIZE=$(find $MONITOR_PARENT_DIR/ -name '*.blow5' -printf "%s\t%p\n" | awk 'BEGIN{sum=0}{sum=sum+$1}END{print sum/(1024*1024*1024)}')
SAVINGS=$(echo $POD5_SIZE - $BLOW5_SIZE | bc)
SAVINGS_PERCENT=$(echo "scale=2; $SAVINGS/$POD5_SIZE*100" | bc)
echo "POD5 size: $POD5_SIZE GB" | tee -a $LOG
echo "BLOW5 size: $BLOW5_SIZE GB" | tee -a $LOG
echo "Savings: $SAVINGS GB ($SAVINGS_PERCENT%)" | tee -a $LOG

echo "Scanning for errors in log files" | tee -a $LOG
find $MONITOR_PARENT_DIR/ -name '*.log' -exec cat {} \; | grep -i "ERROR" | tee -a $LOG
echo "Scanning for warnings in log files" | tee -a $LOG
find $MONITOR_PARENT_DIR/ -name '*.log' -exec cat {} \; | grep -i "WARNING" | tee -a $LOG
echo "[realp2s.sh] exiting" | tee -a $LOG
