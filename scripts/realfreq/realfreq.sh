#!/bin/bash
#================================================================
# HEADER
#================================================================
#% SYNOPSIS
#+    ${SCRIPT_NAME} -m [directory] -g [guppy_bin] -f [reference] -x [reference_index] -e [model] [options ...]
#%
#% DESCRIPTION
#%    Realtime methylation frequency computation of a given sequencing directory.
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
#%    -c [port]                                     Server port for realfreq
#%    -t [time]                                     Timeout in seconds [default: 21600]
#%    -p [processes]                                Maximum number of parallel conversion processes [default: 1]
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
server_port=""
THREADS=

## Handle flags
while getopts "ihnyrm:g:f:x:e:o:l:t:d:f:s:p:c:k:" o; do
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
        c)
            server_port=${OPTARG}
            ;;
        k)
            THREADS=${OPTARG}
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
    if ! $monitor_dir_specified; then echo "[$SCRIPT_NAME] No monitor directory specified!"; fi
	usage
	exit 1
fi

# If GUPPY_BIN not set
if [ -z ${GUPPY_BIN} ]; then
    echo -e $RED"[$SCRIPT_NAME] Guppy binary not set! Set with -g option"$NORMAL
    exit 1
fi

# If reference genome not set
if [ -z ${REF} ]; then
    echo -e $RED"[$SCRIPT_NAME] Reference genome not set! Set with -f option"$NORMAL
    exit 1
fi

# If reference genome index not set
if [ -z ${REFIDX} ]; then
    echo -e $RED"[$SCRIPT_NAME] Reference genome index not set! Set with -x option"$NORMAL
    exit 1
fi

# If model not set
if [ -z ${MODEL} ]; then
    echo -e $RED"[$SCRIPT_NAME] Model not set! Set with -e option"$NORMAL
    exit 1
fi

# If output file not set
if [ -z ${OUTPUT_FILE} ]; then
    OUTPUT_FILE="$MONITOR_PARENT_DIR/freq.tsv"
    echo -e "[$SCRIPT_NAME] Output file not set. Using default $OUTPUT_FILE"
fi

# if threads not set
if [ -z ${THREADS} ]; then
    THREADS=8
    echo -e "[$SCRIPT_NAME] Threads not set. Using default $THREADS"
fi

server_port_flag=
if [ ! $server_port == "" ]; then
    server_port_flag="-s $server_port"
fi

# wait till the monitor directory is available
if [ ! -d $MONITOR_PARENT_DIR ]; then
    echo "[$SCRIPT_NAME] $MONITOR_PARENT_DIR not found. Waiting for it to be created."
    while [ ! -d $MONITOR_PARENT_DIR ]; do
        sleep 1
    done
fi

DUMP_FILE="$MONITOR_PARENT_DIR/dump.tmp"

# set the temporary file and log file
[ -z ${TMP_FILE_PATH} ] && TMP_FILE_PATH=${MONITOR_PARENT_DIR}/realfreq_attempted_list.log
[ -z ${FAILED_LIST} ] && FAILED_LIST=${MONITOR_PARENT_DIR}/realfreq_failed_list.log
[ -z ${LOG} ] && LOG=${MONITOR_PARENT_DIR}/realfreq_script.log
MONITOR_TRACE=${MONITOR_PARENT_DIR}/realfreq_monitor_trace.log              #trace of the monitor for debugging
START_END_TRACE=${MONITOR_PARENT_DIR}/realfreq_start_end_trace.log          #trace for debugging
MONITOR_TEMP=${MONITOR_PARENT_DIR}/realfreq_monitor_temp                    #used internally to communicate with the monitor
test -d ${MONITOR_PARENT_DIR} || { echo "[$SCRIPT_NAME] Monitor directory does not exist!"; exit 1; }
[ -z ${SLOW5TOOLS} ] && export SLOW5TOOLS=slow5tools
[ -z ${BLUECRAB} ] && export BLUECRAB=blue-crab
[ -z ${BUTTERY_EEL} ] && export BUTTERY_EEL=buttery-eel
[ -z ${MINIMAP2} ] && export MINIMAP2=minimap2
[ -z ${SAMTOOLS} ] && export SAMTOOLS=samtools


${SLOW5TOOLS} --version &> /dev/null || { echo -e $RED"[$SCRIPT_NAME] slow5tools not found! Either put slow5tools under path or set SLOW5TOOLS variable, e.g.,export SLOW5TOOLS=/path/to/slow5tools"$NORMAL; exit 1;}
${BLUECRAB} --version &> /dev/null || { echo -e $RED"[$SCRIPT_NAME] blue-crab not found! Either put blue-crab under path or set BLUECRAB variable, e.g.,export BLUECRAB=/path/to/blue-crab"$NORMAL; exit 1;}
command -v ${BUTTERY_EEL} &> /dev/null || { echo -e $RED"[$SCRIPT_NAME] buttery-eel not found! Either put buttery-eel under path or set BUTTERY_EEL variable, e.g.,export BUTTERY_EEL=/path/to/buttery-eel"$NORMAL; exit 1;}
${MINIMAP2} --version &> /dev/null || { echo -e $RED"[$SCRIPT_NAME] minimap2 not found! Either put minimap2 under path or set MINIMAP2 variable, e.g.,export MINIMAP2=/path/to/minimap2"$NORMAL; exit 1;}
${SAMTOOLS} --version &> /dev/null || { echo -e $RED"[$SCRIPT_NAME] samtools not found! Either put samtools under path or set SAMTOOLS variable, e.g.,export SAMTOOLS=/path/to/samtools"$NORMAL; exit 1;}
which inotifywait &> /dev/null || { echo -e $RED"[$SCRIPT_NAME] inotifywait not found! On ubuntu: sudo apt install inotify-tools"$NORMAL; exit 1; }

#== Echo the options ==#
echo "[$SCRIPT_NAME] Current options:"
echo -e "[$SCRIPT_NAME]\tMonitor directory:\t $MONITOR_PARENT_DIR"
echo -e "[$SCRIPT_NAME]\tGuppy binary:\t\t $GUPPY_BIN"
echo -e "[$SCRIPT_NAME]\tReference genome:\t $REF"
echo -e "[$SCRIPT_NAME]\tReference genome index:\t $REFIDX"
echo -e "[$SCRIPT_NAME]\tModel:\t\t\t $MODEL"
echo -e "[$SCRIPT_NAME]\tOutput file:\t\t $OUTPUT_FILE"
echo -e "[$SCRIPT_NAME]\tDump file:\t\t $DUMP_FILE"
echo -e "[$SCRIPT_NAME]\tResuming:\t\t $resuming"
echo -e "[$SCRIPT_NAME]\tRealtime:\t\t $realtime"
echo -e "[$SCRIPT_NAME]\tServer port:\t\t $server_port"
echo -e "[$SCRIPT_NAME]\tProcessed list:\t\t $TMP_FILE_PATH"
echo -e "[$SCRIPT_NAME]\tFailed list:\t\t $FAILED_LIST"
echo -e "[$SCRIPT_NAME]\tLog file:\t\t $LOG"
echo -e "[$SCRIPT_NAME]\tMonitor trace:\t\t $MONITOR_TRACE"
echo -e "[$SCRIPT_NAME]\tStart end trace:\t $START_END_TRACE"
echo -e "[$SCRIPT_NAME]\tIdle time:\t\t $TIME_INACTIVE"
echo -e "[$SCRIPT_NAME]\tMax pipeline processes:\t $MAX_PROC"
echo -e "[$SCRIPT_NAME]\tPipeline threads:\t $THREADS"
#== Begin Run ==#

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
                        test -e $START_END_TRACE && rm $START_END_TRACE # Empty log file
                        test -e $FAILED_LIST && rm $FAILED_LIST # Empty log file
                        test -e $TMP_FILE_PATH && rm $TMP_FILE_PATH # Empty log file
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

# Create folders to copy the results (slow5 files logs)
# test -d $MONITOR_PARENT_DIR/slow5         || mkdir $MONITOR_PARENT_DIR/slow5            || exit 1
# echo "[$SCRIPT_NAME] SLOW5 files will be written to $MONITOR_PARENT_DIR/slow5"
# test -d $MONITOR_PARENT_DIR/slow5_logs    || mkdir $MONITOR_PARENT_DIR/slow5_logs       || exit 1
# echo "[$SCRIPT_NAME] SLOW5 p2s individual logs will be written to $MONITOR_PARENT_DIR/slow5_logs"


# # Start realfreq
# echo "[$SCRIPT_NAME] Starting realfreq" | tee $LOG
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
    touch $DUMP_FILE || { echo "[$SCRIPT_NAME] Could not create dump file $DUMP_FILE. Exiting." | tee $LOG; exit 1; }
    echo "[$SCRIPT_NAME] Non realtime conversion of all files in $MONITOR_PARENT_DIR" | tee -a $LOG
    test -e $TMP_FILE_PATH && rm $TMP_FILE_PATH
    find $MONITOR_PARENT_DIR/ -name "*.pod5" | "$PIPELINE_SCRIPT" -g $GUPPY_BIN -r $REF -i $REFIDX -m $MODEL |& tee -a $LOG | \
    catch_bam | realfreq ${server_port_flag} -d $DUMP_FILE -o $OUTPUT_FILE $REF
    > >(tee -a $TMP_FILE_PATH) 2> >(tee -a $LOG)

else # Else assume realtime analysis is desired

    # Monitor the new file creation in pod5 folder and execute realtime f5-pipeline script
    # Close after timeout met
    if $resuming; then # If resuming option set
        echo "[$SCRIPT_NAME] resuming" | tee -a $LOG
        if [ ! -e $DUMP_FILE ]; then
            echo "[$SCRIPT_NAME] Dump file $DUMP_FILE not found. Exiting." | tee -a $LOG
            exit 1
        fi
        touch $DUMP_FILE || { echo "[$SCRIPT_NAME] Could not create dump file $DUMP_FILE. Exiting." | tee -a $LOG; exit 1; }
        "$SCRIPT_PATH"/monitor/monitor.sh -t $TIME_INACTIVE -f -d ${MONITOR_TEMP} $MONITOR_PARENT_DIR/  |
        "$SCRIPT_PATH"/monitor/ensure.sh -r -d $TMP_FILE_PATH -l ${MONITOR_TRACE}  |
        "$PIPELINE_SCRIPT" -d $TMP_FILE_PATH -l $START_END_TRACE -f $FAILED_LIST -p $MAX_PROC -g $GUPPY_BIN -r $REF -i $REFIDX -m $MODEL |& tee $LOG | \
        catch_bam | realfreq ${server_port_flag} -d $DUMP_FILE -o $OUTPUT_FILE -r -l $TMP_FILE_PATH $REF || { echo "[$SCRIPT_NAME] Realfreq failed processing $BAM_FILE. Try resuming." | tee -a $FAILED_LIST | tee -a $LOG; exit 1; }
        > >(tee -a $TMP_FILE_PATH) 2> >(tee -a $LOG)
    else
        touch $DUMP_FILE || { echo "[$SCRIPT_NAME] Could not create dump file $DUMP_FILE. Exiting." | tee $LOG; exit 1; }
        echo "[$SCRIPT_NAME] running" | tee -a $LOG
        test -e $TMP_FILE_PATH && rm $TMP_FILE_PATH
        "$SCRIPT_PATH"/monitor/monitor.sh -t $TIME_INACTIVE -f -d ${MONITOR_TEMP} $MONITOR_PARENT_DIR/  |
        "$SCRIPT_PATH"/monitor/ensure.sh -d $TMP_FILE_PATH -l ${MONITOR_TRACE}  |
        "$PIPELINE_SCRIPT" -d $TMP_FILE_PATH -l $START_END_TRACE -f $FAILED_LIST -p $MAX_PROC -g $GUPPY_BIN -r $REF -i $REFIDX -m $MODEL |& tee $LOG | \
        catch_bam | realfreq ${server_port_flag} -d $DUMP_FILE -o $OUTPUT_FILE -l $TMP_FILE_PATH $REF
        > >(tee -a $TMP_FILE_PATH) 2> >(tee -a $LOG)
    fi
    if [ ! -e $DUMP_FILE ]; then
        echo "[$SCRIPT_NAME] Dump file $DUMP_FILE not found. Exiting." | tee -a $LOG
        exit 1
    fi
    touch $DUMP_FILE || { echo "[$SCRIPT_NAME] Could not create dump file $DUMP_FILE. Exiting." | tee -a $LOG; exit 1; }
    echo "[$SCRIPT_NAME] No new pod5 files found in last ${TIME_INACTIVE} seconds." | tee -a $LOG
    echo "[$SCRIPT_NAME] converting left overs" | tee -a $LOG
    find $MONITOR_PARENT_DIR/ -name "*.pod5"   |
    "$SCRIPT_PATH"/monitor/ensure.sh -r -d $TMP_FILE_PATH -l ${MONITOR_TRACE}  |
    "$PIPELINE_SCRIPT" -d $TMP_FILE_PATH -l $START_END_TRACE -f $FAILED_LIST -p $MAX_PROC -g $GUPPY_BIN -r $REF -i $REFIDX -m $MODEL |& tee -a $LOG | \
    catch_bam | realfreq ${server_port_flag} -d $DUMP_FILE -o $OUTPUT_FILE -r -l $TMP_FILE_PATH $REF
    > >(tee -a $TMP_FILE_PATH) 2> >(tee -a $LOG)

fi

test -e $FAILED_LIST && echo -e $RED"[$SCRIPT_NAME] $(wc -l $FAILED_LIST) files failed the pipeline. See $FAILED_LIST for the list"$NORMAL | tee -a $LOG
NUMPOD5=$(find $MONITOR_PARENT_DIR/ -name '*.pod5' | wc -l)
NUMBAM=$(find $MONITOR_PARENT_DIR/ -name '*.remora.bam' | wc -l)
if [ ${NUMPOD5} -ne ${NUMBAM} ] ; then
    echo -e $RED"[$SCRIPT_NAME] In $MONITOR_PARENT_DIR, $NUMPOD5 pod5 files, but only $NUMBAM bam files. Check the logs for any failures."$NORMAL | tee -a $LOG
else
    echo "[$SCRIPT_NAME] In $MONITOR_PARENT_DIR, $NUMPOD5 pod5 files, $NUMBAM blow5 files." | tee -a $LOG
fi
# POD5_SIZE=$(find $MONITOR_PARENT_DIR/ -name '*.pod5' -printf "%s\t%p\n" | awk 'BEGIN{sum=0}{sum=sum+$1}END{print sum/(1024*1024*1024)}')
# BAM_SIZE=$(find $MONITOR_PARENT_DIR/ -name '*.remora.bam' -printf "%s\t%p\n" | awk 'BEGIN{sum=0}{sum=sum+$1}END{print sum/(1024*1024*1024)}')
# SAVINGS=$(echo $POD5_SIZE - $BAM_SIZE | bc)
# SAVINGS_PERCENT=$(echo "scale=2; $SAVINGS/$POD5_SIZE*100" | bc)
# echo "POD5 size: $POD5_SIZE GB" | tee -a $LOG
# echo "BAM size: $BAM_SIZE GB" | tee -a $LOG
# echo "Savings: $SAVINGS GB ($SAVINGS_PERCENT%)" | tee -a $LOG

echo "Scanning for errors in log files" | tee -a $LOG
find $MONITOR_PARENT_DIR/ -name '*.log' -exec cat {} \; | grep -i "ERROR" | tee -a $LOG
echo "Scanning for warnings in log files" | tee -a $LOG
find $MONITOR_PARENT_DIR/ -name '*.log' -exec cat {} \; | grep -i "WARNING" | tee -a $LOG
echo "[$SCRIPT_NAME] exiting" | tee -a $LOG
