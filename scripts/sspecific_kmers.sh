#!/bin/bash

THREADS=$(echo "$(nproc) /4" | bc )
THREADS=$([[ ${THREADS} -ge "4" ]] && echo "${THREADS}" || echo "4") #ensure not less than 4
#MEM=$(echo "$(grep MemTotal /proc/meminfo | tr -s ' ' ' ' | cut -f2 -d' ') /4000000" | bc) #set max mem to 1/4 physical

OVERWRITE=false
SUFFIX1="sample-specific"
SUFFIX2="sample-specific"
MIN_FREQ1=2
MIN_FREQ2=2
MIN_FREQ_OUT1=2
MIN_FREQ_OUT2=2
MAX_FREQ_OUT1=9999999
MAX_FREQ_OUT2=9999999
SAMPLE2_FILES=() #when used just as filter, multiple sample2 files can be specified
FILTER_MODE=false

#printf format for messages
FORMAT="%s %-20s %-10s %10s %s\n"
function report {
  printf "${FORMAT}" "`date +"%Y-%m-%d %a %T"`" "[$(basename $0)]" "[${1}]" "${2}" >&2
}

#Exit if any command fails
set -e

##############################
# Parse command line options #
##############################
usage="USAGE: $(basename $0) [-h] [-i sample_1_kmers_db] [-I sample_2_kmers_db]

where:
  -h Show this helpful help text, and exit
  -i <string>  sample 1 input k-mers file name                          [REQUIRED]
  -I <string>  sample 2 input k-mers file name(s)                       [REQUIRED]
  -f <int>     min input frequency threshold for sample 1 k-mers  (default: ${MIN_FREQ1})
  -F <int>     min input frequency* threshold for sample 2 k-mers  (default: ${MIN_FREQ2})
  -q <int>     min output frequency threshold for sample 1 k-mers  (default: ${MIN_FREQ_OUT1})
  -Q <int>     min output frequency threshold for sample 2 k-mers  (default: ${MIN_FREQ_OUT2})
  -x <int>     max output frequency threshold for sample 1 k-mers  (default: ${MAX_FREQ_OUT1})
  -X <int>     max output frequency threshold for sample 2 k-mers  (default: ${MAX_FREQ_OUT2})
  -R           filter sample 1 only (by subracting kmers also present in sample 2 DB(s))
  -s <string>  output name suffix for sample 1 output file (default: ${SUFFIX1})
  -S <string>  output name suffix for sample 2 output file (default: ${SUFFIX2})
  -t <int>     number of threads for kmc_tools (defaults to ${THREADS})
  -O           overwrite existing database


     * - To use different min frequencies (in filterig mode) for individual sample 2 DBs use e.g. -F \"1 10 4\",
       these must be in the same order as the input DBs
  "


while getopts ":hi:I:f:F:q:Q:X:x:Rs:S:t:O" opt; do
#   echo "PARSING -${opt} ${OPTARG}"
  case $opt in
    h) echo "$usage"
       exit;;
    i) SAMPLE1_FILE=${OPTARG};;
    I) SAMPLE2_FILES+=(${OPTARG});;
    f) MIN_FREQ1=${OPTARG};;
    F) MIN_FREQ2=${OPTARG};;
    q) MIN_FREQ_OUT1=${OPTARG};;
    Q) MIN_FREQ_OUT2=${OPTARG};;
    x) MAX_FREQ_OUT1=${OPTARG};;
    X) MAX_FREQ_OUT2=${OPTARG};;
    R) FILTER_MODE=true;;
    s) SUFFIX1=${OPTARG};;
    S) SUFFIX2=${OPTARG};;
    t) THREADS=${OPTARG};;
    O) OVERWRITE=true;;
    ?) printf "Illegal option: '-%s'\n" "${OPTARG}" >&2
       echo "$usage" >&2
       exit 1;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "$usage" >&2
      exit 1;;
  esac
done

if [ -z ${SAMPLE1_FILE} ] || [ -z ${SAMPLE2_FILES} ] ; then
  report "ERROR" "Required option not specified by the user, terminating!!!" >&2
  echo -e "${usage}" >&2
  exit 1
fi

if [ ! -f ${SAMPLE1_FILE} ] && [ ! -f "${SAMPLE1_FILE}.kmc_pre" ]; then
  report "ERROR" "Input DB file ${SAMPLE1_FILE} not found, terminating!" >&2
  exit 1
fi

for f in ${SAMPLE2_FILES[@]}; do
  if [ ! -f ${f} ] && [ ! -f "${f}.kmc_pre" ]; then
    report "ERROR" "Input DB file ${f} not found, terminating!" >&2
    exit 1
  fi
done

if [ ${#SAMPLE2_FILES[@]} -gt 1 ] && [ FILTER_MODE == false   ]; then
  report "ERROR" "Only one file allowed for sample 2 when not in filter mode (label 2 is set as well)" >&2
  exit 1
fi
#######################################################################################
#
#######################################################################################

OUTTMP1=${SAMPLE1_FILE%%.*}
OUTFILE1=${SAMPLE1_FILE%%.*}_${SUFFIX1}.db
BNAME1=${SAMPLE1_FILE##*/}
DIR1=${SAMPLE1_FILE%${BNAME1}}

mkdir -p ${DIR1}logs
PREFIX1=${SAMPLE1_FILE%%.*}
##if filtering only (single output)
if [ -z  ] ; then
  for f in ${SAMPLE2_FILES[@]}; do
    F=$(basename ${f})
    PREFIX2+="_"${F%%.*}
    FFILES+=" "${F%%.*}
  done
else
  F=${SAMPLE2_FILES%%.*}
  PREFIX2="_"${F##*/}
  FFILES+=${F##*/}
fi
#LOGFILE="${DIR1}logs/${PREFIX1##*/}_vs${PREFIX2}.sspecific.log"

#report "INFO" "The logfile is ${LOGFILE}" | tee -a ${LOGFILE}

if [[ -f "${OUTFILE1}.kmc_pre" ]] && [[ "${OVERWRITE}" = false ]]; then
  report "WARNING" "${OUTFILE1} already exists, use the -O flag to overwrite " #| tee -a ${LOGFILE} >&2
  exit 0
fi



##BACKUP PREVIOUS LOGFILE
if [[ -f "${LOGFILE}" ]]; then
  mv ${LOGFILE} ${LOGFILE%.log}."$(date +"%Y-%m-%d-%H:%m")".bak
fi

#echo "All sample2 files: "${SAMPLE2_FILES[@]}

#TODO FFILES NOT REPORTING ALL FILES!!!
report "INFO" "Identifying sample-specific k-mers for: ${SAMPLE1_FILE##*/} vs ${FFILES}" #| tee -a ${LOGFILE}


SAMPLE2_FILE=${SAMPLE2_FILES[0]}
OUTFILE2=${SAMPLE2_FILE%%.*}_${SUFFIX2}.db
if [[ -f "${OUTFILE2}.kmc_pre" ]] && [[ "${OVERWRITE}" = false ]]; then
  report "WARNING" "${OUTFILE2} already exists, use the -O flag to overwrite " #| tee -a ${LOGFILE} >&2
  exit 0
fi

report "INFO" "${SAMPLE1_FILE##*/} MinIn=${MIN_FREQ1}, Min_Out=${MIN_FREQ_OUT1}, MaxOut=${MAX_FREQ_OUT1}" #| tee -a ${LOGFILE}

DB1=${SAMPLE1_FILE%.kmc_???}
DB2=${SAMPLE2_FILE%.kmc_???}


if [ ${FILTER_MODE} == true ]; then
  report "INFO" "${#SAMPLE2_FILES[@]} filter DBs MinIn=${MIN_FREQ2}" #| tee -a ${LOGFILE}

  report "INFO" "Filtering mode... " #| tee -a ${LOGFILE}
  read -r -a MIN_FREQ2_ARR <<< "${MIN_FREQ2}"
  if [[ ${#MIN_FREQ2_ARR[@]} != ${#SAMPLE2_FILES[@]} ]]; then
    report "INFO" "Num filters ${#SAMPLE2_FILES[@]} != ${#MIN_FREQ2_ARR} num freq thresholds..." #| tee -a ${LOGFILE}
    if [[ ${#MIN_FREQ2_ARR[@]} != 1 ]]; then
        report "[ERROR]" "Mismatch between the number of filters ${#SAMPLE2_FILES[@]} and the number of min frequency values ${#MIN_FREQ2_ARR}"
        #| tee -a ${LOGFILE}
        exit 1
    fi
    MIN_FREQ2_ARR=()
    for f in ${SAMPLE2_FILES[@]}; do
      MIN_FREQ2_ARR+=("MIN_FREQ2")
    done
  fi

  COMPLEX=()
  COMPLEX+=("INPUT:")
  COMPLEX+=("\nsetOne=${DB1}" "-ci${MIN_FREQ1}")
  COMMAND="${OUTFILE1} = setOne"
  for i in ${!SAMPLE2_FILES[@]}; do
    f=${SAMPLE2_FILES[${i}]}
    MIN_FREQ2=${MIN_FREQ2_ARR[${i}]}
#    report "INFO" "Still filtering mode... f=${MIN_FREQ2} " | tee -a ${LOGFILE}
    COMPLEX+=("\nset${i}=${f%.kmc_???}" "-ci${MIN_FREQ2}")
    COMMAND+=" - set${i}"
#    report "INFO" "Still filtering mode... " | tee -a ${LOGFILE}
  done

  COMPLEX+=("\nOUTPUT:")
  COMPLEX+=("\n${COMMAND}")
  COMPLEX+=("\nOUTPUT_PARAMS:")
  COMPLEX+=("\n-ci${MIN_FREQ_OUT2}" "-cx${MAX_FREQ_OUT2}")
#  echo -e ${COMPLEX[@]} > complex.tmp

  tput sc #store coursor position, later clear to beginning of line (tput el1) and restore position(tpur rc)
  set -o pipefail && ionice -c 3 kmc_tools -t${THREADS} complex <(echo -e ${COMPLEX[@]} | tee complex) \
  && (tput el1 && tput rc && report "INFO" "Finished filtering k-mers to: ${OUTFILE1} " | tee -a ${LOGFILE} ) \
  || (report "ERROR" "Failed filtering k-mers" | tee -a ${LOGFILE} && exit 1)
else
  report "INFO" "${SAMPLE2_FILE##*/} MinIn=${MIN_FREQ2}, Min_Out=${MIN_FREQ_OUT2}, MaxOut=${MAX_FREQ_OUT2}" | tee -a ${LOGFILE}
  report "INFO" "Running kmc_tools kmers_subtract... "
#tput sc #store coursor position, later clear to beginning of line (tput el1) and restore position(tpur rc)
  set -o pipefail && \
  ionice -c 3 kmc_tools simple \
  ${DB1} -ci${MIN_FREQ1} \
  ${DB2} -ci${MIN_FREQ2} \
  kmers_subtract ${OUTFILE1} -ci${MIN_FREQ_OUT1} -cx${MAX_FREQ_OUT1} \
  reverse_kmers_subtract ${OUTFILE2} -ci${MIN_FREQ_OUT2} -cx${MAX_FREQ_OUT2} \
  && (report "INFO" "Finished identifying sample-specific k-mers, DBs: ${OUTFILE1} ${OUTFILE2}") \
  || (report "ERROR" "Failed identifying sample-specific k-mers" && exit 1)
fi















