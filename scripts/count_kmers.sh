#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#COUNT settings
LOWER=2 #4 # Don't output k-mer with count < lower
MAX_RAM=100

THREADS_KMC=10
MAX_COUNT=1000000
OVERWRITE=false
JELLYFISH_NOT_KMC=false
FASTQ_FILES=()

#printf format for messages
FORMAT="%s %-20s %-10s %10s %s\n"
function report {
  printf "${FORMAT}" "`date +"%Y-%m-%d %a %T"`" "[$(basename $0)]" "[${1}]" "${2}" >&2
}

##############################
# Parse command line options #
##############################
usage="\nUSAGE: $(basename $0) [-h] -k <int> -i <fastq.gz> -b <basename> [options]

where:
  -h Show this helpful help text, and exit
  -k <int>     k-mer lenght                                                               [REQUIRED]
  -i <string>  input either FASTQ or FASTA file name(s) (can be specified multiple times) [REQUIRED]
  -b <string>  specify basename for the output files                                      [REQUIRED]
  -d <string>  output directory (defaults to ./k-mers with k set to the user specified value)
  -D <string>  temporary files directory (defaults to output directory)
  -r           RAM-only mode
  -m <int>     maximal value of a counter (defaults to ${MAX_COUNT}, set to 255 to preserve memory)
  -t <int>     number of threads for k-mer counting (defaults to ${THREADS_KMC})
  -S <int>     max memory in GB (defaults to ${MAX_RAM})
  -L <int>     min frequency of a kmer to be reported (defaults to ${LOWER})
  -U <int>     max frequency of a kmer to be reported (not used by default)
  -O           overwrite existing database


May also try escaped wildcards:
  $(basename $0) -i fastq.r\\?.gz -b memorable_string"



while getopts ":hi:L:U:b:m:d:D:k:S:t:Or" opt; do
  case $opt in
    h) echo "$usage"
       exit;;
    i) FASTQ_FILES+=($OPTARG);;
    L) LOWER=${OPTARG};;
    U) UPPER=${OPTARG};;
    b) BNAME=${OPTARG};;
    d) DIR=${OPTARG};;
    D) TMPDIR=${OPTARG};;
    k) k=${OPTARG};;
    m) MAX_COUNT=${OPTARG};;
    S) MAX_RAM=${OPTARG};;
    t) THREADS_KMC=${OPTARG};;
    O) OVERWRITE=true;;
    r) RAM_ONLY=true;;
    ?) printf "Illegal option: '-%s'\n" "${OPTARG}" >&2
       echo "$usage" >&2
       exit 1;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "$usage" >&2
      exit 1;;
  esac
done


if [ ${#FASTQ_FILES[@]} = 0 ]; then
  report "ERROR" "No input FASTQ file(s) provided by the user, terminating!!!" >&2
  echo -e "${usage}" >&2
  exit 1
fi

if [ -z "${BNAME}" ]; then
  report "ERROR" "Please specify the basename for the output files using -b memorable_basename, terminating!!!" >&2
  echo -e "${usage}" >&2
  exit 1
fi


if [ -z "${k}" ]; then
  report "ERROR" "Please specify the value of k using -k someInt, terminating!!!" >&2
  echo -e "${usage}" >&2
  exit 1
fi

for f in "${FASTQ_FILES[@]}"
do
  if [ ! -f ${f} ]; then
    report "ERROR" "Input file ${f} not found, terminating!"
    exit 1
  fi
done


#IF OUT DIR NOT SPECIFIED USE DEFAULT
if [ -z "${DIR}" ]; then
  DIR="${k}-mers"
fi
if [ -z "${TMPDIR}" ]; then
  TMPDIR=${DIR}
fi
mkdir -p ${DIR} ${TMPDIR}

#ls ${TMPDIR} || exit 1

kmers_db="${DIR}/${k}-mers_${BNAME}.db"

#######################################################################################
# Count k-mers
#######################################################################################

mkdir -p ${DIR}/logs
if [[ "${JELLYFISH_NOT_KMC}" = true ]]; then
  LOGFILE="${DIR}/logs/${k}-mers_${BNAME}_${k}-mer_JELLY_count.log"
else
  LOGFILE="${DIR}/logs/${k}-mers_${BNAME}_KMC_count.log"
fi



if [[ -s "${kmers_db}.kmc_pre" ]] && [[ "${OVERWRITE}" = false ]]; then
  report "WARNING" "${kmers_db} already exists, to re-run kmer counting for k=${k}, sample=${BNAME} use -O flag"
  exit 0
fi

##BACKUP PREVIOUS LOGFILE
if [[ -s "${LOGFILE}" ]]; then
  mv ${LOGFILE} ${LOGFILE%.log}."$(date +"%Y-%m-%d-%H:%m")".bak
fi

##RECORDING IN LOG (1) current help (2) user input
echo -e "$usage \nUSER INPUT: "${ALLARGS[@]}"\n" >> ${LOGFILE}


report "INFO" "${#FASTQ_FILES[@]} input file(s) specified" | tee -a ${LOGFILE}


report "INFO" "Counting ${k}-mers in ${BNAME}, db file: ${kmers_db}" | tee -a ${LOGFILE}

  report "INFO" "MIN_FREQ=${LOWER} MAX_FREQ=${UPPER} MAX_COUNT=${MAX_COUNT}" THREADS=${THREADS_KMC} | tee -a ${LOGFILE}

  UPPER_VAR=""
  if [ ! -z "${UPPER}" ]; then
    UPPER_VAR="-cx${UPPER} "
  fi

  if [ "${FASTQ_FILES[0]##*.}" != "gz" ]; then
    FTYPE="-fm"
  fi

  if [[ "${RAM_ONLY}" = true ]]; then
    R="-r"
  else
    R="-sm"
  fi

  set -o pipefail && ionice -c 3 kmc ${FTYPE} -k${k} -m${MAX_RAM} -ci${LOWER} ${UPPER_VAR} \
  -cs${MAX_COUNT} -t${THREADS_KMC} -v ${R} \
  @<(for f in ${FASTQ_FILES[@]}; do echo ${f}; done) \
  ${kmers_db} ${TMPDIR}  2>&1 | tr -d '*' >> ${LOGFILE} \
  && report "INFO" "Finished counting ${k}-mers for sample ${BNAME}, see ${LOGFILE} for details" | tee -a ${LOGFILE} \
  || (report "ERROR" "Failed counting ${k}-mers for sample ${BNAME}, see ${LOGFILE} for details" | tee -a ${LOGFILE} && exit 1)
