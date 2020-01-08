#!/bin/bash

EXTENDEROPTS=()
OVERWRITE=false
SUFFIX="extended"

Xms="1G"
Xmx="150G"

#printf format for messages
FORMAT="%s %-20s %-10s %10s %s\n"
function report {
  printf "${FORMAT}" "`date +"%Y-%m-%d %a %T"`" "[$(basename $0)]" "[${1}]" "${2}" >&2
}

##############################
# Parse command line options #
##############################
usage="USAGE: $(basename $0) [-h] [-i sample_1_kmers_file]

where:
  -h Show this helpful help text, and exit
  -H Show yakat kextender help, and exit
  -i <string>      input k-mers KMC DB file name                                       [REQUIRED]
  -S <string>      suffix of the output file (default: ${SUFFIX})
  -e \"options\"     passess options to yakat kextend
  -s <size>        set initial Java heap size (defaults to ${Xms})
  -m <size>        set maximum Java heap size (defaults to ${Xmx})
  -f <int>         min frequency for a k-mer to be passed to the extender (default: whatever is present in DB)
  -F <int>         max frequency for a k-mer to be passed to the extender (default: whatever is present in DB)
  -O               overwrite existing output file

Example:
$(basename $0) -i sample1_kmers.gz -e \"--threads 12 --output-fasta --name-prefix sample1\"
"

while getopts ":hHi:s:m:e:S:f:F:O" opt; do
  case $opt in
    h) echo -e "$usage"
       exit;;
    H) ${YAKAT} kextend -h
       exit;;
    i) SAMPLE1_FILE=${OPTARG};;
    f) MIN_FREQ="-ci"${OPTARG};;
    F) MAX_FREQ="-cx"${OPTARG};;
    O) OVERWRITE=true;;
    e) EXT_OPT=${OPTARG};;
    S) SUFFIX=${OPTARG};;
    s) Xms=${OPTARG};;
    m) Xmx=${OPTARG};;
    ?) printf "Illegal option: '-%s'\n" "${OPTARG}" >&2
       echo "$usage" >&2
       exit 1;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "$usage" >&2
      exit 1;;
  esac
done

if [ -z ${SAMPLE1_FILE} ] ; then
  report "ERROR" "Required option not specified by the user, terminating!!!" >&2
  echo -e "${usage}" >&2
  exit 1
fi


if [ ! -f ${SAMPLE1_FILE} ]; then
  report "ERROR" "Input file ${SAMPLE1_FILE} not found, terminating!" >&2
  exit 1
fi


BNAME1=${SAMPLE1_FILE##*/}
DIR1=${SAMPLE1_FILE%${BNAME1}}
OUTBASE=${BNAME1%%.*}_${SUFFIX}
OUTFILE=${DIR1}${OUTBASE}.fa

if [[ -s "${OUTFILE}" ]] && [[ "${OVERWRITE}" = false ]] ; then
  report "WARNING" "${OUTFILE} already exists, use the -O flag to overwrite "
else

  report "INFO" "Running k-mer extender for ${SAMPLE1_FILE}"

  INPUTDB=${SAMPLE1_FILE%.kmc_???}
  # KEXTEND="java -XX:+UseGCOverheadLimit -XX:+UseNUMA -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_CPUS_PER_TASK:-16} -Xss150M -Xms${Xms} -Xmx${Xmx} -jar ${YAKAT} kextend"
   

  set -o pipefail && kmc_dump ${MIN_FREQ} ${INPUTDB} /dev/stdout \
  | ${YAKAT} kextend --JVM "-XX:+UseGCOverheadLimit -XX:+UseNUMA -XX:+UseParallelGC -XX:ParallelGCThreads=${SLURM_CPUS_PER_TASK:-16} -Xss150M -Xms${Xms} -Xmx${Xmx}" ${EXT_OPT} 1> ${OUTFILE} \
  && report "INFO" "Finished extending k-mers " \
  || (report "ERROR" "Failed  extending k-mers" && rm  ${OUTFILE} && exit 1)
fi