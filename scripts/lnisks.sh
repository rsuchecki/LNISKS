#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


#UPDATE AS NECESSARY
export YAKAT="${DIR}/yakat.jar"
export PATH="${PATH}:"$(pwd)"/bin"

#RESOURCES
THREADS=$(echo "$(nproc) /4" | bc )
THREADS=$([[ ${THREADS} -ge "4" ]] && echo "${THREADS}" || echo "4") #ensure not less than 4
MEM=$(echo "$(grep MemTotal /proc/meminfo | tr -s ' ' ' ' | cut -f2 -d' ') /4000000" | bc) #set max mem to 1/4 physical
#MEMSTART=$(echo "${MEM} /10" | bc ) #initial mem for JVMm, aiming for 1/10 of max set
#MEMSTART=$([[ ${MEMSTART} -ge "1" ]] && echo "${MEMSTART}" || echo "1") #ensure not less than one GB


OUTDIR=output
OVERWRITE_FROM_STEP=99 #not that many steps so not applied
STOP_AFTER_STEP=7 ##further extension of paired seeds disabled by default
SKIP_FIRST_STEPS=0

###
MT_FILES=()
MT_NAME="mutant"
MT_HET=false
MT_MIN_FREQ_IN=2
MT_MIN_FREQ_OUT=5
MT_MAX_FREQ_OUT=9999999
MT_INFER_MIN_FREQ=false
MT_FILTER_MIN_FREQ=2
FILTER_MT_FILES=()

WT_FILES=()
WT_NAME="wildtype"
WT_HET=false
WT_MIN_FREQ_IN=2
WT_MIN_FREQ_OUT=5
WT_MAX_FREQ_OUT=9999999
WT_INFER_MIN_FREQ=false
WT_FILTER_MIN_FREQ=2
FILTER_WT_FILES=() ##Not needed in typical analysis


HT_FILES=()
HT_NAME="hets"
HT_MIN_FREQ_IN=2
HT_MIN_FREQ_OUT=5
HT_MAX_FREQ_OUT=9999999
HT_INFER_MIN_FREQ=false
HT_FILTER_MIN_FREQ=2
FILTER_HT_FILES=()

#Further extensions of paired seeds

KEXTEND_K_STEP=10
INFER_MIN_FREQ_EXTEND=false

#Long, unpaired
MIN_LONG=250


if [ -z ${COLUMNS} ]; then
  COLUMNS=160
fi


#Exit if any command fails
set -e

#printf format for messages
FORMAT="%s %-20s %-10s %10s %s\n"
function report {
  printf "${FORMAT}" "`date +"%Y-%m-%d %a %T"`" "[$(basename $0)]" "[${1}]" "${2}" >&2
}

##############################
# Parse command line options #
##############################


usage="USAGE: $(basename $0) [-h] -k <int> -j <int> -M <fastq.gz> -W <fastq.gz> [options]
[Input/Output settigns]
  -k <int>        k-mer lenght                                             [REQUIRED]
  -j <int>        j-mer length, j<k or even j<<k ^^                        [Step skipped if not given]
  -o <string>     output directory (default: ${OUTDIR})
  -M <string>     sample 1 input fastq.gz file* [REQUIRED if not skipping step 1 (-s n, n>0)]
  -W <string>     sample 2 input fastq.gz file* [REQUIRED if not skipping step 1 (-s n, n>0)]
  -m <string>     sample 1 name/label (default: ${MT_NAME})
  -w <string>     sample 2 name/label (default: ${WT_NAME})
  -Q <int>        sample 1 min input k-mer frequency (default: ${MT_MIN_FREQ_IN})
  -q <int>        sample 2 min input k-mer frequency (default: ${WT_MIN_FREQ_IN})
  -P <int>        sample 1 min output k-mer frequency (default: ${MT_MIN_FREQ_OUT})
  -p <int>        sample 2 min output k-mer frequency (default: ${WT_MIN_FREQ_OUT})
  -X <int>        sample 1 max output k-mer frequency (default: ${MT_MAX_FREQ_OUT})
  -x <int>        sample 2 max output k-mer frequency (default: ${WT_MAX_FREQ_OUT})
  -L <int>        Output long, unpaired seeds of no less than <int> bp (default: ${MIN_LONG})^^
[Simple procedure for inferring k-mer frequency - nonsense if input is RNA-Seq data]
  -I              sample 1 infer min out k-mer frequency from distribution
  -i              sample 2 infer min out k-mer frequency from distribution
[Further filtering options for k-mer DBs]
  -F <string>     filter DB(s)*^ for sample 1
  -f <string>     filter DB(s)*^ for sample 2
  -Y <int>        min frequency** for sample 1 filter DB(s) (default: ${MT_FILTER_MIN_FREQ})
  -y <int>        min frequency** for sample 2 filter DB(s) (default: ${WT_FILTER_MIN_FREQ})
[Pairing seeds and calling SNPs]
  -d <float>      min identity (0.0<=id<=1.0) required when clustering/paring seeds (default: 1-4/(2k-1))
  -D <float>      min identity (0.0<=id<=1.0) required when calling snps (default: 1-2/(2k-1))
[Further extensions of paired seeds - disabled by default, see -S]
  -B              k-merize all input reads (no baiting), suitable for small input datasets and/or big memory machines
  -J              infer min k-mer frequency for further extensions from distribution
  -b <int>        baiting k-mer length (default equals j-mer length)
  -n <int>        min extension k-mer size (default equals k)
  -N <int>        max extension k-mer size (default equals 2k)
  -e <int>        extension k-mer size step (default: ${KEXTEND_K_STEP})
[Runtime, skipping, stopping]
  -C <int>        Print width for some reports, recommended use: -C \$COLUMNS (defaults to ${COLUMNS})
  -t <int>        number of threads for parallelized tasks (defaults to max(4,nproc/4: ${THREADS})
  -T <string>     temporary files directory for KMC
  -E <int>        max physical memory in GB to be used (defaults to 1/4 of physical mem: ${MEM})
  -O <int>        overwrite existing output files starting from step:
                   1 k-mer counting
                   2 k-mer min frequency inferring [skipped by default, use {-I, -i} to run ]
                   3 sample-specific k-mer identification
                   4 custom filtering of sample-specific k-mers
                   5 k-mer extension to seeds (unitigs)
                   6 restriction of extended seeds to those sharing j-mers
                   7 matching extended seeds and SNP calling
                   8 baiting reads which match the paired sequences (step ignored if -B is used)
                   9 kmerizing baited (or all if -B was used) input reads
                  10 extending matched seeds (ideally beyond the sample-specific 2k-1 bp)
  -s <int>        Skip first <int> steps of the pipeline (steps listed above, default: ${SKIP_FIRST_STEPS})
  -S <int>        Stop pipeline after step <int> (steps listed above, default: ${STOP_AFTER_STEP})

   * - Can be specified multiple times, and/or use escaped wildcards e.g. -M filename_R\?.fastq.gz
  ^^ - The set of unpaired sequences will be incomplete if using j-mers for speeding up pairing
  *^ - Filters are KMC DBs which can be generated using count_kmers.sh wrapper around KMC
  ** - To use different min frequencies for individual filters use e.g. -Y \"1 10 4\",
       these must be in the same order as the input filters
"
#  -g <string>     additional heterozygous sample name/label (default: ${HT_NAME}) [DISABLED]
#  -G <string>     additional heterozygous sample input fastq.gz file*             [DISABLED]
#  -Z              add pre-filtered hets to sample 1                               [DISABLED]
#  -z              add pre-filtered hets to sample 2                               [DISABLED]



ALLARGS=()
while getopts ":hiIzk:j:o:M:m:W:w:Q:q:P:p:X:x:F:f:Y:y:d:D:BJb:n:N:e:t:T:C:E:O:L:s:S:" opt; do
  ALLARGS+=("-"${opt} ${OPTARG})
#  echo "PARSING -${opt} ${OPTARG}"
  case $opt in
    h) echo "$usage"
       exit;;
    k) k=${OPTARG};;
    j) j=${OPTARG};;
    o) OUTDIR=${OPTARG};;
    M) MT_FILES+=(${OPTARG});;
    m) MT_NAME=${OPTARG};;
    W) WT_FILES+=(${OPTARG});;
    w) WT_NAME=${OPTARG};;
    Z) MT_HET=true;;
    Q) MT_MIN_FREQ_IN=${OPTARG};;
    P) MT_MIN_FREQ_OUT=${OPTARG};;
    X) MT_MAX_FREQ_OUT=${OPTARG};;
    I) MT_INFER_MIN_FREQ=true;;
    q) WT_MIN_FREQ_IN=${OPTARG};;
    p) WT_MIN_FREQ_OUT=${OPTARG};;
    x) WT_MAX_FREQ_OUT=${OPTARG};;
    i) WT_INFER_MIN_FREQ=true;;
    F) FILTER_FILES_MT+=(${OPTARG});;
    f) FILTER_FILES_WT+=(${OPTARG});;
    Y) MT_FILTER_MIN_FREQ=${OPTARG};;
    y) WT_FILTER_MIN_FREQ=${OPTARG};;
    d) MIN_ID1=${OPTARG};;
    D) MIN_ID2=${OPTARG};;
    B) KEXTEND_NO_BAIT=true;;
    J) INFER_MIN_FREQ_EXTEND=true;;
    b) KEXTEND_K_BAIT=${OPTARG};;
    n) KEXTEND_K_MIN=${OPTARG};;
    N) KEXTEND_K_MAX=${OPTARG};;
    e) KEXTEND_K_STEP=${OPTARG};;
    C) COLUMNS=${OPTARG};;
    t) THREADS=${OPTARG};;
    T) TMPDIR=${OPTARG};; ##TODO USE mkdtemp
    E) MEM=${OPTARG};;
    L) MIN_LONG=${OPTARG};;
    O) OVERWRITE_FROM_STEP=${OPTARG};;
    s) SKIP_FIRST_STEPS=${OPTARG};;
    S) STOP_AFTER_STEP=${OPTARG};;
    ?) printf "Illegal option: '-%s'\n" "${OPTARG}" >&2
       echo "$usage" >&2
       exit 1;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "$usage" >&2
      exit 1;;
  esac
done

if [ -z ${k} ] ; then
  report "ERROR" "Required option (-k) not specified by the user, terminating!!!" >&2
  exit 1
fi

if [ ! -z ${j} ]  && [ ${j} -ge ${k} ]; then
  report "ERROR" "Consider reducing the value of j. When j >= k none of the extended seeds would be considered for matching, so no variants would be called" >&2
  exit 1
fi

if [ -z ${KEXTEND_K_BAIT} ]; then
  if [ ! -z ${j} ]; then
    KEXTEND_K_BAIT=${j}
  else
    KEXTEND_K_BAIT=${k}
  fi
fi
if [ -z ${KEXTEND_K_MIN} ]; then
  KEXTEND_K_MIN=${k}
fi
if [ -z ${KEXTEND_K_MAX} ]; then
  KEXTEND_K_MAX=$((${k}+${k}))
fi

for arg in ${k} ${OVERWRITE_FROM_STEP} ${SKIP_FIRST_STEPS} ${STOP_AFTER_STEP} ${KEXTEND_K_BAIT} ${KEXTEND_K_MIN} ${KEXTEND_K_MAX} ${KEXTEND_K_STEP} ${COLUMNS}; do
  numeric='^[0-9]+$'
  if ! [[ "${arg}" =~ $re ]] ; then
    report "ERROR" "Numeric value expected, offending argument: ${arg}" >&2
    exit 1
  fi
done



if [ ${SKIP_FIRST_STEPS} -le 0 ]; then
  if [ -z ${WT_FILES} ] || [ -z ${MT_FILES} ] ; then
    report "ERROR" "Required option not specified by the user, terminating!!!" >&2
#    echo -e "$usage \nUSER INPUT: "${ALLARGS[@]}"\n"
#    echo -e "${usage}" >&2
    exit 1
  fi
fi


#Deal with dir input
#TMP_FILE_LIST=()
#for f in ${MT_FILES[@]}; do
#  if [ -d ${f} ]; then
#    TMP_FILE_LIST+=$(readlink -f ${f}/* | fgrep -i fastq.gz | tr '\n' ' ')
#  else
#    TMP_FILE_LIST+="${f} "
#  fi
#done
#MT_FILES=${TMP_FILE_LIST[@]}

#TMP_FILE_LIST=()
#for f in ${WT_FILES[@]}; do
#  if [ -d ${f} ]; then
#    TMP_FILE_LIST+=$(readlink -f ${f}/* | fgrep -i fastq.gz | tr '\n' ' ')
#  else
#    TMP_FILE_LIST+=${f}
#  fi
#done
#WT_FILES=${TMP_FILE_LIST[@]}

for f in ${WT_FILES[@]} ${MT_FILES[@]} ${FILTER_MT_FILES[@]} ${FILTER_WT_FILES[@]} ${HT_FILES[@]}; do
  if [ ! -f ${f} ]; then
    report "ERROR" "Input file ${f} not found, terminating!" >&2
    exit 1
  fi
done

if [ ${MT_HET} == true ] || [ ${WT_HET} == true ]; then
  if [ ${#HT_FILES[@]} -eq 0 ]; then
    report "ERROR" "No het files specified  (offending arguments: -z -Z), terminating!" >&2
    exit 1
  fi
fi

if [ ${#HT_FILES[@]} -gt 0 ]; then
  if [ ${MT_HET} == false ] && [ ${WT_HET} == false ]; then
    report "ERROR" "Het input given without specifying -z or -Z, terminating!" >&2
    exit 1
  fi
fi

if [[ ${MT_NAME} == *['!'@#\$%^\&*()+/.]* ]]; then
  report "ERROR" "Special characters not allowed in sample name (offending argument: -m ${MT_NAME}), terminating!" >&2
  exit 1
fi

if [[ ${WT_NAME} == *['!'@#\$%^\&*()+/.]* ]]; then
  report "ERROR" "Special characters not allowed in sample name (offending argument: -m ${WT_NAME}), terminating!" >&2
  exit 1
fi

WD="${OUTDIR}/${k}-mers"
mkdir -p ${WD}/logs

LOGFILE="${WD}/logs/${k}-mers_lnisks_"$(date +"%Y-%m-%d-%a-%Hh%Mm%Ss").log
exec > >(tee -a ${LOGFILE})
exec 2> >(tee -a ${LOGFILE} >&2)

##RECORDING IN LOG (1) current help (2) user input
echo -e "${0} ${@}\n" >> ${LOGFILE}
echo -e "$usage \nUSER INPUT: "${ALLARGS[@]}"\n" >> ${LOGFILE}


if [ ! -z ${TMPDIR} ]; then
  SETTMP="-D ${TMPDIR}"
fi


report "INFO" "Pipeline initiated"
STEP=1


if [ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]; then
  OV="-O"
fi

#######################################################################################
#XXX COUNT KMERS
#######################################################################################
if [ "${SKIP_FIRST_STEPS}" -ge "${STEP}" ]; then
  report "INFO" "Skipping step ${STEP} at users request"
else
  report "INFO" "Running step ${STEP} (k-mer counting)"
  # WT
  FILES=()
  for F in ${WT_FILES[@]}; do
    FILES+=("-i")
    FILES+=(${F})
  done
  if [ "${F##*.}" == "kmc_pre" ] || [ "${F##*.}" == "kmc_suf" ]; then
    report "INFO" "Not counting k-mers in  ${F} which is a KMC k-mer DB "
    STOP_AFTER_STEP=7
    if [ "${#WT_FILES[@]}" -gt "1" ]; then
      report "ERROR" "No more than one KMC DB allowed as sample input" >&2
      exit 1
    fi
    ln -sr ${F%.*}.kmc_pre ${WD}/${k}-mers_${WT_NAME}.db.kmc_pre
    ln -sr ${F%.*}.kmc_suf ${WD}/${k}-mers_${WT_NAME}.db.kmc_suf
  else
    set -o pipefail && ${DIR}/count_kmers.sh -k ${k} -L ${WT_MIN_FREQ_IN} -S ${MEM} -m 255 \
    ${FILES[@]} -b ${WT_NAME} -t ${THREADS} -d ${WD} ${SETTMP} ${OV} || exit 1
  fi

  # MT
  FILES=()
  for F in ${MT_FILES[@]}; do
    FILES+=("-i")
    FILES+=(${F})
  done
  if [ "${F##*.}" == "kmc_pre" ] || [ "${F##*.}" == "kmc_suf" ]; then
    report "INFO" "Not counting k-mers in  ${F} which is a KMC k-mer DB "
    STOP_AFTER_STEP=7
    if [ "${#MT_FILES[@]}" -gt "1" ]; then
      report "ERROR" "No more than one KMC DB allowed as sample input" >&2
      exit 1
    fi
    ln -sr ${F%.*}.kmc_pre ${WD}/${k}-mers_${MT_NAME}.db.kmc_pre
    ln -sr ${F%.*}.kmc_suf ${WD}/${k}-mers_${MT_NAME}.db.kmc_suf
  else
    set -o pipefail && ${DIR}/count_kmers.sh -k ${k} -L ${MT_MIN_FREQ_IN} -S ${MEM} -m 255 \
    ${FILES[@]} -b ${MT_NAME} -t ${THREADS} -d ${WD} ${SETTMP} ${OV} || exit 1
  fi

#  # HT
#  if [ ${#HT_FILES[@]} -gt 0 ]; then
#    FILES=()
#    for F in ${HT_FILES[@]}; do
#      FILES+=("-i")
#      FILES+=(${F})
#    done
#    set -o pipefail && ${DIR}/count_kmers.sh -k ${k} -L ${HT_MIN_FREQ_IN} -S ${MEM} -m 255 \
#    ${FILES[@]} -b ${HT_NAME} -t ${THREADS} -d ${WD} ${OV} || exit 1
#  fi
fi


if [ "${STOP_AFTER_STEP}" -eq "${STEP}" ]; then
  report "INFO" "Stopped after step ${STEP}"
  exit 0
fi
((STEP++))
if [ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]; then
  OV=true
  #"-O"
fi


#DBs generated
WT_TO_SS=${WD}/${k}-mers_${WT_NAME}.db.kmc_pre
MT_TO_SS=${WD}/${k}-mers_${MT_NAME}.db.kmc_pre
#HT_TO_SS=${WD}/${k}-mers_${HT_NAME}.db.kmc_pre

#######################################################################################
#XXX ESTABLISH k-mer FREQUENCY CUTOFFS
#######################################################################################
## iterate through histogram, storing frequency f, if f(line) > f(line-1) use f(line) as the minimum f allowed
mkdir -p tmp

if [ "${SKIP_FIRST_STEPS}" -ge "${STEP}" ]; then
  report "INFO" "Skipping step ${STEP} at users request"
else
  if [ "${MT_INFER_MIN_FREQ}" == true ]; then
    MT_HISTO=${WD}/${k}-mers_${MT_NAME}.histogram
    report "INFO" "Attempting to estimate the minimum frequency for ${MT_NAME}..."

    if [[ -s "${MT_HISTO}" ]] && [[ "${OV}" != true ]]; then
      report "WARNING" "${MT_HISTO} already exists, use -O ${STEP} to overwrite"
    else
      kmc_tools -hp histogram ${MT_TO_SS%.kmc_pre} ${MT_HISTO} || \
      (report "ERROR" "Failed to estimate the minimum frequency for ${MT_NAME} k-mers to be included" && exit 1)
    fi
    MT_MIN_FREQ_OUT=$(awk -v prev=99999999999 '{if ($2>prev){print $1-1; exit}; prev=$2}' ${MT_HISTO})
    MAX_PRINT1=$(awk -vMIN=${MT_MIN_FREQ_OUT} -vPREV=0 '{if ($1>MIN && $2<prev){print $1-1; exit}; prev=$2}' ${MT_HISTO})
    #MAX_PRINT1=$(sort -k2,2nr ${MT_HISTO} | head -1 | cut -f1)
    ./scripts/plot_histogram.sh ${MT_HISTO} ${COLUMNS} | tee ${MT_HISTO}.plot | awk -vM="${MAX_PRINT1}" 'NR<=(M*3)'
    report "INFO" "Estimated minimum frequency for ${MT_NAME} k-mers to be included: ${MT_MIN_FREQ_OUT}, see ${MT_HISTO}.plot"
  fi

  if [ "${WT_INFER_MIN_FREQ}" == true ]; then
    WT_HISTO=${WD}/${k}-mers_${WT_NAME}.histogram
    report "INFO" "Attempting to estimate the minimum frequency for ${WT_NAME}..."
    if [[ -s "${WT_HISTO}" ]] && [[ "${OV}" != true ]]; then
      report "WARNING" "${WT_HISTO} already exists, use -O ${STEP} to overwrite"
    else
      kmc_tools -hp histogram ${WT_TO_SS%.kmc_pre} ${WT_HISTO} || \
      (report "ERROR" "Failed to estimate the minimum frequency for ${WT_NAME} k-mers to be included" && exit 1)
    fi
    WT_MIN_FREQ_OUT=$(awk -v prev=99999999999 '{if ($2>prev){print $1-1; exit}; prev=$2}' ${WT_HISTO})
    #MAX_PRINT2=$(sort -k2,2nr ${WT_HISTO} | head -1 | cut -f1)
    MAX_PRINT2=$(awk -vMIN=${WT_MIN_FREQ_OUT} -vPREV=0 '{if ($1>MIN && $2<prev){print $1-1; exit}; prev=$2}' ${WT_HISTO})
    ./scripts/plot_histogram.sh ${WT_HISTO} ${COLUMNS} | tee ${WT_HISTO}.plot | awk -vM="${MAX_PRINT2}" 'NR<=(M*2)'
    report "INFO" "Estimated minimum frequency for ${WT_NAME} k-mers to be included: ${WT_MIN_FREQ_OUT}, see ${WT_HISTO}.plot"
  fi
fi

if [ "${STOP_AFTER_STEP}" -eq "${STEP}" ]; then
  report "INFO" "Stopped after step ${STEP}"
  exit 0
fi
((STEP++))
if [ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]; then
  OV="-O"
fi

########################################################################################
## XXX IDENTIFY SAMPLE SPECIFIC k-mers
########################################################################################
if [ "${SKIP_FIRST_STEPS}" -ge "${STEP}" ]; then
  report "INFO" "Skipping step ${STEP} at users request"
else
  report "INFO" "Running step ${STEP} (sample-specific) "

## TODO consider doing sspecific and custom filtering in one 'kmc complex' operation

  set -o pipefail &&  ${DIR}/sspecific_kmers.sh \
  -i ${WT_TO_SS} -f ${WT_MIN_FREQ_IN} -q ${WT_MIN_FREQ_OUT} -x ${WT_MAX_FREQ_OUT} \
  -I ${MT_TO_SS} -F ${MT_MIN_FREQ_IN} -Q ${MT_MIN_FREQ_OUT} -X ${MT_MAX_FREQ_OUT} \
  -s "MINUS_${MT_NAME}" -S "MINUS_${WT_NAME}" \
  -t ${THREADS} ${OV} || exit 1

  sed -i 's/[[:cntrl:]].*//' ${LOGFILE} #remove kmc progress lines, which makes it necessary to re-set exec
  exec > >(tee -a ${LOGFILE})
  exec 2> >(tee -a ${LOGFILE} >&2)


fi

if [ "${STOP_AFTER_STEP}" -eq "${STEP}" ]; then
  report "INFO" "Stopped after step ${STEP}"
  exit 0
fi
((STEP++))
if [ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]; then
  OV="-O"
fi
########################################################################################
## XXX FILTER SAMPLE SPECIFIC k-mers USING ONE OR MORE CUSTOM FILTER(S)
## TODO consider doing sspecific and custom filtering in one 'kmc complex' operation
########################################################################################
WT_TOEXTEND=${WT_TO_SS%%.*}_MINUS_${MT_NAME}.db.kmc_pre
MT_TOEXTEND=${MT_TO_SS%%.*}_MINUS_${WT_NAME}.db.kmc_pre


if [ "${SKIP_FIRST_STEPS}" -ge "${STEP}" ]; then
  report "INFO" "Skipping step ${STEP} at users request"
else
  report "INFO" "Running step ${STEP} (custom filters)"
    #WT
  if [ ${#FILTER_FILES_WT[@]} -gt 0 ]; then
    FILES=()
    for F in ${FILTER_FILES_WT[@]}; do
      FILES+=("-I")
      FILES+=(${F})
    done
    set -o pipefail && ${DIR}/sspecific_kmers.sh -R -i ${WT_TOEXTEND} \
    -f ${WT_MIN_FREQ_IN} -F \"${WT_FILTER_MIN_FREQ}\" -q ${WT_MIN_FREQ_OUT} \
    ${FILES[@]} -s filtered -t ${THREADS} ${OV} || exit 1
  fi
  #MT
  if [ ${#FILTER_FILES_MT[@]} -gt 0 ]; then
    FILES=()
    for F in ${FILTER_FILES_MT[@]}; do
      FILES+=("-I")
      FILES+=(${F})
    done

    set -o pipefail && ${DIR}/sspecific_kmers.sh -R -i ${MT_TOEXTEND} \
    -f ${MT_MIN_FREQ_IN} -F "${MT_FILTER_MIN_FREQ}" -q ${MT_MIN_FREQ_OUT} \
    ${FILES[@]} -s filtered -t ${THREADS} ${OV} || exit 1
  fi
  sed -i 's/[[:cntrl:]].*//' ${LOGFILE} #remove kmc progress lines, which makes it necessary to re-set exec
  exec > >(tee -a ${LOGFILE})
  exec 2> >(tee -a ${LOGFILE} >&2)
fi

#if filtering ON, use filtered dataset for extending
if [ ${#FILTER_FILES_WT[@]} -gt 0 ]; then
  WT_TOEXTEND=${WT_TOEXTEND%%.*}_filtered.db.kmc_pre
fi
if [ ${#FILTER_FILES_MT[@]} -gt 0 ]; then
  MT_TOEXTEND=${MT_TOEXTEND%%.*}_filtered.db.kmc_pre
fi

if [ "${STOP_AFTER_STEP}" -eq "${STEP}" ]; then
  report "INFO" "Stopped after step ${STEP}"
  exit 0
fi
((STEP++))
if [ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]; then
  OV="-O"
fi
########################################################################################
## XXX Extend k-mers
########################################################################################

if [ "${SKIP_FIRST_STEPS}" -ge "${STEP}" ]; then
  report "INFO" "Skipping step ${STEP} at users request"
else
  report "INFO" "Running step ${STEP} (kmer extender)"
  #WT
  set -o pipefail && ${DIR}/extend_kmers.sh -f ${WT_MIN_FREQ_OUT} -i ${WT_TOEXTEND} -m ${MEM}G \
  -e "--threads ${THREADS} --fasta-out --fasta-id-prefix ${WT_NAME}" ${OV} || exit 1
  #MT
  set -o pipefail && ${DIR}/extend_kmers.sh -f ${MT_MIN_FREQ_OUT} -i ${MT_TOEXTEND} -m ${MEM}G \
  -e "--threads ${THREADS} --fasta-out --fasta-id-prefix ${MT_NAME}" ${OV} || exit 1
#  -s ${MEMSTART}G
fi

if [ "${STOP_AFTER_STEP}" -eq "${STEP}" ]; then
  report "INFO" "Stopped after step ${STEP}"
  exit 0
fi
((STEP++))
if [ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]; then
  OV="-O"
fi

########################################################################################
## XXX Reduce extended seeds sets to those that share j-mers - potentailly alignable
########################################################################################
if [ "${SKIP_FIRST_STEPS}" -ge "${STEP}" ]; then
  report "INFO" "Skipping step ${STEP} at users request"
elif [ -z "${j}" ]; then
  report "INFO" "Skipping step ${STEP} as -j not specified by the user"
  MATCHINGSUFFIX=extended.fa
  WT_MATCHING=${WT_TOEXTEND%%.*}_${MATCHINGSUFFIX}
  MT_MATCHING=${MT_TOEXTEND%%.*}_${MATCHINGSUFFIX}
else
  report "INFO" "Running step ${STEP} "
  MATCHINGSUFFIX=j-filtered.fa
  WT_MATCHING=${WT_TOEXTEND%%.*}_${MATCHINGSUFFIX}
  MT_MATCHING=${MT_TOEXTEND%%.*}_${MATCHINGSUFFIX}
  if [[ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]]  || [[ ! -s ${WT_MATCHING} ]]; then
    #WT
    report "INFO" "Attempting to j-filter ${WT_NAME}, j=${j}"

    set -o pipefail && paste - - < ${WT_TOEXTEND%%.*}_extended.fa \
    | java -Xms${MEM}G -Xmx${MEM}G -jar ${YAKAT} kmatch \
    --k-mer-length ${j} --k-mers ${MT_TOEXTEND%%.*}_extended.fa  \
    --threads ${THREADS} --print-user-settings \
    | tr '\t' '\n' > ${WT_MATCHING} \
    && report "INFO" "Finished j-filtering ${WT_NAME}" \
    || (report "ERROR" "Failed j-filtering ${WT_NAME}" && exit 1)
  else
    report "WARNING" "File ${WT_MATCHING} already exists, use -O ${STEP} to overwrite"
  fi
  #MT
  if [[ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]]  || [[ ! -s ${MT_MATCHING} ]]; then
    report "INFO" "Attempting to j-filter ${MT_NAME}, j=${j}"
    set -o pipefail && paste - - < ${MT_TOEXTEND%%.*}_extended.fa \
    | java -Xms${MEM}G -Xmx${MEM}G -jar ${YAKAT} kmatch \
    --k-mer-length ${j} --k-mers ${WT_TOEXTEND%%.*}_extended.fa  \
    --threads ${THREADS} --print-user-settings \
    | tr '\t' '\n' > ${MT_MATCHING} \
    && report "INFO" "Finished j-filtering ${MT_NAME}"  \
    || (report "ERROR" "Failed j-filtering ${MT_NAME}" && exit 1)
  else
    report "WARNING" "File ${MT_MATCHING} already exists, use -O ${STEP} to overwrite"
  fi
fi

if [ "${STOP_AFTER_STEP}" -eq "${STEP}" ]; then
  report "INFO" "Stopped after step ${STEP}"
  exit 0
fi
((STEP++))
if [ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]; then
  OV="-O -o"
fi


########################################################################################
## XXX match extended
########################################################################################
if [ "${SKIP_FIRST_STEPS}" -ge "${STEP}" ]; then
  report "INFO" "Skipping step ${STEP} at users request"
else
  report "INFO" "Running step ${STEP} "

  if [ -z ${MIN_ID1} ]; then
     MIN_ID1=$(awk -v k=$k "BEGIN{print 1-4/(2*k-1)}");
  fi
  report "INFO" "Minimum identity for clustering/pairing set at ${MIN_ID1:0:4} " #XXX to FLOOR: ${MIN_ID:0:4}
  if [ -z ${MIN_ID2} ]; then
     MIN_ID2=$(awk -v k=$k "BEGIN{print 1-2/(2*k-1)}");
  fi
  report "INFO" "Minimum identity for snp calling set at ${MIN_ID2:0:4} " #XXX to FLOOR: ${MIN_ID:0:4}

  if [[ "${WT_NAME}" > "${MT_NAME}" ]]; then
    REV="-R " #ensure SNPs are reported wt->mut not the other way around
  fi


  set -o pipefail && ${DIR}/match_and_call.sh ${REV} -d ${MIN_ID1:0:4} -D ${MIN_ID2:0:4} \
  -i ${WT_MATCHING} \
  -I ${MT_MATCHING} \
  -n ${WT_NAME} \
  -N ${MT_NAME} \
  -C ${COLUMNS} \
  -k ${k} -L ${MIN_LONG} -t ${THREADS} -E ${MEM} ${OV} || exit 1
fi

if [ "${STOP_AFTER_STEP}" -eq "${STEP}" ]; then
  report "INFO" "Stopped after step ${STEP}."
  report "INFO" "If generated contigs (seeds) are not long enough for relaible annotation, re-run with -s 7 -S 10 to extend further."
  exit 0
fi
((STEP++))




if [ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]; then
  OV="-o" #overwrite baiting and the remaining steps
elif [ "${OVERWRITE_FROM_STEP}" -le "$((${STEP}+1))" ]; then
  OV="-O" #overwrite k-mer counting and the remaining steps
elif [ "${OVERWRITE_FROM_STEP}" -le "$((${STEP}+2))" ]; then
  OV="-R" #overwrite extending and the remaining steps
fi

#########################################################################################
### XXX further extend matched k-mers to allow higher specificity for annotation
#########################################################################################
if [ "${SKIP_FIRST_STEPS}" -ge "${STEP}" ]; then
  report "WARNING" "Skipping step ${STEP} implies skipping the next two steps as well"
elif [ "${SKIP_FIRST_STEPS}" -ge "$((${STEP}+1))" ]; then
  report "INFO" "Skipping step "$((${STEP}+1))" implies skipping the ext step as well"
elif [ "${SKIP_FIRST_STEPS}" -ge "$((${STEP}+2))" ]; then
  report "INFO" "Skipping step "$((${STEP}+2))" at users request"
else

  report "INFO" "Running steps ${STEP},"$((${STEP}+1))","$((${STEP}+2))" (use paired seeds to identify matching reads, kmerize, extend) "


  WT_PAIRED=${WT_MATCHING%.fa}.matched.fa  #${WD}/${k}-mers_${WT_NAME}.matched.fa #${WD}/${k}-mers_${WT_NAME}.paired.fa
  MT_PAIRED=${MT_MATCHING%.fa}.matched.fa #${WD}/${k}-mers_${MT_NAME}.matched.fa

  MTTEMP=${MT_MATCHING%.fa}
  BOTH_PAIRED=${WT_MATCHING%.fa}_${MTTEMP##*/}.matched.fa


  fgrep -A1 --no-group-separator ${MT_NAME} ${BOTH_PAIRED} > ${MT_PAIRED}
  fgrep -A1 --no-group-separator ${WT_NAME} ${BOTH_PAIRED} > ${WT_PAIRED}

 ###XXX TODO unpaired (if needed)

#  MT_UNPAIRED=${WD}/${k}-mers_${MT_NAME}.unmatched.fa
#  WT_UNPAIRED=${WD}/${k}-mers_${WT_NAME}.unmatched.fa

#  UNPAIRED=${WD}/${k}-mers_${WT_NAME}_${k}-mers_${MT_NAME}.unmatched.fa

#  fgrep -A1 --no-group-separator ${MT_NAME} ${UNPAIRED} > ${MT_UNPAIRED}
#  fgrep -A1 --no-group-separator ${WT_NAME} ${UNPAIRED} > ${WT_UNPAIRED}


  if [[ ${KEXTEND_NO_BAIT} == true ]]; then
    NOBAIT="-B"
  fi

  if [ "${INFER_MIN_FREQ_EXTEND}" == true ]; then
    INFER_FREQ_EXTEND="-I"
  fi

  for SM in MT WT; do
    tFILES=${SM}"_FILES"[@]
    tFREQ=${SM}"_MAX_FREQ_OUT"
    for SUBSET in PAIRED; #UNPAIRED; do
    do
      #Using "indirect expansion"
      tFNAME=${SM}"_"${SUBSET} #build variable name
      report "INFO" "Bait and extend ${!tFNAME}" #Get the value
      FILES=()
      for F in ${!tFILES}; do
        FILES+=("-i")
        FILES+=(${F})
      done
      set -o pipefail && ${DIR}/targeted_extend.sh \
      -s ${!tFNAME} -F ${!tFREQ} -k ${KEXTEND_K_BAIT} \
      -m ${KEXTEND_K_MIN} -M ${KEXTEND_K_MAX} -S ${KEXTEND_K_STEP} \
      ${FILES[@]} ${NOBAIT} ${INFER_FREQ_EXTEND} -E ${MEM} \
      -t ${THREADS} ${OV} || exit 1
    done
  done
fi
#((STEP++))
#((STEP++))

#if [ "${STOP_AFTER_STEP}" -eq "${STEP}" ]; then
#  report "INFO" "Stopped at user's request"
#  exit 0
#fi
#((STEP++))
#if [ "${OVERWRITE_FROM_STEP}" -le "${STEP}" ]; then
#  OV="-O"
#fi


#########################################################################################
### XXX
#########################################################################################
#if [ "${SKIP_FIRST_STEPS}" -ge "${STEP}" ]; then
#  report "INFO" "Skipping step ${STEP} at users request"
#else
#  report "INFO" "Running step ${STEP} (...)"
#fi


########################################################################################
## XXX/TODO Verify? Compare to other data? Annotate?
########################################################################################

report "INFO" "Pipeline run completed"






