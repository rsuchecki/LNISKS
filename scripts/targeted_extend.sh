#!/bin/bash
set -eo pipefail #http://redsymbol.net/articles/unofficial-bash-strict-mode/

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

export YAKAT="yakat"
# export PATH="${PATH}:"$(pwd)"/bin"



#THREADS=$(echo "$(nproc) /4" | bc )
#THREADS=$([[ ${THREADS} -ge "4" ]] && echo "${THREADS}" || echo "4") #ensure not less than 4
#MEM=$(echo "$(grep MemTotal /proc/meminfo | tr -s ' ' ' ' | cut -f2 -d' ') /4000000" | bc) #set max mem to 1/4 physical

#printf format for messages
FORMAT="%s %-20s %-10s %10s %s\n"
function report {
  printf "${FORMAT}" "`date +"%Y-%m-%d %a %T"`" "[$(basename $0)]" "[${1}]" "${2}" >&2
}


#echo $0 $*

# exit 1

MEM=$(echo "$(grep MemTotal /proc/meminfo | tr -s ' ' ' ' | cut -f2 -d' ') /4000000" | bc) #set max mem to 1/4 physical
#MEMSTART=$(echo "${MEM} /10" | bc ) #initial mem for JVMm, aiming for 1/10 of max set
#MEMSTART=$([[ ${MEMSTART} -ge "1" ]] && echo "${MEMSTART}" || echo "1") #ensure not less than one GB

#defaults
OVERWRITE1=false
OVERWRITE2=false
K_STEP=10
K_MIN=51
K_MAX=81
INFER_MIN_FREQ=false
MIN_KMER_FREQ=2
MAX_KMER_FREQ=255
BAIT=true

FASTQ_FILES=()

OVERWRITE1=false
OVERWRITE2=false
OVERWRITE3=false

COLUMNS=${COLUMNS:-160}

THREADS=$(echo "$(nproc) /4" | bc )
THREADS=$([[ ${THREADS} -ge "4" ]] && echo "${THREADS}" || echo "4") #ensure not less than 4
##############################
# Parse command line options #
##############################
usage="USAGE: $(basename $0) [-h] [-s input.fasta] [-i input_fastq and/or -d input_fastq_dir] [options]

where:
  -h Show this helpful help text, and exit
  -s <string>     FASTA file with seequences to be used for baiting and then extended   [REQUIRED]
  -d <string>     input FASTQ directory                                    [-d and/or -i REQUIRED]
  -i <string>     input FASTQ files                                        [-d and/or -i REQUIRED]
  -k <int>        baiting k-mer length                                        [REQUIRED unless -B]
  -B              no baiting, k-merize all input reads
  -m <int>        extending k minimum (defaults to ${K_MIN})
  -M <int>        extending k maximum (defaults to ${K_MAX})
  -S <int>        extending k step (defaults to ${K_STEP})
  -I              infer min frequency of a kmer to be used for extending - only makes sense with -B
  -f <int>        min frequency of a kmer to be used for extending (defaults to ${MIN_KMER_FREQ})
  -F <int>        max frequency of a kmer to be used for extending (defaults to ${MAX_KMER_FREQ})
  -t <int>        number of threads for parallelized tasks (defaults to max(4,(nproc/4=${THREADS}))
  -E <int>        max physical memory in GB to be used (defaults to 1/4 of physical mem)
  -r              run KMC in ram-only mode
  -o              overwrite existing matched FASTQ reads, implies -O, -R, ignored if -B
  -O              overwrite existing k-mer DBs, implies -R
  -R              overwrite existing extended seeds
  -C <int>        Print width for some reports, recommended use: -C \$COLUMNS (defaults to ${COLUMNS})
"


while getopts ":hs:d:i:k:Bm:M:S:f:F:t:E:oORIr" opt; do
  case $opt in
    h) echo "$usage"
       exit;;
    s) SEEDS_FASTA=${OPTARG};;
    k) K_BAIT=${OPTARG};;
    B) BAIT=false;;
    d) FASTQ_DIR=${OPTARG};;
    i) FASTQ_FILES+=(${OPTARG});;
    m) K_MIN=${OPTARG};;
    M) K_MAX=${OPTARG};;
    S) K_STEP=${OPTARG};;
    I) INFER_MIN_FREQ=true;;
    f) MIN_KMER_FREQ=${OPTARG};;
    F) MAX_KMER_FREQ=${OPTARG};;
    t) THREADS=${OPTARG};;
    E) MEM=${OPTARG};;
    r) KMC_RAM_ONLY=true;;
    o) OVERWRITE1=true; OVERWRITE2=true; OVERWRITE3=true;;
    O) OVERWRITE2=true; OVERWRITE3=true;;
    R) OVERWRITE3=true;;
    C) COLUMNS=${OPTARG};;
    ?) printf "Illegal option: '-%s'\n" "${OPTARG}" >&2
       echo "$usage" >&2
       exit 1;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "$usage" >&2
      exit 1;;
  esac
done

report "INFO" "Overwrite settigns: (1)=${OVERWRITE1} (2)=${OVERWRITE2} (3)=${OVERWRITE3}"

if [ ${BAIT} == false ] || [ ${INFER_MIN_FREQ} == true ]; then
  report "ERROR" "Inferring k-mer frequency from a set of baited reads is not supported. The -I flag can only be used in conjunction with -B." >&2
  exit 1
fi

if [ -z ${SEEDS_FASTA} ] || ([ -z ${K_BAIT} ] && [ ${BAIT} == true ] ); then
  report "ERROR" "Required option not specified by the user, terminating!!!" >&2
  echo -e "${usage}" >&2
  exit 1
fi

if [ -z ${FASTQ_DIR} ] && [ ${#FASTQ_FILES[@]} -eq 0 ]; then
  report "ERROR" "Must specify -d and/or -i !!!" >&2
  echo -e "${usage}" >&2
  exit 1
fi



for f in ${SEEDS_FASTA} ${FASTQ_FILES[@]} ; do
  if [ ! -f ${f} ]; then
    report "ERROR" "Input file ${f} not found, terminating!" >&2
    exit 1
  fi
done

if [[ ${KMC_RAM_ONLY} == true ]]; then
  KMC_IN_RAM="-r"
fi


#if [[ ${FILTER_NAME} == *['!'@#\$%^\&*()+/.]* ]]; then
#  report "ERROR" "Special characters not allowed in filter name (offending argument: -m ${FILTER_NAME}), terminating!" >&2
#  exit 1
#fi

WD="${SEEDS_FASTA%/*}"
#in case seeds.fa is in pwd
if [ ! -d ${WD} ]; then
  WD="."
fi
INFILE="${SEEDS_FASTA##*/}"
OUTFILE0=${INFILE%.fa}


#if [ "${MAX_LEN}" -eq "${MAX_LEN_DEFAULT}" ]; then
#  OUTFILE=${WD}/${OUTFILE0%.fasta}.GE_${MIN_LEN}bp.unpaired.fa
#else
#  OUTFILE=${WD}/${OUTFILE0%.fasta}.further-extended.fa #GE_${MIN_LEN}bp.LE_${MAX_LEN}bp.unpaired.fa
  mkdir -p ${WD}/DBs
  OUTBASE=${WD}/DBs/${OUTFILE0%.fasta}.catch_k
  OUTFILE=${WD}/${OUTFILE0%.fasta}.catch_k${K_BAIT}.k_${K_MIN}_${K_STEP}_${K_MAX}.fasta
#fi
#mkdir -p ${WD}/logs
#LOGFILE="${WD}/logs/${INFILE}_get_long_unpaired"$(date +"%Y-%m-%d-%a-%Hh%Mm%Ss").log
#exec > >(tee -a ${LOGFILE})
#exec 2> >(tee -a ${LOGFILE} >&2)


CATCH_SE=${OUTBASE}${K_BAIT}_SE.fastq.gz
CATCH_R1=${OUTBASE}${K_BAIT}_R1.fastq.gz
CATCH_R2=${OUTBASE}${K_BAIT}_R2.fastq.gz


report "INFO" "Targeted extending initiated, k=({k_i, K_{i-1}+${K_STEP},...,k_j} | ${K_MIN}<=k<=${K_MAX}), MIN_KMER_FREQ=${MIN_KMER_FREQ}, MAX_KMER_FREQ=${MAX_KMER_FREQ}"

if [[ ${BAIT} == true ]]; then
  report "INFO" "Baiting reads, k=${K_BAIT}"




  #######################################################################################
  #XXX  BAIT READS
  #######################################################################################

# BAIT SE
  if [ ! -f ${CATCH_SE} ] || [ ${OVERWRITE1} == true ]; then
    SE=$((ls ${FASTQ_FILES[@]} | fgrep -i SE.fastq.gz | tr '\n' ' ') || echo '' )
    if [[ ! -z ${FASTQ_DIR} ]]; then
      FASTQ_DIR="${FASTQ_DIR}/*"
      SEd=$(readlink -f ${FASTQ_DIR} | grep -iE f(ast)?q.gz | grep -viE '(1|2)\.f(ast)?q\.gz' | tr '\n' ' ')
      report "INFO" "INFILES: ${SEd}"
    # else
    #   report "INFO" "No FASTQ dir?"
    fi
    if [[ ! -z ${SE} ]] || [[ ! -z ${SEd} ]]; then
      set -o pipefail && pigz -dcp2 ${SE} ${SEd} | paste - - - - \
      |  ${YAKAT} kmatch --JVM "-Xms${MEM}G -Xmx${MEM}G" \
        --k-mer-length ${K_BAIT} --k-mers <(tr -d '-' < ${SEEDS_FASTA}) \
        --threads ${THREADS} --print-user-settings \
      | tr '\t' '\n' | pigz -9cp8 > ${CATCH_SE} \
      || exit 1
    # else
    #   report "INFO" "No FASTQ dir OR FASTQ files?"
    fi
  fi
# BAIT PE
  if [ ! -f ${CATCH_R1} ] || [ ! -f ${CATCH_R2} ] || [ ${OVERWRITE1} == true ]; then
    R1=$((ls ${FASTQ_FILES[@]} | grep -iE '1\.f(ast)?q\.gz' 2> /dev/null | tr '\n' ' ') || echo '')
    R2=$((ls ${FASTQ_FILES[@]} | grep -iE '2\.f(ast)?q\.gz' 2> /dev/null | tr '\n' ' ') || echo '')
    if [[ ! -z ${FASTQ_DIR} ]]; then
      FASTQ_DIR="${FASTQ_DIR}/*"
      R1d=$(readlink -f ${FASTQ_DIR} | grep -iE '1\.f(ast)?q\.gz' 2> /dev/null | tr '\n' ' ')
      R2d=$(readlink -f ${FASTQ_DIR} | grep -iE '1\.f(ast)?q\.gz' 2> /dev/null | tr '\n' ' ')
    fi
    if [[ ( ! -z ${R1}  &&  ! -z ${R2} ) || ( ! -z ${R1d}  &&  ! -z ${R2d} ) ]]; then
      set -o pipefail && paste \
        <(pigz -dcp2 ${R1} ${R1d} | paste - - - -) \
        <(pigz -dcp2 ${R2} ${R2d} | paste - - - -) \
      |  ${YAKAT} kmatch --JVM "-Xms${MEM}G -Xmx${MEM}G" \
        --k-mer-length ${K_BAIT} --k-mers <(tr -d '-' < ${SEEDS_FASTA})  \
        --threads ${THREADS} --print-user-settings \
      | tee >(cut -f 1-4 -d$'\t' | tr '\t' '\n' | pigz -9cp8 > ${CATCH_R1}) \
      | cut -f 5-8 -d$'\t' | tr '\t' '\n' | pigz -9cp8 > ${CATCH_R2} \
      || exit 1
    fi
  else
    report "WARNING" "One or more of the matched FASTQ files already exist(s), use the -o flag to overwrite " | tee -a ${LOGFILE}
  fi
  if [ -f ${CATCH_R1} ]; then
    TO_KMERIZE+="-i ${CATCH_R1} "
  fi
  if [ -f ${CATCH_R2} ]; then
    TO_KMERIZE+="-i ${CATCH_R2} "
  fi
  if [ -f ${CATCH_SE} ]; then
    TO_KMERIZE+="-i ${CATCH_SE} "
  fi

else # NO BAITING, JUST KMERIZE INPUT READS
  TO_KMERIZE=()
  for f in $(ls ${FASTQ_FILES[@]}); do
    TO_KMERIZE+=("-i" "${f}")
  done
  if [[ ! -z ${FASTQ_DIR} ]]; then
    for F in ${FASTQ_DIR}/*; do
      TO_KMERIZE+=("-i" "${F}")
    done
  fi
fi


#######################################################################################
#XXX  k-merize for each k
#######################################################################################

for k in $(seq ${K_MIN} ${K_STEP} ${K_MAX} ); do
  FILE=${WD}/DBs/${k}-mers_${OUTFILE0%.fasta}.catch_k_${K_BAIT}.db
#  echo $(ls ${FILE}.kmc_pre)
#  echo "OV2 ${OVERWRITE2} ${FILE}"

  if [ ! -f ${FILE}.kmc_pre ] || [ ${OVERWRITE2} == true ]; then
    report "INFO" "Counting ${k}-mers in reads..." ${TO_KMERIZE[@]}
    set -o pipefail && ${DIR}/count_kmers.sh -k ${k} -b ${OUTBASE#${WD}/DBs/}_${K_BAIT} -d ${WD}/DBs -m ${MAX_KMER_FREQ} \
    -t ${THREADS} -S ${MEM} -L ${MIN_KMER_FREQ} -U ${MAX_KMER_FREQ} -O ${TO_KMERIZE[@]} ${KMC_IN_RAM} \
    || exit 1
  else
    report "WARNING" "${FILE} already exists, use the -O flag to overwrite " | tee -a ${LOGFILE}
  fi
done
report "INFO" "Finished counting k-mers in reads"


if [ ! -f ${OUTFILE} ] || [ ${OVERWRITE3} == true ]; then
  #######################################################################################
  #XXX  KEXTEND
  #######################################################################################
  # KEXTEND="java -Xms${MEM}G -Xmx${MEM}G -Xss10M -jar ${YAKAT} kextend"

  # DUMP=()

# echo "$0 $@"
  for k in $(seq ${K_MIN} ${K_STEP} ${K_MAX} ); do
    DB=${WD}/DBs/${k}-mers_${OUTFILE0%.fasta}.catch_k_${K_BAIT}.db
    if [ "${INFER_MIN_FREQ}" == true ]; then
      HISTO=${DB}.histogram
      report "INFO" "Attempting to estimate the minimum frequency from ${DB}..."
      kmc_tools -hp transform ${DB} histogram ${HISTO} || report "ERROR" "Failed to estimate the minimum frequency for ${DB} k-mers to be included"
      MIN_KMER_FREQ=$(awk -v prev=99999999999 '{if ($2>prev){print $1-1; exit}; prev=$2}' ${HISTO})
      # Each of the following two lines causes slient exit and/or loop break
      # MAX_PRINT=$(awk -v MIN="${MIN_KMER_FREQ}" '$1>=MIN' "${HISTO}" | sort -k2,2nr | head -1 | cut -f1)
      # ${DIR}/plot_histogram.sh ${HISTO} ${COLUMNS} | tee ${HISTO}.plot && report true #| awk -vM="${MAX_PRINT}" 'NR<=(M*2)' >&2
      report "INFO" "Estimated minimum frequency of k-mers to be included: ${MIN_KMER_FREQ}, see ${HISTO}.plot"
    fi
    report "INFO" "Dumping ${k}-mers from caugth reads > kextend "
    # DUMP+=("<(kmc_dump" "-ci${MIN_KMER_FREQ}" "-cx${MAX_KMER_FREQ}" "${DB}" "/dev/stdout)")
    ionice -c 2 -n 7 kmc_dump -ci${MIN_KMER_FREQ} -cx${MAX_KMER_FREQ} ${DB} /dev/stdout  #| pigz -9cp8 > ${DB}.dump.gz
  done \
  | ${YAKAT} kextend --JVM "-Xms${MEM}G -Xmx${MEM}G -Xss10M" \
      --seed-file <(tr -d '-' < ${SEEDS_FASTA}) --threads ${THREADS} --fasta-out > ${OUTFILE}
  report "INFO" "Finished extending seeds:  ${OUTFILE}"
else
  report "WARNING" "${OUTFILE} already exists, use the -R flag to overwrite "
fi