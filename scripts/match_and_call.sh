#!/bin/bash
SCRIPTSDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

MEM=$(echo "$(grep MemTotal /proc/meminfo | tr -s ' ' ' ' | cut -f2 -d' ') /4000000" | bc) #set max mem to 1/4 physical

THREADS=16

OVERWRITE1=false
OVERWRITE2=false

SUFFIX="matched"

MIN_LEN_UNCLUSTERED=200
MIN_ID1=0.97
MIN_ID2=0.98
MIN_ALN_LEN=1

REVERSE_LEX_ORDER=false;

if [ -z ${COLUMNS} ]; then
  COLUMNS=200
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
usage="USAGE: $(basename $0) [-h] [-i sample_1_kmers_file] [-I sample_2_kmers_file] [-k kmer_size]

where:
  -h Show this helpful help text, and exit
  -i <string>  sample 1 FASTA file name                                       [REQUIRED]
  -I <string>  sample 2 FASTA file name                                       [REQUIRED]
  -k <int>     k-mer size                                                     [REQUIRED]
  -d <float>   min identity (0.0<=id<=1.0) required when clustering/paring seeds (default: ${MIN_ID1})
  -D <float>   min identity (0.0<=id<=1.0) required when calling SNPs (default: ${MIN_ID2})
  -R           report SNPs in reverse lexicographical order of samples labels
  -L <int>     minimum length required to output an unclustered sequence (defaults to ${MIN_LEN_UNCLUSTERED})
  -s <string>  output name suffix (default: ${SUFFIX})
  -C <int>        Print width for some reports, recommended use: -C \${COLUMNS} (defaults to ${COLUMNS})
  -t <int>     number of threads (defaults to ${THREADS})
  -E <int>     max physical memory in GB to be used (defaults to 1/4 of physical mem: ${MEM})
  -O           overwrite existing clustering output files
  -o           overwrite existing SNP-calls output files "

#  -m <int>     min alignment length required when calling SNPs (default: ${MIN_ALN_LEN}) TODO ?

#echo $@

while getopts ":hi:I:s:C:t:d:D:m:L:k:E:ROo" opt; do
  case $opt in
    h) echo "$usage"
       exit;;
    i) SAMPLE1_FILE=${OPTARG};;
    I) SAMPLE2_FILE=${OPTARG};;
    s) SUFFIX=${OPTARG};;
    C) COLUMNS=${OPTARG};;
    d) MIN_ID1=${OPTARG};;
    D) MIN_ID2=${OPTARG};;
    m) MIN_ALN_LEN=${OPTARG};;
    L) MIN_LEN_UNCLUSTERED=${OPTARG};;
    t) THREADS=${OPTARG};;
    k) k=${OPTARG};;
    E) MEM=${OPTARG};;
    R) REVERSE_LEX_ORDER=true;;
    O) OVERWRITE1=true;;
    o) OVERWRITE2=true;;
    ?) printf "Illegal option: '-%s'\n" "${OPTARG}" >&2
       echo "$usage" >&2
       exit 1;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "$usage" >&2
      exit 1;;
  esac
done

#if [ -z ${SAMPLE1_FILE} ] || [ -z ${SAMPLE2_FILE} ] || [ -z ${k} ] ; then
if [ -z ${SAMPLE1_FILE} ] || [ -z ${k} ] ; then
  report "ERROR" "Required option not specified by the user, terminating!!!" >&2
  echo -e "${usage}" >&2
  exit 1
fi

if [ ! -f ${SAMPLE1_FILE} ]; then
  report "ERROR" "Input file ${SAMPLE1_FILE} not found, terminating!" >&2
  exit 1
fi

if [ ! -s ${SAMPLE1_FILE} ]; then
  report "ERROR" "Input file ${SAMPLE1_FILE} is empty, terminating!" >&2
  exit 1
fi

#######################################################################################
#
#######################################################################################

LABEL1=$(head -n 1 ${SAMPLE1_FILE} | tr -d '>' | cut -f1 -d' ' | sed 's/_[0-9]*$//')
LABEL2=$(head -n 1 ${SAMPLE2_FILE} | tr -d '>' | cut -f1 -d' ' | sed 's/_[0-9]*$//')

BNAME1=${SAMPLE1_FILE##*/}
DIR1=${SAMPLE1_FILE%${BNAME1}}

if [ ! -z ${SAMPLE2_FILE} ]; then
  if [ ! -f ${SAMPLE2_FILE} ]; then
    report "ERROR" "Input file ${SAMPLE2_FILE} not found, terminating!" >&2
    exit 1
  fi
  BNAME2=${SAMPLE2_FILE##*/}
  OUTBASE=${BNAME1%%.*}_${BNAME2%%.*}.${SUFFIX}
  CAT=${DIR1}${BNAME1%%.*}_${BNAME2%%.*}.fa
else
  OUTBASE=${BNAME1%%.*}.${SUFFIX}
fi

UNCLUSTERD=${DIR1}${OUTBASE%.${SUFFIX}}.unmatched.fa
CLUSTERSFILE=${DIR1}${OUTBASE}.msa
SNPS_LIST=${CLUSTERSFILE}.snps


#######################################################################################
## XXX VSEARCH,  if output does not exist or overwrite flag -O is set
#######################################################################################
if [[ -f "${CLUSTERSFILE}" ]] && [[ "${OVERWRITE1}" == false ]]; then
  report "WARNING" "${CLUSTERSFILE} already exists, use the -O flag to overwrite " #| tee -a ${LOGFILE}
else
  report "INFO" "Matching extended seeds from ${SAMPLE1_FILE##*/} ${SAMPLE2_FILE##*/}" #| tee -a ${LOGFILE}


    cat ${SAMPLE1_FILE} ${SAMPLE2_FILE} > ${CAT}
    set -o pipefail && vsearch --threads ${THREADS} --strand both --cluster_fast ${CAT} --id ${MIN_ID1} \
    --msaout ${CLUSTERSFILE} --qmask none --iddef 3  --fasta_width 0 --minseqlength ${k} 2>&1 \
    | grep -v "^$" | awk '{ print strftime("%Y-%m-%d %a %H:%M:%S"), "[VSEARCH cluster]  ",$0; fflush(); }' \
    && report "INFO" "Finished VSEARCH clustering"  \
    || (report "ERROR" "Failed  VSEARCH clustering" && exit 1)

#    iddef:
#    0 CD-HIT definition: (matching columns) / (shortest sequence length).
#    1 edit distance: (matching columns) / (alignment length).
#    2 edit distance excluding terminal gaps (same as --id: (matching columns) / (alignment length - terminal gaps)).
#    3 counting each extended gap (internal or terminal) as a single difference: 1.0 - [(mismatches + gaps)/(longest seq len)]
fi


#######################################################################################
## XXX POSTPROCESS clustering output and call SNPs
#######################################################################################
if [[ -f "${SNPS_LIST}" ]] && [[ "${OVERWRITE2}" == false ]]; then
  report "WARNING" "${SNPS_LIST} already exists, use the -o flag to overwrite "
else
  report "INFO" "Identifying variants between clustered sequences in ${CLUSTERSFILE}"

  if [[ "${REVERSE_LEX_ORDER}" == true ]]; then
    REV="--reverse-lex-order "
    rev="-r "
  fi

  set -o pipefail && java  -Xms${MEM}G -Xmx${MEM}G -jar ${YAKAT} vclusters \
  --sample-ids ${LABEL1} ${LABEL2} ${REV} \
  --min-inter-identity ${MIN_ID2} \
  --max-indel-length 99999999999 \
  --clusters-msa ${CLUSTERSFILE} \
  --out-unclustered-fasta ${UNCLUSTERD} \
  --out-unclustered-min-len ${MIN_LEN_UNCLUSTERED} \
  --out-clusters-fasta ${CLUSTERSFILE%.msa}.fa \
  --out-clusters-msa ${CLUSTERSFILE}.fa \
  --print-user-settings \
  --stdout-redirect ${SNPS_LIST} \
  && report "INFO" "Finished parsing VSEARCH clusters" \
  || (report "ERROR" "Failed parsing VSEARCH clusters" && exit 1)


  #XXX Extract FASTA file per sample (with the clustered sequences) for further extending
   set -o pipefail && sed -e 's/^\-*//g' -e 's/\-*$//g' ${CLUSTERSFILE}.fa \
   | tee ${CLUSTERSFILE%.msa}.fa \
   | paste - - \
   | tee >(grep -E ">Cluster_[0-9]+_${LABEL1}.*" | tr '\t' '\n' > ${SAMPLE1_FILE%.fa}.matched.fa) \
   | grep -E ">Cluster_[0-9]+_${LABEL2}.*" | tr '\t' '\n' > ${SAMPLE2_FILE%.fa}.matched.fa \
   && report "INFO" "Separated paired seeds ${SAMPLE1_FILE%.fa}.matched.fa ${SAMPLE2_FILE%.fa}.matched.fa " \
  || (report "ERROR" "Failed Separating paired seeds" && exit 1)


###XXX histogram of alignment lengths
  report "INFO" "Cluster alignment lenghts histogram: "
  set -o pipefail && ./scripts/vcluster_length_distribution.sh ${SNPS_LIST} ${COLUMNS}


    #XXX frequency of numbers of SNPs per CLUSTER
  report "INFO" "Frequencies of SNPs numbers per cluster: "
  set -o pipefail && echo -e "num_snps\tfrequency" \
  && tail -n+2 ${SNPS_LIST} | cut -f1 | sort | uniq -c | awk '{print $1}' | sort -n | uniq -c | awk '{print $2"\t"$1}' \
  | ./scripts/number_format.sh | column -t

  #XXX classes A,B,C based on distance of SNPs from ends
  cat ${SNPS_LIST} | awk -v k=${k} -v o1=${SNPS_LIST}.A -v o2=${SNPS_LIST}.B -v o3=${SNPS_LIST}.C  \
  '{if($1=="ClusterId") {
      print > o1; print > o2; print > o3
    } else if(($2==2*k-1 || $4==2*k-1 || $6==2*k-1) && $7==k && $7==$2-k+1) {
      print > o1
    } else if( $7==k || $7==$2-k+1) {
      print > o2
    } else print >o3
  }' \
  && report "INFO" "Finished allocating SNPs to classes A, B, C based on distance of SNP(s) from sequence ends"  \
  || (report "ERROR" "Failed allocating SNPs to classes A, B, C based on distance of SNP(s) from sequence ends"  && exit 1)

  fgrep -wf <(cut -f1 ${SNPS_LIST} | sort | uniq -d) ${SNPS_LIST} | awk -v k=${k} 'BEGIN{OFS="\t"}; { \
  if($1==prev1) {
    split(previous,prev,"\t");
    if(prev[7]==k && $7==$2-k+1) {
      print previous; print $0
    }
  }
  previous=$0; prev1=$1}' > ${SNPS_LIST}.D \
  && report "INFO" "Finished allocating double-SNPs to class D based on distance of SNP(s) from sequence ends, NOTE: D is a subset of B"  \
  || (report "ERROR" "Failed allocating double-SNPs to class D based on distance of SNP(s) from sequence ends" )


  #PRODUCE SOME STATS
  for f in ${SNPS_LIST}.{A,B,C,D}
  do
    LINES=$(tail -n+2 ${f} | wc -l)
    report "INFO" "${LINES} SNPs in ${f}"
    (echo -e "${LABEL1}\n${LABEL2}" | sort ${rev} | tr '\n' '\t' && echo "Frequency" \
    && ${SCRIPTSDIR}/count_snps_all.sh ${f} ) | column -t
  done

fi