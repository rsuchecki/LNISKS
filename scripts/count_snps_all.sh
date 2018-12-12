#!/bin/bash
#wrap with something like
#(echo -e "${BNAME1%%.*}\n${BNAME2%%.*}" | sort | tr '\n' '\t' && echo "Frequency" \
#  && ${SCRIPTSDIR}/count_snps_all.sh < ${PAIRS_SNPs} | cut -f1-3) | csvlook -t | tee -a ${LOGFILE}
#OR (to invert column order)
#(echo -ne "wt\tmut\t" && echo "Frequency" && ../../scripts/count_snps.sh < <(awk 'BEGIN{OFS="\t"}; {print $1,$3,$2,$4,$6,$5}' blah.snps) | sed 's/\t$//g' ) | csvlook -t

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

for snp in 'A\t-' 'A\tC' 'A\tG' 'A\tT' 'C\t-' 'C\tG' 'C\tT' '-\tG' 'G\tT' '-\tT'; do
  echo -ne "${snp}\t"
  for f in ${@}; do
    awk '{print $8"\t"$9}' ${f} \
    | sed -e 's/^C\tA/G\tT/g' -e 's/^G\tA/C\tT/g' -e 's/^G\tC/C\tG/g' -e 's/^T\tG/A\tC/g' -e 's/^T\tC/A\tG/g' -e 's/^T\tA/A\tT/g' -e 's/^T\t\-/A\t\-/g' -e 's/^G\t\-/C\t\-/g' -e 's/^\-\tA/\-\tT/g' -e 's/^\-\tC/\-\tG/g' \
    | fgrep -cf <(echo -ne "${snp}")  | tr '\n' '\t'
  done
  echo
done | ${DIR}/number_format.sh | sed 's/\t$//g'
