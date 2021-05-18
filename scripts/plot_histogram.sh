#!/bin/bash

#INPUT is a single file - a list of SNPs produced by parsing VSEARCH clustering MSA

if [ -z "${1}" ] || [ -z "${2}" ]; then
  echo Two argunets required,
  echo   filename.histogram
  echo   '${COLUMNS} variable'
  exit 1
fi

awk -v cols=$((${2}-30)) '
  NR==FNR {
    max=$2>max?$2:max;
  };
  NR!=FNR {
    div=max<=cols?1:cols/max;
    x=$2*div;
    ceil=int(x)>=x?int(x):int(x)+1;
    bar=$2>0?gensub(/ /, "#", "g", sprintf("%*s", ceil, "")):" "
    print $1,$2,bar
  }' ${1} ${1}  | ./scripts/number_format.sh | column -t