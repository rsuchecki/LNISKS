#!/bin/bash

#INPUT is a single file - a list of SNPs produced by parsing VSEARCH clustering MSA 

if [ -z "${1}" ] || [ -z "${2}" ]; then
  echo Two argunets required:
  echo   filename.snps
  echo   '${COLUMNS} variable'
  exit 1
fi


MAX_FREQ=$(cut -f2 ${1} | sort -n | uniq -c | sort -nr |  head -1 | awk '{print $1}') 
MAX_LEN=$(cut -f2 ${1} | sort -n | sort -nr |  head -1 | awk '{print $1}') 

#echo ${MAX_FREQ} ${MAX_LEN} 

PAD=$((${#MAX_FREQ}+${#MAX_LEN}+5)) #5=(2*2 base gap)+1 end base


tail -n+2 ${1} | cut -f2 | sort -n | uniq -c | awk -v max=${MAX_FREQ} -v cols=${2} -v pad=${PAD} \
'BEGIN{OFS="\t"; 
  if(max<=cols-pad) {
    div=1
  } else {
    div=(cols-pad)/max
  }; 
  print "len\tfq\thistogram_scaled_to_"(cols-pad)
};
{
  x=$1*div;
  ceil=int(x)>=x?int(x):int(x)+1; 
  print $2,$1,gensub(/ /, "#", "g", sprintf("%*s", ceil, ""))
}' \
| ./scripts/number_format.sh | column -t 


#eg
#len  fq   histogram_scaled_to_189
#103  1    #
#104  1    #
#105  1    #
#107  2    ##
#108  1    #
#109  238  #############################################################################################################################################################################################
#116  1    #
#119  1    #
#121  2    ##
#125  1    #
#133  1    #
#152  1    #
#156  1    #
#159  1    #
#162  1    #
#168  1    #
#169  1    #
#170  1    #
#174  1    #
#183  1    #
#186  1    #
#187  1    #
#191  1    #
#202  1    #
#230  2    ##
#293  1    #

