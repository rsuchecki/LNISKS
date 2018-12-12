#!/bin/bash
while read line
do
  for word in ${line}; do 
    if [[ "${word}" =~ ^-?[0-9]+$ ]]; then
      printf "%'d\t" ${word}; 
    else
      printf "%s\t" ${word};
    fi
  done
  echo
done < "${1:-/dev/stdin}"


