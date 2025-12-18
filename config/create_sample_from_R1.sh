for fq1 in $(cat SP1.txt); do
  if [[ "$fq1" == *_1.fastq.gz ]]; then
    fq2=${fq1/_1.fastq.gz/_2.fastq.gz}
    prefix=$(basename "$fq1" _1.fastq.gz)
  elif [[ "$fq1" == *_1.fq.gz ]]; then
    fq2=${fq1/_1.fq.gz/_2.fq.gz}
    prefix=$(basename "$fq1" _1.fq.gz)
  else
    echo "Unrecognized file format: $fq1" >&2
    continue
  fi
  echo -e "${prefix}\t${fq1}\t${fq2}"
done
