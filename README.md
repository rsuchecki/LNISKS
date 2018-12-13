
# Table of Contents <!-- omit in toc -->
- [Dependencies](#dependencies)
- [Quick-start example](#quick-start-example)

## Dependencies

* [pigz](https://zlib.net/pigz/) e.g. on debian `sudo apt install pigz`
* [KMC3](https://github.com/refresh-bio/KMC) should be available on PATH or downloaded into `bin/` as follows:
  ```
  mkdir -p bin \
    && cd bin \
    && wget https://github.com/refresh-bio/KMC/releases/download/v3.1.0/KMC3.1.0.linux.tar.gz \
    && tar xzvf KMC3.1.0.linux.tar.gz \
    && rm KMC3.1.0.linux.tar.gz \
    && cd ..
  ```
* [VSEARCH](https://github.com/torognes/vsearch) should be available on PATH or or downloaded into `bin/` as follows:
  ```
  mkdir -p bin \
    && cd bin \
    && wget https://github.com/torognes/vsearch/releases/download/v2.10.2/vsearch-2.10.2-linux-x86_64.tar.gz \
    && tar xzvf vsearch-2.10.2-linux-x86_64.tar.gz \
    && mv vsearch-2.10.2-linux-x86_64/bin/vsearch . \
    && rm -r vsearch-2.10.2-linux-x86_64* \
    && cd ..
  ```
* [`yakat`](https://github.com/rsuchecki/yakat) - `yakat.jar` should be placed in `scripts/`
  ```
  wget --directory-prefix scripts/ https://github.com/rsuchecki/yakat/releases/download/v0.9/yakat.jar
  ```





## Quick-start example

To provide a minimal test set for the pipeline we used [ART read simulator](https://doi.org/10.1093/bioinformatics/btr708) within [RNF framework](https://doi.org/10.1093/bioinformatics/btv524) to generate just over 12,000 paired-end reads from mitochondria assemblies of Arabidopsis and Rice.

The following snippet will run the pipeline, using the two sets of simulated reads. In addition to the required parameters specifying the k-mer size (`-k 16`) the input read sets (`-M, -W`) and labels for these datasets (`-m, -w`), we also specify the number of threads to be used (`-t 2`). In addition we want the minimum frequency of k-mers used for each of the samples to be detected from the respective distributions (`-I -i`). Finally, we use `-C $COLUMNS` to ensure that ascii plots make best use of available terminal width.

```sh
./scripts/lnisks.sh -k 16 \
  -M example/A_thaliana_TAIR10_Mt_ArtIllumina_reads.\?.fq.gz \
  -W example/O_sativa_IRGSP-1.0_Mt_ArtIllumina_reads.\?.fq.gz \
  -m A_thaliana \
  -w O_sativa \
  -t 2 \
  -I -i \
  -C $COLUMNS
```

The output to the terminal is rather verbose, mind you the pipeline was mostly used on large and complex datasets were we found the verbosity useful.
All this information with additional detail on parameters used can be found in the log files in the  output directory, in this case under `output/16-mers/logs/`.