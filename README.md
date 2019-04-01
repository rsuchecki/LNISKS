
# Table of Contents <!-- omit in toc -->
- [About LNISKS](#about-lnisks)
  - [How to cite](#how-to-cite)
- [Dependencies](#dependencies)
- [Quick-start](#quick-start)
  - [Example data](#example-data)
  - [Example run](#example-run)
- [Execution on real data](#execution-on-real-data)
  - [The choice of k-mer size](#the-choice-of-k-mer-size)
  - [Call prioritization](#call-prioritization)

## About LNISKS

LNISKS (longer needle in a scanter _k_-stack) is a high-throughput pipeline developed for reference free mutation identification.
LNISKS extends concepts from [NIKS (needle in a _k_-stack)](dx.doi.org/10.1038/nbt.2515) but does not share code with it.
It enables quick and reliable identification of mutations in large and complex crop genomic datasets without a reference genome assembly.

Through the use of state-of-the-art tools such as [KMC3](https://github.com/refresh-bio/KMC) and [VSEARCH](https://github.com/torognes/vsearch) as well as algorithms implemented in [`yakat`](https://github.com/rsuchecki/yakat) toolkit, LNISKS can process datasets consisting of billions of reads within hours.

For example, given hexaploid wheat datasets of 19X and 23X for wild-type and mutant bulks, respectively, the pipeline completed in 2 hours and 33 minutes on an allocation od of 32G of RAM and 16 logical cores (Xeon E5-2699 v3 2.30GHz CPUs).
For comparison, sequential, single-threaded reading, decompression and discarding of the same input data takes 2 hours.

### How to cite

**LNISKS: Reference-free mutation identification for large and complex crop genome**. *Radosław Suchecki, Ajay Sandhu, Stéphane Deschamps, Victor Llaca, Petra Wolters, Nathan S. Watson-Haigh, Margaret Pallotta, Ryan Whitford, Ute Baumann*
bioRxiv 580829; doi: https://doi.org/10.1101/580829

## Dependencies

After you download the latest release or clone this repository go to LNISKS directory, and install/download the following:

* [pigz](https://zlib.net/pigz/), e.g. on debian
  ```
  sudo apt install pigz
  ```
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
  wget --directory-prefix scripts/ https://github.com/rsuchecki/yakat/releases/download/v0.9.1/yakat.jar
  ```


## Quick-start

### Example data

To provide a minimal test set for the pipeline we used [ART read simulator](https://doi.org/10.1093/bioinformatics/btr708) within [RNF framework](https://doi.org/10.1093/bioinformatics/btv524) to generate just over 12,000 paired-end reads from mitochondria assemblies of Arabidopsis and Rice.


### Example run

The following snippet will run the pipeline, using the two sets of simulated reads.
In addition to the required parameters specifying the k-mer size (`-k 16`) the input read sets (`-M, -W`) and labels for these datasets (`-m, -w`), we also specify the number of CPU threads to be used (`-t 2` but the more the merrier).
We want the minimum frequency of *k*-mers used for each of the samples to be detected from the respective distributions (`-I -i`).
Finally, we use `-C $COLUMNS` to ensure that ascii plots make best use of the available terminal width.

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

The output to the terminal is rather verbose, but the pipeline was mostly used on large and complex datasets were we found this verbosity useful.
All this information with additional detail on parameters used can be found in the log files in the  output directory, in this case under `output/16-mers/logs/`.

## Execution on real data

The intended application of LNISKS is to Bulked Segregant data, but in principle it should work well for identifying homozygous mutations between two sets of reads.
The example run above is a good indication on how you could run it on real data, however

* be wary of `-i` and  `-I` settings as the minimum frequency estimation is rather basic and likely only applicable to whole genome sequencing data
* select a higher *k*, 16 just happens to work reasonably well for the tiny example data sets

### The choice of k-mer size

Longer *k* -mers are more likely to be unique within a genome than shorter *k*-mers, but as *k* increases (up to the read length), so does the sequencing coverage required for *k*-mers to occur with sufficient frequency to be distinguishable from low frequency *k*-mers containing sequencing errors.

We want *k* to be close to highest possible value for which neither of the distributions for the two input sets is truncated.
We start by running the pipeline for several values of `-k`, but to speed things up we can stop before identification of sample-specific k-mers by applying `-S 2`.
We then investigate the generated *k*-mer histograms.

Using our example dataset, we run the pipeline with `-k 10` and `-k 20` and observe that at `k=10` the number of k-mers occurring 2,3,4,... times goes down before going up again, while at `k=20` the initial slope is missing.
The distribution at `k=20` is truncated so the sequencing coverage is too low for this k-mer size and some k-mers representing real sequence have been lost.
```
./scripts/plot_histogram.sh output/20-mers/20-mers_A_thaliana.histogram 80 | head -15
2    39,370  ################################
3    55,548  #############################################
4    62,803  ##################################################
5    57,318  ##############################################
6    43,635  ###################################
7    29,447  ########################
8    17,304  ##############
9    10,156  #########
10   5,651   #####
11   3,189   ###
12   1,722   ##
13   1,049   #
14   709     #
15   439     #
16   279     #
```

```
./scripts/plot_histogram.sh output/10-mers/10-mers_A_thaliana.histogram 80 | head -15
2    18,706  ########################################
3    14,591  ###############################
4    19,342  #########################################
5    22,533  ################################################
6    23,806  ##################################################
7    21,528  ##############################################
8    17,861  ######################################
9    14,465  ###############################
10   11,964  ##########################
11   9,939   #####################
12   8,478   ##################
13   7,052   ###############
14   6,292   ##############
15   5,366   ############
16   4,657   ##########
```
Note that by default we exclude all *k*-mers occurring just once as these are likely to arise from sequencing errors in the input reads.

Also note that for approx 20X data sets from common wheat the optimal *k*-mer size was 54.

### Call prioritization

One way to prioritize variant calls is to look in more detail into how well supported they are. For example, let's take the highest confidence, category A calls, although the full set of calls or each remaining category (B, C, D) could be processed in the same way.

```bash
CALLS=output/16-mers/16-mers_O_sativa_MINUS_A_thaliana_extended_16-mers_A_thaliana_MINUS_O_sativa_extended.matched.msa.snps.A
```

We then take the k-mer databases for each of the input data sets (A_thaliana, O_sativa), give them short labels (At and Os), and dump all k-mers into `yakat snpmers`.

```bash
(
  echo Os; bin/kmc_dump output/16-mers/16-mers_O_sativa.db /dev/stdout; \
  echo At; bin/kmc_dump output/16-mers/16-mers_A_thaliana.db /dev/stdout; \
) | java -jar scripts/yakat.jar snpmers \
    --k-mer-length 16 \
    --lnisks-snps ${CALLS} \
    --print-user-settings \
    --out-calls /dev/null \
    --out-calls-with-evidence output/16-mers/snpmers.A.tsv
```

In a real word scenario, each bulk would likely comprise of multiple individuals which could be input by looping through individual k-mer databases. For each individual, its name/id needs to be printed immediately before the k-mers are dumped from the corresponding database - this would first require to compute individual k-mer database using `kmc` either directly or via the [`count_kmers.sh` wrapper](https://github.com/rsuchecki/LNISKS/blob/master/scripts/count_kmers.sh) as used by LNISKS pipeline.


The head of output table generated by the example above should look like this:


| Ref:Pos        | Al1 | Cov1 | Al2 | Cov2 | Os_call | Os_cov1 | Os_cov2 | Os_freq1 | Os_freq2 | At_call | At_cov1 | At_cov2 | At_freq1 | At_freq2 |
|----------------|-----|------|-----|------|---------|---------|---------|----------|----------|---------|---------|---------|----------|----------|
| Cluster_528:16 | G   | 16   | A   | 16   | G       | 16      | 0       | 7        | 0        | A       | 0       | 16      | 0        | 6        |
| Cluster_530:16 | T   | 16   | C   | 16   | T       | 16      | 0       | 6        | 0        | C       | 0       | 16      | 0        | 5        |
| Cluster_531:16 | C   | 16   | A   | 16   | C       | 16      | 0       | 6        | 0        | A       | 0       | 16      | 0        | 4        |
| Cluster_532:16 | A   | 16   | C   | 16   | A       | 16      | 0       | 3        | 0        | C       | 0       | 16      | 0        | 6        |
| Cluster_533:16 | A   | 16   | C   | 16   | A       | 16      | 0       | 5        | 0        | C       | 0       | 16      | 0        | 6        |
| Cluster_534:16 | G   | 16   | C   | 16   | G       | 16      | 0       | 6        | 0        | C       | 0       | 16      | 0        | 3.5      |
| Cluster_535:16 | T   | 16   | C   | 16   | T       | 16      | 0       | 5        | 0        | C       | 0       | 16      | 0        | 6        |
| Cluster_536:16 | C   | 16   | T   | 16   | C       | 16      | 0       | 4        | 0        | T       | 0       | 16      | 0        | 7        |
| Cluster_537:16 | T   | 16   | A   | 16   | T       | 16      | 0       | 3        | 0        | A       | 0       | 16      | 0        | 7        |

Column headers can be deciphered as follows:

* From each LNISKS call we get:
  * Al1 / Al2 - called bases/alleles
  * Cov1 / Cov2 - total number of possible k-mers supporting respective base/allele
* From each set of k-mers (here labelled {At, OS}) we obtain:
  * _call - base call based on all input k-mers
  * _cov - number of k-mers supporting a base/allele
  * _freq - median frequency o k-mers supporting that base/allele

These can be used for filtering or sorting calls.