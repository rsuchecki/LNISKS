
## Dependencies

* [pigz](https://zlib.net/pigz/)
* [KMC](https://github.com/refresh-bio/KMC)
* [VSEARCH](https://github.com/torognes/vsearch)
* [`yakat`](https://github.com/rsuchecki/yakat) (`yakat.jar` should be placed in `scripts/`)


## Quick-start example

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