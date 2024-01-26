FROM openjdk:11.0-jre

RUN apt-get update \
    && apt-get install -y \
      pigz \
      gawk \
      bc \
      bsdmainutils \
    && rm -rf /var/lib/apt/lists/*

ENV KMC_VER=3.2.1
ENV VSEARCH_VER=2.17.0
ENV YAKAT_VER=0.9.5

WORKDIR /usr/bin

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends procps \
  && apt-get -qq -y autoremove \
  && apt-get autoclean \
  && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
  && wget https://github.com/refresh-bio/KMC/releases/download/v${KMC_VER}/KMC${KMC_VER}.linux.tar.gz \
  && tar xzvf KMC${KMC_VER}.linux.tar.gz --strip-components=1 \
  && rm KMC${KMC_VER}.linux.tar.gz \
  && wget https://github.com/torognes/vsearch/releases/download/v${VSEARCH_VER}/vsearch-${VSEARCH_VER}-linux-x86_64.tar.gz \
  && tar xzvf vsearch-${VSEARCH_VER}-linux-x86_64.tar.gz \
  && mv vsearch-${VSEARCH_VER}-linux-x86_64/bin/vsearch . \
  && rm -r vsearch-${VSEARCH_VER}-linux-x86_64* \
  && wget https://github.com/rsuchecki/yakat/releases/download/v${YAKAT_VER}/yakat \
  && chmod +x yakat


#WORKDIR /LNISKS

COPY scripts .

#COPY example /example

WORKDIR /

#ENTRYPOINT ["lnisks.sh"]
#CMD ["-h"]
#CMD ["lnisks.sh"]

#docker run -v "$PWD":"$PWD" -w "$PWD"  rsuchecki/lnisks:latest lnisks.sh -k 16 -M example/A_thaliana_TAIR10_Mt_ArtIllumina_reads.\?.fq.gz   -W example/O_sativa_IRGSP-1.0_Mt_ArtIllumina_reads.\?.fq.gz   -m A_thaliana   -w O_sativa   -t 2   -I -i   -C $COLUMNS
