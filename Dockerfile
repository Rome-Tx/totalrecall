FROM debian:bullseye as source

# install required packages from repositories
#RUN sed -i -e "s/ main[[:space:]]*\$/ main contrib non-free/" /etc/apt/sources.list
RUN DEBIAN_FRONTEND=noninteractive apt-get update -y && apt-get upgrade -y && apt-get install -y --no-install-recommends python3 python3-biopython python3-pysam python3-scipy python3-sklearn libhts-dev cython3 gawk sed procps psmisc util-linux python3-pip python3-dev gcc snakemake xvfb xauth default-jre ncbi-blast+ last-align fml-asm zstd

# add the code, the genome file and the IGV files
RUN mkdir -p /opt/totalrecall/last && mkdir /code
ADD files/igv_211.tgz /opt
COPY --chmod=755 files/samtools /usr/local/bin/
ADD totalrecall.src.tgz /code

# install/setup data files
RUN cd /code && pip3 install . && \
    #install -m 755 contrib/cromwell_monitoring_script.sh /usr/bin && \
    mv Snakefile /opt/totalrecall && mv snakemake /opt/totalrecall && \
    mv classifier/* /opt/totalrecall && \
    mv transposon_data/hs/* /opt/totalrecall && cd /opt/totalrecall/last && \
    lastdb -R 10 transposon ../sequence.fa && \
    lastdb -c -R 10 signatures ../signatures.fa && \
    lastdb -c -R 10 first_pass ../signatures_1p.fa && \
    rm -rf /code

# place the variable part in the very end
ARG genome_tar="files/hg38.tgz"
ADD "$genome_tar" /opt/totalrecall

# copy all the files over to the final image in order to squash it
FROM scratch
COPY --from=source / /

# environment variables
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV LC_COLLATE C
ENV TERM xterm

# default directory
WORKDIR /

# default interactive shell
CMD /bin/bash
