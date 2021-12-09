################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.10.3

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.10.3"
LABEL software="rna_splicing-nf"
LABEL software.version="0.1"
LABEL about.summary="Container image containing all requirements for rna_slicing-nf"
LABEL about.home="https://github.com/neonicoo/RNA-splicing-workflow"
LABEL about.documentation="https://github.com/neonicoo/RNA-splicing-workflow/README.md"
LABEL about.license_file="https://github.com/neonicoo/RNA-splicing-workflow/LICENSE.txt"
LABEL about.license="MIT License"
LABEL authors="Son-Michel DINH, Aurélie FISCHER, Maëlys MARRY, Nicolas MENDIBOURE"

################## MAINTAINER ######################
MAINTAINER neonicoo <**nicolas.mendiboure@insa-lyon.fr**>

################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda config --set channel_priority strict
RUN conda env create -n rna-splicing -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/rnaseq-nf/bin:$PATH
