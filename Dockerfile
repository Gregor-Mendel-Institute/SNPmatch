FROM continuumio/miniconda3:latest

RUN conda update -n base -c defaults conda
## installation for having ps -- requirement for Nextflow
RUN apt-get update && apt-get install -y procps 
COPY environment.yml /
RUN mkdir /lustre
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/snpmatch/bin:$PATH
