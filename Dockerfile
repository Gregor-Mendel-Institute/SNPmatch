FROM continuumio/miniconda3:latest

RUN conda update -n base -c defaults conda
COPY environment.yml /
RUN mkdir /lustre
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/snpmatch/bin:$PATH
