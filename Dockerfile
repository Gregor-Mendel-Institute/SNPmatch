FROM python:2

WORKDIR /usr/src/app

RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN mkdir /lustre
RUN pip install numpy
RUN pip install --no-cache-dir git+https://github.com/Gregor-Mendel-Institute/SNPmatch.git@3.0.1

