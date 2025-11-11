FROM continuumio/miniconda3

WORKDIR /app

# Create the environment:
COPY cellPool_py3.yml .
COPY cellPool_cellpose.yml .
COPY cellPool_R.yml .
COPY retrieveCellposeURLs.py .
COPY retrieveCellposeURLs.sh .
RUN apt-get -y update
RUN apt-get install -y ffmpeg
RUN apt-get install -y wget
RUN chmod +x ./retrieveCellposeURLs.sh
#RUN conda env create -f cellPool_R.yml
RUN conda env create -f cellPool_cellpose.yml
RUN conda env create -f cellPool_py3.yml
SHELL ["/bin/bash", "--login", "-c"]
RUN exec ./retrieveCellposeURLs.sh
RUN wget https://github.com/opensourcedesign/fonts/raw/master/ubuntu-font-family-0.80/UbuntuMono-R.ttf
RUN wget https://github.com/marrlab/BaSiC/archive/refs/heads/master.zip
RUN apt install unzip
RUN unzip /app/master.zip
RUN rm /app/master.zip
RUN mkdir /app/cellpool
RUN mkdir /app/cellpool_procedures
COPY cellpool /app/cellpool/
COPY cellpool_procedures /app/cellpool_procedures/
RUN chmod -R +rw /app/cellpool/
RUN chmod -R +rw /app/cellpool_procedures/
