
FROM debian:bookworm

ENV DEBIAN_FRONTEND=noninteractive


#RUN apt install apt-transport-https ca-certificates
RUN echo "deb http://mirrors.ustc.edu.cn/debian bookworm main contrib non-free non-free-firmware" | tee -a /etc/apt/sources.list
RUN echo "deb http://mirrors.ustc.edu.cn/debian bookworm-updates main contrib non-free non-free-firmware" | tee -a /etc/apt/sources.list
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y build-essential
RUN apt-get install -y liblapack-dev
RUN apt-get install -y libblas-dev
RUN apt-get install -y gfortran
RUN apt-get install -y autoconf
RUN apt-get install -y automake
RUN apt-get install -y cmake
RUN apt-get install -y curl
RUN apt-get install -y git
RUN apt-get install -y libtool
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y pkg-config
RUN apt-get install -y unzip
RUN apt-get install -y wget
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y python3-dev
RUN apt-get install -y python3-pip
RUN apt-get install -y python3-venv
RUN apt-get install -y python3-wheel
RUN apt-get install -y --no-install-recommends software-properties-common
RUN apt-get install -y --no-install-recommends dirmngr
RUN apt-get install -y --no-install-recommends wget autoconf


RUN apt update

RUN apt-get install -y --no-install-recommends  r-base

RUN R -e "install.packages('changepoint', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "install.packages('data.table', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "install.packages('randomForest', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "install.packages('smoother', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "install.packages('caret', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "install.packages('BiocManager', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "BiocManager::install('rhdf5')"


RUN UNAME=$(uname -m) && \
    case $UNAME in \
        aarch64*) \
            URL="https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-aarch64.sh" ; \
        ;; \
        x86_64*) \
            URL="https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh" ; \
        ;; \
        *) \
            echo "Unsupported architecture: $UNAME" && exit 1 ; \
        ;; \
    esac && \
    wget $URL && \
    bash Miniconda3-latest-Linux-${UNAME%}.sh -b -p /opt/conda && \
    /opt/conda/bin/conda init

RUN /opt/conda/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
    /opt/conda/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
    /opt/conda/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
    /opt/conda/bin/conda config --set show_channel_urls yes


RUN mkdir /DEMINERS-main
COPY ./* /DEMINERS-main
RUN R CMD INSTALL /DEMINERS-main/DecodeR_0.1.0.tar.gz

# RUN cd /DEMINERS-main && \
#     /opt/conda/bin/conda create -y -n densecall && \
#     /opt/conda/bin/conda init &&\
#     . /root/.bashrc && \
#     conda activate densecall

RUN cd /DEMINERS-main && \
    /opt/conda/bin/conda create -y -n densecall python=3.9 seaborn=0.13.2 && \
    /opt/conda/bin/conda init &&\
    . /root/.bashrc && \
    conda activate densecall && \
    pip install --upgrade pip && \
    pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    python setup.py develop && \
    cd ont-bonito-0.7.3 && \
    python setup.py develop


WORKDIR /DEMINERS


CMD ["/bin/bash"]