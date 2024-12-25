FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y


RUN apt install -y build-essential liblapack-dev libblas-dev gfortran autoconf automake cmake curl git libtool libcurl4-openssl-dev pkg-config unzip wget zlib1g-dev python3-dev python3-pip python3-venv python3-wheel
RUN apt-get install -y  --no-install-recommends software-properties-common dirmngr
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
RUN  add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get update
RUN apt install -y --no-install-recommends r-base
RUN R -e "install.packages(c('changepoint', 'data.table', 'randomForest', 'smoother', 'caret', 'BiocManager'))" && \
    R -e "BiocManager::install('rhdf5')"

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    /opt/conda/bin/conda init
    

RUN /opt/conda/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ && \
/opt/conda/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
/opt/conda/bin/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
/opt/conda/bin/conda config --set show_channel_urls yes

RUN mkdir /DEMINERS-main
COPY ./* /DEMINERS-main
RUN R CMD INSTALL /DEMINERS-main/DecodeR_0.1.0.tar.gz

RUN cd /DEMINERS-main && \ 
/opt/conda/bin/conda create -y -n densecall python=3.9  seaborn=0.13.2 && . ~/.bashrc &&  \
conda activate /opt/conda/envs/densecall && \
    pip install --upgrade pip && pip install -r requirements.txt -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple && \
    python setup.py develop && cd ont-bonito-0.7.3 && \
    python setup.py develop

RUN cd / 
WORKDIR /DEMINERS

CMD ["/bin/bash"]