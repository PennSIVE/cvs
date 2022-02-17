FROM mimosa_base

# R dependencies
RUN echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /usr/local/lib/R/etc/Rprofile.site && \
    Rscript -e "install.packages(c('devtools', 'rlist', 'xml', 'argparser', 'fslr'))" && \
    Rscript -e "source('https://neuroconductor.org/neurocLite.R'); neuro_install(c('ANTsRCore', 'ANTsR', 'extrantsr', 'NiftiArray', 'oasis', 'mimosa', 'rtapas'), release = 'stable')"
# install MASS
ENV PATH=/opt/dramms-1.5.1/bin:/opt/mass-1.1.1/bin:/dramms-1.5.1-source/build/bundle/bin:${PATH}
RUN apt-get update && apt-get install -y cmake-curses-gui build-essential zlib1g-dev git git-lfs && \
    curl -LO https://github.com/ouyangming/DRAMMS/raw/master/dramms-1.5.1-source.tar.gz && \
    tar -xzf dramms-1.5.1-source.tar.gz && rm -rf dramms-1.5.1-source.tar.gz && \
    cd dramms-1.5.1-source/build && cmake -D CMAKE_INSTALL_PREFIX:STRING=/opt/dramms-1.5.1 . && make && \
    cd / && git clone https://github.com/CBICA/MASS.git && \
    cd MASS && git checkout 1.1.1 && git lfs fetch --all && \
    mkdir mass-1.1.1-build && cd mass-1.1.1-build && \
    cmake -D CMAKE_INSTALL_PREFIX:STRING=/opt/mass-1.1.1 -D SCHEDULER:STRING=NONE .. && \
    make && make install && chmod -R 775 /opt/mass-1.1.1 && \
    cd / && rm -rf MASS && \
    apt-get --purge -y autoremove git git-lfs curl
# python dependencies + bids validator
COPY requirements.txt .
RUN apt-get update && \
    apt-get install -y python3 python3-pip nodejs npm && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    npm install -g bids-validator && \
    python3 -m pip install -r requirements.txt


ENV PYTHONPATH=""
ENV ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

COPY . /

ENTRYPOINT ["/run.py"]

