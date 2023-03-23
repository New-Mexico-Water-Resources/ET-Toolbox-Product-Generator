FROM condaforge/mambaforge

ENV APP_ROOT /app

# update Ubuntu package manager
RUN apt-get update
RUN apt-get -y install gcc make
# install fish shell
RUN apt-add-repository ppa:fish-shell/release-3; apt-get -y install fish; chsh -s /usr/local/bin/fish; mamba init fish

RUN mamba install -v -y  "python=3.9"

RUN mamba install -v -y astropy \
    beautifulsoup4 \
    cmake \
    descartes \
    ephem \
    geopandas \
    h5py \
    imageio \
    imageio-ffmpeg \
    jupyter \
    keras \
    matplotlib \
    netcdf4 \
    nose \
    pip \
    pygeos \
    pygrib \
    pyhdf \
    pyresample \
    pysolar \
    pystac-client \
    python-dateutil \
    requests \
    scikit-image \
    tensorflow \
    termcolor \
    untangle \
    urllib3 \
    xmltodict \
    xtensor \
    xtensor-python \
    wget

RUN mamba install -v -y  "gdal>3.1" "rasterio>1.0.0" "setuptools<58" "shapely<2.0.0"

RUN pip install \
    mgrs \
    pybind11_cmake \
    pycksum \
    pygrib \
    sentinelsat \
    spacetrack

# install app
RUN mkdir ${APP_ROOT}
WORKDIR ${APP_ROOT}
ADD . ${APP_ROOT}

RUN make install-docker
