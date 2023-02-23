FROM condaforge/mambaforge

ENV APP_ROOT /app

# update Ubuntu package manager
RUN apt-get update
# install fish shell
RUN apt-add-repository ppa:fish-shell/release-3; apt-get -y install fish; chsh -s /usr/local/bin/fish; mamba init fish
# install dependencies
RUN mamba install -v -y \ 
    astropy \
    beautifulsoup4 \
    cmake \
    descartes \
    ephem \
    "gdal>3.1" \
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
    pybind11 \
    pygeos \
    pyhdf \
    pyresample \
    pysolar \
    pystac-client \
    "python=3.9" \
    python-dateutil \
    "rasterio>1.0.0" \
    requests \
    scikit-image \
    "setuptools<58" \
    "shapely<2.0.0" \
    tensorflow \
    termcolor \
    untangle \
    urllib3 \
    xmltodict \
    xtensor \
    xtensor-python \
    wget
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

RUN python setup.py install
