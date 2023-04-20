FROM condaforge/mambaforge

ENV APP_ROOT /app

# update Ubuntu package manager
RUN apt-get update
RUN apt-get -y install gcc
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

#RUN mamba install -v -y  "python=3.9"
#
#RUN mamba install -v -y astropy
#RUN mamba install -v -y beautifulsoup4
#RUN mamba install -v -y cmake
#RUN mamba install -v -y descartes
#RUN mamba install -v -y ephem
#RUN mamba install -v -y geopandas
#RUN mamba install -v -y h5py
#RUN mamba install -v -y imageio
#RUN mamba install -v -y imageio-ffmpeg
#RUN mamba install -v -y jupyter
#RUN mamba install -v -y keras
#RUN mamba install -v -y matplotlib
#RUN mamba install -v -y netcdf4
#RUN mamba install -v -y nose
#RUN mamba install -v -y pip
#RUN mamba install -v -y pygeos
#RUN mamba install -v -y pygrib
#RUN mamba install -v -y pyhdf
#RUN mamba install -v -y pyresample
#RUN mamba install -v -y pysolar
#RUN mamba install -v -y pystac-client
#RUN mamba install -v -y python-dateutil
#RUN mamba install -v -y requests
#RUN mamba install -v -y scikit-image
#RUN mamba install -v -y tensorflow
#RUN mamba install -v -y termcolor
#RUN mamba install -v -y untangle
#RUN mamba install -v -y urllib3
#RUN mamba install -v -y xmltodict
#RUN mamba install -v -y xtensor
#RUN mamba install -v -y xtensor-python
#RUN mamba install -v -y wget
#
#RUN mamba install -v -y  "gdal>3.1" "rasterio>1.0.0" "setuptools<58" "shapely<2.0.0"

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
