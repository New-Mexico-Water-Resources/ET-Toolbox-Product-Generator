.ONESHELL:
SHELL=bash
VERSION := $(shell cat ETtoolbox/version.txt)

default:
	make install

version:
	$(info ETtoolbox Collection 2 pipeline version ${VERSION})

mamba:
ifeq ($(word 1,$(shell conda run -n base conda list mamba | grep mamba)),mamba)
	@echo "mamba already installed"
else
	-conda deactivate; conda install -y -c conda-forge mamba
endif

mamba-docker:
ifeq ($(word 1,$(shell conda list mamba | grep mamba)),mamba)
	@echo "mamba already installed"
else
	-conda install -y -c conda-forge mamba
endif

create-blank-env:
	$(info creating blank ETtoolbox environment)
	-conda run -n base mamba create -n ETtoolbox

update-env-mamba:
	-conda run -n ETtoolbox mamba env update --file ETtoolbox.yml
#	-conda activate ETtoolbox; mamba env update --file ETtoolbox.yml
#	-source activate ETtoolbox; mamba env update --file ETtoolbox.yml

environment:
	make mamba
	-conda deactivate; pip install pybind11_cmake
	make create-blank-env
	make update-env-mamba

refresh-env:
	make remove
	make environment

environment-docker:
	pip install pybind11_cmake
	mamba env update -n base --file ETtoolbox_docker.yml

clean:
	$(info cleaning build)
	-rm -rvf build
	-rm -rvf dist
	-rm -rvf *.egg-info
	-rm -rvf CMakeFiles
	-rm CMakeCache.txt

uninstall:
	$(info uninstalling ETtoolbox package)
	-conda run -n ETtoolbox pip uninstall ETtoolbox -y
#	-conda activate ETtoolbox; pip uninstall ETtoolbox -y
#	-source activate ETtoolbox; pip uninstall ETtoolbox -y

unit-tests:
	$(info running unit tests)
	conda run -n ETtoolbox nosetests -v -w tests
#	conda activate ETtoolbox; nosetests -v -w tests
#	source activate ETtoolbox; nosetests -v -w tests

unit-tests-docker:
	nosetests -v -w tests

setuptools:
	-conda run -n ETtoolbox python setup.py install
#	-conda activate ETtoolbox; python setup.py install
#	-source activate ETtoolbox; python setup.py install

install-package:
	$(info installing ETtoolbox package)
#	conda run -n ETtoolbox python setup.py install
	-make setuptools
	make clean
	make unit-tests

install-package-docker:
	python setup.py install
	make clean
	make unit-tests-docker

install:
	make environment
	make clean
	make uninstall
	make install-package

install-docker:
	make clean
	make install-package-docker

remove:
	conda run -n base conda env remove -n ETtoolbox

reinstall-hard:
	make remove
	make install

reinstall-soft:
	make uninstall
	make install-package

docker-build:
	docker build -t ETtoolbox .

docker-build-mamba:
	docker build --target mamba -t ETtoolbox .

docker-build-environment:
	docker build --target environment -t ETtoolbox .

docker-shell:
	docker run -it ETtoolbox bash
