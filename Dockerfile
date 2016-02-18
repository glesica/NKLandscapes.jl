FROM ubuntu:14.04

ARG version=releases

RUN apt-get -y update && \
    apt-get -y install software-properties-common && \
    apt-add-repository ppa:staticfloat/julia${version} && \
    apt-add-repository ppa:staticfloat/julia-deps && \
    apt-get -y update && \
    apt-get -y install git julia

VOLUME /opt/src

WORKDIR /opt/src

CMD julia --color=yes --check-bounds=yes -e 'Pkg.clone(pwd(), "NKLandscapes"); Pkg.build("NKLandscapes"); Pkg.test("NKLandscapes")'
