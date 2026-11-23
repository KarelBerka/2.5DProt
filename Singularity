Bootstrap: docker
From: registry.fedoraproject.org/fedora:32

%help
    This is a container for 2DProt project
    See https://gitlab.com/jhutar/bilkoviny/
    Email ivarekova@centrum.cz

%labels
    Homepage https://gitlab.com/jhutar/bilkoviny/
    Author ivarekova@centrum.cz
    Maintainer jhutar@seznam.cz
    Version 0.1

%post
    dnf -y install \
        bc \
        dnf-utils \
        dotnet-sdk-3.1 \
        dssp \
        ImageMagick \
        pymol \
        python2-numpy \
        python3-biopython \
        python3-numpy \
        python3-pip \
        python3-pyyaml \
        python3-requests \
        python3-scikit-learn \
        python3-svgwrite \
        python3-tabulate \
        python3-wheel \
        unzip
    rm -rf /var/cache/dnf/*
    mkdir /bilkoviny
    chown -R 1000:1000 /bilkoviny
    pip3 install \
        ete3 \
        matplotlib \
        numba \
        numpy \
        requests \
        six \
        svgwrite

%test
    cd /bilkoviny || exit 1
    ###mono cytos.exe --help >/dev/null; [ $? -eq 255 ] || exit 1
    mkdssp --version >/dev/null || exit 1
    ###python -m unittest discover || exit 1

%environment
    export LC_ALL=C

%runscript
    fortune
