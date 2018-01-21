FROM ubuntu:xenial

# Install dependencies of the os
RUN apt-get update && apt-get install -y software-properties-common && rm -rf /var/lib/apt/lists/*

# Install dependencies of LibsDyogen
RUN add-apt-repository universe && apt-get update &&\
 echo 'Install core dependencies' &&\
 apt-get install -y python2.7 git cython &&\
 echo 'Install marginal dependencies' &&\
 apt-get install -y python-numpy &&\
 echo 'Install marginal dependencies from the universe deposit' &&\
 apt-get install -y python-matplotlib python-scipy &&\
 rm -rf /var/lib/apt/lists/*

# Install LibsDyogen in INSTALL_DIR
ENV INSTALL_DIR "/home/${USER}/Dev"
RUN mkdir -p ${INSTALL_DIR}
WORKDIR ${INSTALL_DIR} 
RUN git clone https://github.com/DyogenIBENS/LibsDyogen &&\
 bash LibsDyogen/cythonisePyxFiles.sh LibsDyogen

# Install plugged-in softwares
# NB:
#   zip is used for unzipping, in homology teams only
#   build-essential provides make, used in homology teams and i-ADHoRe 3.0
RUN apt-get update && apt-get install -y wget zip build-essential && rm -rf /var/lib/apt/lists/*

# Install homology teams
# NB: homolgy teams depends on gcc compiler for make
WORKDIR ${INSTALL_DIR}
RUN wget http://euler.slu.edu/~goldwasser/homologyteams/homologyteams-1.1.zip &&\
 unzip homologyteams-1.1.zip &&\
 cd homologyteams-1.1/src &&\
 make
# To plug homologyteams to LibsDyogen
# sed -i "/PATH_HOMOLOGYTEAMS_BIN =/c\PATH_HOMOLOGYTEAMS_BIN = \"${INSTALL_DIR}/homologyteams-1.1/src/homologyteams\"" ${PATH_LIBSDYOGEN}/utils/myGeneTeams.py

# Install i-ADHoRe 3.0
# Follow the INSTALL file (i-adhore-3.0.01/INSTALL) provided by the ADHoRe team
# NB: needs the universe deposit for mpi: add-apt-repository universe
WORKDIR ${INSTALL_DIR}
RUN apt-get update &&\
 echo 'Install core dependencies' &&\
 apt-get install -y cmake g++ &&\
 echo 'Install marginal dependencies' &&\
 apt-get install -y libpng-dev zlib1g-dev &&\
 echo 'Install marginal dependencies from the universe deposit' &&\
 apt-get install -y mpi &&\
 rm -rf /var/lib/apt/lists/* &&\
 wget http://bioinformatics.psb.ugent.be/downloads/psb/i-adhore/i-adhore-3.0.01.tar.gz &&\
 tar -zxvf i-adhore-3.0.01.tar.gz &&\
 rm i-adhore-3.0.01.tar.gz &&\
 cd i-adhore-3.0.01 &&\
 mkdir build &&\
 cd build &&\
 cmake .. &&\
 make
# You do not need to install it, skip the make install
# To plug i-adhore to LibsDyogen
# sed -i "/PATH_ADHORE_BIN =/c\PATH_ADHORE_BIN = \"${INSTALL_DIR}/i-adhore-3.0.01/build/src/i-adhore\"" ${PATH_LIBSDYOGEN}/utils/myADHoRe.py

# Install Cyntenator
WORKDIR ${INSTALL_DIR}
# download the cyntenator files (pointed by https://www.bioinformatics.org/cyntenator/wiki/Main/HomePage)
RUN apt-get update &&\
 echo 'Install core dependencies' &&\
 apt-get install -y g++ &&\
 rm -rf /var/lib/apt/lists/* &&\
 wget -r -np -nH --cut-dirs=3 -R index.html https://bbc.mdc-berlin.de/svn/bioinformatics/Software/cyntenator/ &&\
 cd cyntenator &&\
 g++ -Wno-deprecated cyntenator.cpp localign.cpp genome.cpp flow.cpp species_tree.cpp -o cyntenator
# To plug cyntenator to LibsDyogen
# sed -i "/PATH_CYNTENATOR_BIN =/c\PATH_CYNTENATOR_BIN = \"${INSTALL_DIR}/cyntenator/cyntenator\"" ${PATH_LIBSDYOGEN}/utils/myCyntenator.py

# Put this at the end since LibsDyogen/enum.py conflicts with python3.5
# If ENV PATH_LIBSDYOGEN "${INSTALL_DIR}/LibsDyogen" is just after the installation step of LibsDyogen:
#    Traceback (most recent call last):
#      File "/usr/bin/add-apt-repository", line 11, in <module>
#        from softwareproperties.SoftwareProperties import SoftwareProperties, shortcut_handler
#      File "/usr/lib/python3/dist-packages/softwareproperties/SoftwareProperties.py", line 49, in <module>
#        from xml.sax.saxutils import escape
#      File "/usr/lib/python3.5/xml/sax/saxutils.py", line 6, in <module>
#        import os, urllib.parse, urllib.request
#      File "/usr/lib/python3.5/urllib/request.py", line 88, in <module>
#        import http.client
#      File "/usr/lib/python3.5/http/__init__.py", line 1, in <module>
#        from enum import IntEnum
#      File "/home/Dev/LibsDyogen/enum.py", line 66
#        raise NotImplementedError, \
#                                 ^
#    SyntaxError: invalid syntax
ENV PATH_LIBSDYOGEN "${INSTALL_DIR}/LibsDyogen"
ENV PYTHONPATH "${PYTHONPATH}:${PATH_LIBSDYOGEN}"

# Install PhylDiag
WORKDIR ${INSTALL_DIR}
# Clone the PhylDiag deposit
RUN git clone https://github.com/DyogenIBENS/PhylDiag &&\
 chmod +x PhylDiag/src/*.py &&\
 chmod +x PhylDiag/src/analysis/*.py &&\
 chmod +x PhylDiag/src/postprocessing/*.py &&\
 cp PhylDiag/src/phylDiag.py /usr/local/bin/phylDiag.py &&\
 cp PhylDiag/src/phylDiagViewer.py /usr/local/bin/phylDiagViewer.py

# Docker containers mix stdout with stderr
# Use this to only output stdout (stream 1)
RUN echo '#!/usr/bin/env bash\n\
date=$(date +%Y-%m-%d-%H-%M-%S)\n\
$@ 2>> ./${date}.stderr\n'\
>> ./only_stdout.sh &&\
 chmod +x ./only_stdout.sh &&\
 mv ./only_stdout.sh /usr/local/bin/only_stdout
# Use this to only output stderr (stream 2), usefull for debugging
RUN echo '#!/usr/bin/env bash\n\
date=$(date +%Y-%m-%d-%H-%M-%S)\n\
$@ 1>> ./${date}.stdout\n'\
>> ./only_stderr.sh &&\
 chmod +x ./only_stderr.sh &&\
 mv ./only_stderr.sh /usr/local/bin/only_stderr

# HOWTO to use this Dockerfile :

# First of all build the docker image
#   docker build -t phyldiagc .

# The easiest way to use this container is to:
# 1) Define  an input folder (INDIR) and an output folder (OUTDIR) for a 2-way sharing of data between host and container
#       INDIR="/home/${USER}/Dev/PhylDiag/data" && OUTDIR="/tmp"
# 2) Start a bash console into the container
#       docker run -v ${INDIR}:/IN -v ${OUTDIR}:/OUT -ti phyldiagc bash
# 3) Execute commands using data in IN/ folder and create results in the OUT/ folder, e.g.:
#       PhylDiag/src/phylDiag.py /IN/Homo.sapiens.genome.bz2 /IN/Mus.musculus.genome.bz2 /IN/Euarchontoglires.families.bz2 > /OUT/res.sbs

# Otherwise you can also use some features of the docker container from outside

# Use phylDiag.py
#   output: wherever you want with the ">" at the end of the command
#   input: ${INDIR}/Homo.sapiens.genome.bz2 ${INDIR}/Mus.musculus.genome.bz2 ${INDIR}/Euarchontoglires.families.bz2
# 1) Set the INDIR (inputdirectory), for example:
#   INDIR="/home/${USER}/Dev/PhylDiag/data"
# 2) Run phylDiag.py
#   docker run -v ${INDIR}:/IN -ti phyldiagc only_stdout PhylDiag/src/phylDiag.py /IN/Homo.sapiens.genome.bz2 /IN/Mus.musculus.genome.bz2 /IN/Euarchontoglires.families.bz2 > ./HM_MM.sbs
# DEBUG) output STDERR for debugging purposes
#   docker run -v ${INDIR}:/IN -ti phyldiagc only_stderr PhylDiag/src/phylDiag.py /IN/Homo.sapiens.genome.bz2 /IN/Mus.musculus.genome.bz2 /IN/Euarchontoglires.families.bz2 > /tmp/err.log

# Use phylDiagImageViewer.py
#   output: ${OUTDIR}/image.svg
#   input: ${INDIR}/Homo.sapiens.genome.bz2 ${INDIR}/Mus.musculus.genome.bz2 ${INDIR}/Euarchontoglires.families.bz2
# 1) Set the INDIR (inputdirectory) and OUTDIR (output directory), for example:
#   INDIR="/home/${USER}/Dev/PhylDiag/data" && OUTDIR="/tmp/"
# 2) Run phyldiagViewer.py
#   docker run -v ${INDIR}:/IN -v ${OUTDIR}:/OUT -ti phyldiagc only_stderr phylDiagViewer.py /IN/Homo.sapiens.genome.bz2 /IN/Mus.musculus.genome.bz2 /IN/Euarchontoglires.families.bz2 /OUT/image.svg > /tmp/err.log
# NB: you can debug phylDiagImageViewer by looking into /tmp/err.log