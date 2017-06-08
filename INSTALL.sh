#!/bin/bash

# check that libsDyogen is in the PYTHONPATH
exec 3<<'HERE'
try:
    from utilsa import *
    print 0
except ImportError:
    print 1
    raise ImportError('LibsDyogen is not in the PYTHONPATH')
HERE
res=$(python /dev/fd/3)
echo ${res}

if [ ${res} = "0" ]
then
        echo "LibsDyogen is properly set" >&2
else
	# res = 1
	echo "please install LibsDyogen first" >&2
	echo "either run 'bash https://github.com/DyogenIBENS/LibsDyogen/blob/master/INSTALL.sh' :" >&2
	echo "Or read https://github.com/DyogenIBENS/LibsDyogen/blob/master/README.txt" >&2
	exit
	# TODO, automatical installation of LibsDyogen
fi

PATH_PARENT_PHYLDIAG="/home/${USER}/Libs"
mkdir -p ${PATH_PARENT_PHYLDIAG}
cd ${PATH_PARENT_PHYLDIAG}
PATH_PHYLDIAG=${PATH_PARENT_PHYLDIAG}/PhylDiag
echo "Installing PhylDiag into ${PATH_PHYLDIAG}" >&2
# Clone the PhylDiag deposit
git clone https://github.com/DyogenIBENS/PhylDiag ${PATH_PHYLDIAG}
# If necessary give execution rights
chmod +x ${PATH_PHYLDIAG}/src/*.py
chmod +x ${PATH_PHYLDIAG}/src/analysis/*.py
chmod +x ${PATH_PHYLDIAG}/src/postprocessing/*.py
