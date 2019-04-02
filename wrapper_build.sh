#!/bin/bash

MQF_DIR="./ThirdParty/MQF/build"
if [ -d ${MQF_DIR} ] 
then
    echo "Building the wrapper..." 
else
    echo "Building MQF static lib"
    mkdir -p ${MQF_DIR}
    cd ${MQF_DIR}
    cmake ..
    make
    cd ../../../
    echo "Building the wrapper..."
fi


rm -rf -f build/ __pycache__/ *cxx *pyc kProcessor.py *so
swig -c++ -python -py3 kProcessor.i
python3 setup.py build_ext --inplace