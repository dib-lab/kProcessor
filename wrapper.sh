#!/usr/bin/env bash
cd build && make && cd ..

rm -rf build/temp* build/lib.linux* dist/* __pycache__/ *cxx *pyc swig_interfaces/kProcessor_wrap.cpp kProcessor.py *so kProcessor.egg-info/ swig_interfaces/kProcessor_wrap.cpp dist/*

$(which python) setup.py bdist_wheel

cd dist/

$(which python) -m pip uninstall kProcessor -y

echo $(which python)

$(which python) -m pip install --no-cache-dir -U kProcessor*cp*.whl

cd ..

#rm -rf -f __pycache__/ build/temp* build/lib.linux* *cxx *pyc swig_interfaces/kProcessor_wrap.cpp kProcessor.py *so kProcessor.egg-info/ # swig_interfaces/kProcessor_wrap.cpp #dist/*
