# Getting Started

## Installation  

### Install from Pypi

> You will need to have the [pip](https://pypi.org/project/pip/) python package manager installed.

`pip install kProcessor`

### Install Manually

#### Requirements

`sudo apt-get install g++ swig cmake python3-dev zlib1g-dev libghc-bzlib-dev python3-distutils`

#### Clone

```bash
git clone --recursive --branch devel https://github.com/dib-lab/kProcessor.git kProcessor
cd kProcessor
```

#### Build kProcessorLib

```bash
mkdir build
cd build/
cmake ..
make
cd ..
```

#### Install kProcessor Python Package

```bash
python3 setup.py install
```

## Quick Start
