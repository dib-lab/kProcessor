# kProcessor

## Installation  

### Requirements

`sudo apt-get install g++ swig cmake python3-dev zlib1g-dev libghc-bzlib-dev python3-distutils` |

### Clone

```bash
git clone --recursive --branch development https://github.com/dib-lab/kProcessor.git kProcessor
cd kProcessor
```

### Build kProcessorLib

```bash
mkdir build
cd build/
cmake ..
make
cd ..
```

### Install kProcessor Python Package

```bash
python3 setup.py install
```

### Optional (Build wheel package)

```bash
python3 setup.py bdist_wheel
```
