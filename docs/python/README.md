# kProcessor
<a href="https://github.com/drtamermansour/nu-ngs02/blob/master/kProcessor/Tutorial.ipynb" rel="Jupyter Notebook Example">![Foo](https://img.shields.io/badge/Jupyter-notebook-green)</a>
## Installation  

### Requirements

`sudo apt-get install g++ swig cmake python3-dev zlib1g-dev libghc-bzlib-dev python3-distutils`

## Install from PyPi

`pip install kProcessor`

## Install Manually

### Clone

```bash
git clone --recursive --branch devel https://github.com/dib-lab/kProcessor.git kProcessor
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

---

## Documentation

- [kDataFrame](./kDataFrame.markdown)
- [colored_kDataFrame](./coloredKDataFrame.markdown)
- [kDataFrame Iterator](./kDataFrameIterator.markdown)
- [kmerDecoder](./kmerDecoder.markdown)
- [kmer counting](./kmerCounting.markdown)
- [Indexing](./indexing.markdown)
- [set functions](./setFunctions.markdown)
