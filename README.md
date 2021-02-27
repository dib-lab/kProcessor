<p align="center">
  <img src="https://i.imgur.com/YPtoUI7.png" alt="Logo"/>

</p>
<h1 align="center"> @DIB-LAB/kProcessor </h1>
<p align="center">
<a href="https://travis-ci.org/dib-lab/kProcessor"><img alt="PyPI - Python Version" src="https://travis-ci.org/dib-lab/kProcessor.svg?branch=master"></a>
<a href=""><img alt="Open Issues" src="https://img.shields.io/github/issues-raw/dib-lab/kProcessor" height="20"/></a> <a href="https://kprocessor.readthedocs.io/en/latest/"><img alt="Read the Docs" src="https://img.shields.io/readthedocs/kprocessor"></a> <a href="https://github.com/dib-lab/kProcessor/blob/master/LICENSE"><img alt="GitHub" src="https://img.shields.io/github/license/dib-lab/kProcessor"></a> <a href="https://pypi.org/project/kProcessor/#files"><img alt="PyPI - Wheel" src="https://img.shields.io/pypi/wheel/kprocessor"></a> <a href=""><img alt="GitHub release (latest by date)" src="https://img.shields.io/github/v/release/dib-lab/kProcessor"></a> <a href=""><img alt="Maintained" src="https://img.shields.io/badge/Maintained%3F-yes-green.svg" height="20"/></a> <a href="https://pypi.org/project/kProcessor"><img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/kprocessor"></a>
</p>

<details>
<summary>📖 Table of Contents</summary>
<br />

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#table-of-contents)

## ➤ Table of Contents

- [➤ Table of Contents](#-table-of-contents)
- [➤ Introduction](#-introduction)
- [➤ Quick Installation (pip)](#-quick-installation-pip)
- [➤ Build from source](#-build-from-source)
  - [Clone](#clone)
  - [Install dependencies](#install-dependencies)
  - [Build](#build)
    - [CMake options](#cmake-options)
    - [**Build The kProcessor Library**](#build-the-kprocessor-library)
    - [**Build Everything**](#build-everything)
- [➤ Manually build the Python bindings](#-manually-build-the-python-bindings)
  - [Generate bindings](#generate-bindings)
- [➤ Contributors](#-contributors)
- [➤ License](#-license)

</details>


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#introduction)

## ➤ Introduction

**kProcessor** is a C++ API with a handy Python interface that enables easy handling of sequence kmer content. kProcessor stores kmers with their associated metadata in a virtual data structure called kDataframe. By default, kDataFrame stores the kmers with their counts in the input dataset. However, kDataFrame supports adding multiple other columns to store more information about these kmers in different data types. Users can easily merge multiple kDataFrames or apply different set functions (e.g. union, intersect and difference) on a group of kDataFrames. Colored kDataFrame is another core virtual data structure in kProcessor that allows the indexing of the kmers in a multi-sequence reference input. It is composed of a kDataFrame that replaces the kmer count with a key (aka color). This key connects the kmer to all sequences associated with this kmer.


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#quick_installation)

## ➤ Quick Installation (pip)

```bash
python -m pip install kProcessor
```

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#build_source)

## ➤ Build from source

### Clone

```bash
git clone https://github.com/dib-lab/kProcessor.git
cd kProcessor/
git submodule update --init --recursive
```


### Install dependencies

```bash
sudo apt-get install g++ swig cmake python3-dev zlib1g-dev libghc-bzlib-dev python3-distutils libboost-all-dev
```

### Build

#### CMake options

|      Description      |     Option     | Default |
|:---------------------:|:--------------:|:-------:|
|      Build tests      |   BUILD_TESTS  |   OFF   |
| Build kProcessor apps |   BUILD_APPS   |   OFF   |
|     Build usecases    | BUILD_USECASES |   OFF   |
|  Build documentation  |   BUILD_DOCS   |   OFF   |
|   Build everything!   |    BUILD_ALL   |   OFF   |

#### **Build The kProcessor Library**

```bash=
# Run CMake configure
cmake -Bbuild

# Run make with parallel execution.
cmake --build build -j4 # -j4 = execute 4 recipes simultaneously.
```

#### **Build Everything**

```bash=
cmake -Bbuild -DBUILD_ALL=1

cmake --build build -j4 # -j4 = execute 4 recipes simultaneously.
```

<details><summary><b>Build Docs Only</b></summary>

Output directory: `build/doxygen/html`

```bash
cmake -Bbuild -DBUILD_DOCS=1
cmake --build build --target GenerateDocs
```

</details>

<br>

<details><summary><b>Build Tests Only</b></summary>

```bash
cmake -Bbuild -DBUILD_TESTS=1
cmake --build build -j4

# Run tests
cd build/tests/kProcessorLibTests
./testKprocessorLib
```

</details>

<br>

<details><summary><b>Build Usecases and Apps</b></summary>

```bash
cmake -Bbuild -DBUILD_USECASES=1 -DBUILD_APPS=1
cmake --build build -j4
```

</details>




[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#manual_build_python)

## ➤ Manually build the Python bindings

Python bindings are generated using [SWIG](https://github.com/swig/swig). It's **recommended** to install `swig=4.0.2` using [Conda](https://anaconda.org/conda-forge/swig/).

You can build the python bindings by executing `build_wrapper.sh`, or you can follow the next steps.

### Generate bindings

1. First, you need to follow the instructions in the [Build from source](#build_source).
2. While `pwd=kProcessor` run: `python setup.py bdist_wheel`.
3. Install the generated wheel package using: `cd dist && python -m pip install kProcessor*.whl`.


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#contributors)

## ➤ Contributors
	

| [<img alt="Mostafa Shokrof" src="https://avatars3.githubusercontent.com/u/5207616?s=400&v=4" width="100">](https://github.com/shokrof) | [<img alt="You?" src="https://avatars2.githubusercontent.com/u/7165864?s=460&&v=4" width="100">](https://github.com/mr-eyes) | [<img alt="Tamer Mansour" src="https://avatars3.githubusercontent.com/u/6537740?s=400&&v=4" width="100">](https://github.com/drtamermansour) |
|:--------------------------------------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------:| -------------------------------------------------------------------------------------------------------------------------------------------- |
|                                             [Mostafa Shokrof](https://github.com/shokrof)                                              |                                       [Mohamed Abuelanin](https://github.com/mr-eyes)                                        | [Tamer Manosur](https://github.com/drtamermansour)                                                                          |

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#license)

## ➤ License
	
Licensed under [BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).
