#!/usr/bin/env python

KPROCESSOR = r"""
  _    _____                                        
 | |  |  __ \                                       
 | | _| |__) | __ ___   ___ ___  ___ ___  ___  _ __ 
 | |/ /  ___/ '__/ _ \ / __/ _ \/ __/ __|/ _ \| '__|
 |   <| |   | | | (_) | (_|  __/\__ \__ \ (_) | |   
 |_|\_\_|   |_|  \___/ \___\___||___/___/\___/|_|                                                                                                        
"""

from distutils.command.build import build
from distutils.spawn import find_executable
import sys
import os

print(sys.argv)

if sys.version_info[:2] < (3, 5) or sys.version_info[:2] > (3, 6):
    raise RuntimeError("Python version == 3.6 required.")

if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

try:
    with open('README.md') as f:
        readme = f.read()
except IOError:
    readme = ''

SOURCES = [
    'src/algorithms.cpp',
    'src/kDataFrame.cpp',
    'src/Utils/kmer.cpp',
    'src/HashUtils/hashutil.cpp',
    'src/KmerDecoder/FastqReader.cpp',
    'ThirdParty/MQF/src/gqf.cpp',
    'ThirdParty/MQF/src/utils.cpp',
    'swig_interfaces/kProcessor.i',
]

if not find_executable('swig'):
    sys.exit("Error:  Building this module requires 'swig' to be installed")


INCLUDES = [
    'include/kProcessor',
    'ThirdParty/CLI',
    'ThirdParty/MQF/include',
    'ThirdParty/seqan/include',
    'ThirdParty/TBB/include',
]

COMPILE_ARGS = [
    '-Wall',
    "-lgomp",
    '-Wextra',
    '-std=c++14',
    '-fPIC',
    '-fopenmp',
    '-W',
    '-Wall',
    '-pedantic',
    '-lrt',
    '-lpthread',
    '-lbz2',
    '-lz',
    '-DSEQAN_HAS_ZLIB=1',
    '-DSEQAN_HAS_BZIP2=1',
    '-DSEQAN_HAS_OPENMP=1',
    "-DSEQAN_HAS_EXECINFO=1"
]

LINK_ARGS = [
    "-fopenmp",
    "-lgomp",
    "-lbz2",
]

LIBRARIES_DIRS = [
    '/usr/lib/x86_64-linux-gnu/'
]

LIBRARIES = [
    'z',
]

SWIG_OPTS = [
    '-DSWIGWORDSIZE64',
    '-c++',
    '-py3',
    '-outdir',
    '.',
    '-Isrc'
]

class CustomBuild(build):
    sub_commands = [
        ('build_ext', build.has_ext_modules),
        ('build_py', build.has_pure_modules),
        ('build_clib', build.has_c_libraries),
        ('build_scripts', build.has_scripts),
    ]


kProcessor_module = Extension('_kProcessor',
                          sources=SOURCES,
                          include_dirs=INCLUDES,
      #                    library_dirs=LIBRARIES_DIRS,
                          libraries=LIBRARIES,
                          extra_compile_args=COMPILE_ARGS,
                          extra_link_args=LINK_ARGS,
                          swig_opts=SWIG_OPTS,                          
                          )

classifiers = [
    "License :: OSI Approved :: Apache Software License",
    'Development Status :: 3 - Alpha',
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
]

setup(name='kProcessor',
      version='0.1',
      author="M. Abuelanin",
      author_email='mabuelanin@gmail.com',
      description="""kProcessor Python interface""",
      ext_modules=[kProcessor_module],
      py_modules=['kProcessor'],
      url='https://github.com/dib-lab/kProcessor',
      python_requires='>=3.6, <3.7',
      cmdclass={'build': CustomBuild},
      license='BSD 3-Clause',
      long_description_content_type='text/markdown',
      long_description=readme,
      classifiers=classifiers,
      include_package_data=True,
      project_urls={
        'Bug Reports': 'https://github.com/dib-lab/kProcessor/issues',
        'Say Thanks!': 'https://saythanks.io/to/mr-eyes',
        'Source': 'https://github.com/dib-lab/kProcessor',
        },
      )
