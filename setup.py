#!/usr/bin/env python

from distutils.command.build import build
from distutils.spawn import find_executable
import sys
import os
import subprocess
import errno
from glob import glob

KPROCESSOR = r"""
  _    _____                                        
 | |  |  __ \                                       
 | | _| |__) | __ ___   ___ ___  ___ ___  ___  _ __ 
 | |/ /  ___/ '__/ _ \ / __/ _ \/ __/ __|/ _ \| '__|
 |   <| |   | | | (_) | (_|  __/\__ \__ \ (_) | |   
 |_|\_\_|   |_|  \___/ \___\___||___/___/\___/|_|                                                                                                        
"""

if sys.version_info[:2] < (3, 6):
    raise RuntimeError("Python version >=3.6")

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

if os.path.islink("KP_BUILD"):
    os.unlink("KP_BUILD")

if os.path.exists("build/libkProcessor.a"):
    os.symlink("build", "KP_BUILD")


def get_version():
    return "2.0.2"

def check_exist(dirs):
    ALL_EXIST = True
    not_found_files = list()
    for directory in dirs:
        if not (os.path.isdir(directory)):
            print(f"[ERROR] | DIR: {directory} does not exist.", file=sys.stderr)
            ALL_EXIST = False
            not_found_files.append(directory)

    if not ALL_EXIST:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), ",".join(not_found_files))


SOURCES = [
    'swig_interfaces/kProcessor.i',
]

if not find_executable('swig'):
    sys.exit("Error:  Building this module requires 'swig' to be installed")

INCLUDES = [
    'ThirdParty/ntCard/include',
    'include/kProcessor',
    'ThirdParty/MQF/include',
    'ThirdParty/sdsl-lite/include',
    'ThirdParty/kmerDecoder/include',
    'ThirdParty/kmerDecoder/lib/parallel-hashmap',
    'ThirdParty/kmerDecoder/lib/kseq/include',
    'ThirdParty/Blight',
    'ThirdParty/mum-hash',
    'ThirdParty/KMC/kmc_api',
    'ThirdParty/caches/include',
]

check_exist(INCLUDES)

LINK_ARGS = [
    "-fopenmp",
    "-lgomp",
    "-lbz2",
    "-lz",
    "-ldl",
]

kp_build_dir = "KP_BUILD"

LIBRARIES_DIRS = [
    f"{kp_build_dir}",
    # "ThirdParty/KMC/kmc_api",
    "ThirdParty/Blight",
    f"{kp_build_dir}/ThirdParty/MQF/src",
    "ThirdParty/ntCard",
    f"{kp_build_dir}/ThirdParty/sdsl-lite/lib",
    f"{kp_build_dir}/ThirdParty/kmerDecoder",
    f"{kp_build_dir}/ThirdParty/MQF/ThirdParty/stxxl/lib",

]

check_exist(LIBRARIES_DIRS)

RUNTIME_LIBRARIES_DIRS = [
    'ThirdParty/KMC/kmc_api',
]

LIBRARIES = [
    # 'KMCAPI',
    'kProcessor',
    'blight',
    'sdsl',
    'MQF',
    'ntcard',
    'kmerDecoder',
    'stxxl_debug',
]

SWIG_OPTS = [
    '-c++',
    '-py3',
    '-outdir',
    '.',
    '-Isrc',
    '-doxygen',
]


class CustomBuild(build):
    sub_commands = [
        ('build_ext', build.has_ext_modules),
        ('build_py', build.has_pure_modules),
        ('build_clib', build.has_c_libraries),
        ('build_scripts', build.has_scripts),
    ]


# BLIGHT_HEADERS = glob("ThirdParty/Blight/*h")

# print(BLIGHT_HEADERS)

kProcessor_module = Extension('_kProcessor',
                              # runtime_library_dirs=RUNTIME_LIBRARIES_DIRS,
                              library_dirs=LIBRARIES_DIRS,
                              libraries=LIBRARIES,
                              sources=SOURCES,
                              include_dirs=INCLUDES,
                            #   includes=BLIGHT_HEADERS,
                              extra_link_args=LINK_ARGS,
                              extra_compile_args=["-O3", "-Ofast", "-std=c++17", "-fPIC"],
                              swig_opts=SWIG_OPTS,
                              )

classifiers = [
    "License :: OSI Approved :: Apache Software License",
    'Development Status :: 3 - Alpha',
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
]

commit_hash_short_name = subprocess.getoutput("git rev-parse --short HEAD").split()[0]
branch_name = subprocess.getoutput("git rev-parse --abbrev-ref HEAD").split()[0]

setup(name='kProcessor',
      version=get_version(),
      author="Tamer Mansour, Mostafa Shokrof, Mohamed Abuelanin",
      author_email='drtamermansour@gmail.com, mostafa.shokrof@gmail.com, mabuelanin@gmail.com',
      description="""kProcessor Python interface""",
      ext_modules=[kProcessor_module],
      py_modules=['kProcessor'],
      url='https://github.com/dib-lab/kProcessor',
      python_requires='>=3.6',
      cmdclass={'build': CustomBuild},
      license='BSD 3-Clause',
      long_description_content_type='text/markdown',
      long_description=readme,
      classifiers=classifiers,
      include_package_data=True,
      project_urls={
          'Bug Reports': 'https://github.com/dib-lab/kProcessor/issues',
          'Source': 'https://github.com/dib-lab/kProcessor',
      },
      )

if os.path.exists("build/libkProcessor.a") and os.path.islink("KP_BUILD"):
    os.unlink("KP_BUILD")
