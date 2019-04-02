from distutils.core import setup, Extension
import os
from glob import glob

os.environ["CC"] = "g++"

src_files = glob(os.path.join('src', '*', '*cpp'))
src_files += glob("src/kDataFrame.cpp")
src_files += ["kProcessor_wrap.cxx"]

includes = ["include/kProcessor", "ThirdParty/CLI",
            "ThirdParty/MQF/include", "ThirdParty/seqan/include"]

kProcessor_module = Extension(name='_kProcessor',
                              library_dirs=['ThirdParty/MQF/build/',
                                            "/usr/lib/x86_64-linux-gnu/"],
                              libraries=["lMQF", "z"],
                              sources=src_files,
                              extra_compile_args=['-Wall', "-lgomp", '-Wextra', '-std=c++14', '-fPIC', '-fopenmp',
                                                  '-W', '-Wall', '-pedantic', '-lrt', '-lpthread', '-lbz2', '-lz', '-DSEQAN_HAS_ZLIB=1',
                                                  '-DSEQAN_HAS_BZIP2=1', '-DSEQAN_HAS_OPENMP=1', "-DSEQAN_HAS_EXECINFO=1"],

                              extra_link_args=["-fopenmp", "-lgomp", "-lbz2"],
                              include_dirs=includes,
                              language='c++'

                              )

setup(name='kProcessor',
      version='0.1',
      author="M. Abuelanin",
      description="""kProcessor pythonic interface""",
      ext_modules=[kProcessor_module],
      py_modules=["kProcessor"],
      )
