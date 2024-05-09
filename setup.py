import glob
from setuptools import setup, Extension
import pysam
from Cython.Build import cythonize

module1 = Extension(name='TRHelper._TRHelper',
                    sources = ['helper.c'],
                    extra_compile_args=["-O3"])

module2 = Extension(name='TRHelper._pysam_ext',
                     sources=['split_bam.pyx'],
                     extra_compile_args=["-O3"],
                     extra_link_args=pysam.get_libraries(),
                     include_dirs=pysam.get_include(),
                     define_macros=pysam.get_defines())
module2 = cythonize(module2)

setup (name = 'TotalReCall',
       version = '0.17',
       description = 'Code for calling somatic transduction events '
           'from case/control paired BAM files',
       packages = ['TRHelper'],
       scripts = glob.glob("scripts/*.py") + glob.glob("scripts/*.sh"),
       ext_modules = [module1] + module2)
