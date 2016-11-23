# from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("NPTFit.NPLL", ["NPTFit/NPLL.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.PLL", ["NPTFit/PLL.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.incgamma_fct", ["NPTFit/incgamma_fct.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.x_m", ["NPTFit/x_m.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3'])
]

setup(
    name='NPTFit',
    version='0.0.1.dev1',
    description='A Python package for Non-Poissoian Template Fitting',
    url='https://github.com/bsafdi/NPTFit',
    author='Benjamin Safdi',
    author_email='bsafdi@princeton.edu',
    license='MIT',
    install_requires=[
            'numpy',
            'matplotlib',
            'scipy',
            'healpy',
            'Cython',
            'pymultinest',
            'jupyter',
            'corner',
            'mpmath',
        ],

    packages=['NPTFit'],
    ext_modules = cythonize(extensions)
)
