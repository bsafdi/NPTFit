from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("NPTFit.npll", ["NPTFit/npll.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.pll", ["NPTFit/pll.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.incgamma_fct", ["NPTFit/incgamma_fct.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.x_m", ["NPTFit/x_m.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3'])
]

setup(
    name='NPTFit',
    version='0.1.1',
    description='A Python package for Non-Poissonian Template Fitting',
    url='https://github.com/bsafdi/NPTFit',
    author='Siddharth Mishra-Sharma',
    author_email='smsharma@princeton.edu',
    license='MIT',
    install_requires=[
            'numpy',
            'matplotlib',
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
