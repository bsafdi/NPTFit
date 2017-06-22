from __future__ import print_function

import logging
from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

# Define an easy function for calling the shell
def call_shell(args):
    import subprocess
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return p.returncode, out, err

class BadExitStatus(Exception):
    pass

# Call the pkg-config command and either return a list of flags, or a list with the '-L' flag parts removed.
def call_pkgconfig(args, remove_flags=False):
    code, out, _ = call_shell(['pkg-config'] + args)
    if code != 0:
        raise BadExitStatus()
    flags = out.strip().split(' ')
    if len(flags) == 1 and flags[0] == '':
        return []
    if remove_flags:
        return [ f[2:] for f in flags if len(f) > 2 ]
    else:
        return flags

# Try to find GSL with pkg-config
def find_gsl(libraries=[], library_dirs=[], extra_link_args=[], include_dirs=[], extra_compile_args=[]):
    default = dict(libraries = ['gsl', 'gslcblas', 'm'])
    try:
        result = dict(
            libraries = call_pkgconfig(['--libs-only-l', 'gsl'], remove_flags=True) + libraries,
            library_dirs = call_pkgconfig(['--libs-only-L', 'gsl'], remove_flags=True) + library_dirs,
            extra_link_args = call_pkgconfig(['--libs-only-other', 'gsl']) + extra_link_args,
            include_dirs = call_pkgconfig(['--cflags-only-I', 'gsl'], remove_flags=True) + include_dirs,
            extra_compile_args = call_pkgconfig(['--cflags-only-other', 'gsl']) + extra_compile_args
            )
    except OSError:
        return default
    except BadExitStatus:
        return default

    # remove empty args
    result = { key: value for key, value in result.items() if len(value) > 0 }
    return result


extensions = [
    Extension("NPTFit.npll", ["NPTFit/npll.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.pll", ["NPTFit/pll.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.incgamma_fct_p", ["NPTFit/incgamma_fct_p.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.x_m", ["NPTFit/x_m.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3']),
    Extension("NPTFit.incgamma_fct", ["NPTFit/incgamma_fct.pyx"], **find_gsl(
            include_dirs = [numpy.get_include()],
            extra_compile_args = ["-ffast-math",'-O3']))
]

setup_args = {'name':'NPTFit',
    'version':'0.2',
    'description':'A Python package for Non-Poissonian Template Fitting',
    'url':'https://github.com/bsafdi/NPTFit',
    'author':'Siddharth Mishra-Sharma',
    'author_email':'smsharma@princeton.edu',
    'license':'MIT',
    'install_requires':[
            'numpy',
            'matplotlib',
            'healpy',
            'Cython',
            'pymultinest',
            'jupyter',
            'corner',
            'mpmath',
        ]}


# Attempt GSL compilation; if this fails, do standard compilation.

try:
    print("Attempting GSL compilation...")
    setup(packages=['NPTFit'],
        ext_modules = cythonize(extensions),
        **setup_args
    )
    print("GSL compilation successful!")

except:
    print("GSL compilation failed! Attempting mpmath compilation...")
    setup(packages=['NPTFit'],
        ext_modules = cythonize(extensions[:-1]),
        **setup_args
    )
    print("mpmath compilation successful!")
