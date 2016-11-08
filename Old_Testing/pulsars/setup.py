#from distutils.core import setup
#from Cython.Build import cythonize

#setup(
#    ext_modules = cythonize("pulsars/*.pyx")
#)



from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

# extensions = [Extension("pulsars/*.pyx",include_dirs=[numpy.get_include()])]
#setup(
#     ext_modules=[
#                  Extension("pulsars/*.pyx",
#                            include_dirs=[numpy.get_include()]),
#                  ],
#     )

# Or, if you use cythonize() to make the ext_modules list,
# include_dirs can be passed to setup()

setup(
	# ext_modules = cythonize(extensions)
	ext_modules=cythonize("*.pyx"),
	include_dirs=[numpy.get_include()]
     )

#["pulsars/likelihood_psf.c","pulsars/special.c"]