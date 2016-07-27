from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("*", ["*.pyx"],
        include_dirs=[numpy.get_include()], extra_compile_args=["-ffast-math",'-O3'])
]
setup(
    #name = "My hello app",
    ext_modules = cythonize(extensions),
)
