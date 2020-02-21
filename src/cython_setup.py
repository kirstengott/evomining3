from distutils.core import Extension, setup
from Cython.Build import cythonize

# define an extension that will be cythonized and compiled
ext = Extension(name="ai_compute", sources=["ai_compute.pyx"])
setup(ext_modules=cythonize(ext))
