from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize, build_ext
import os
import numpy

NVCC='nvcc'

def customize_compiler(self):
    """Modified from: https://github.com/rmcgibbo/npcuda-example/tree/master/cython"""
    
    global NVCC
    
    default_compiler_so = self.compiler_so
    super = self._compile
    
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        self.set_executable('compiler_so', NVCC)  # TODO:
        cc_args_o = cc_args
        cc_args += extra_postargs[NVCC]
        super(obj, src, ext, cc_args, [], pp_opts)
        self.compiler_so = default_compiler_so
    self._compile = _compile

class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler(self.compiler)
        build_ext.build_extensions(self)


examples_extension = Extension(
    name="pyreconstruction",
    sources=["pyreconstruction.pyx"],
    libraries=["beamforming","armadillo","gomp","cudart_static","rt"],
    library_dirs=["../../lib","/usr/local/cuda-10.2/lib64"],
    include_dirs=["../../include", numpy.get_include()],
    language='c++',
    extra_compile_args={ NVCC: ["-O3", "-DARMA_ALLOW_FAKE_GCC", "-Xcompiler", "-fopenmp", "-Xcompiler", "-fPIC"]}
)
setup(
    name="pyreconstruction",
    ext_modules=cythonize([examples_extension], compiler_directives={'language_level' : "3"}), 
    cmdclass={'build_ext': custom_build_ext}
)
