from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import subprocess

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cfg = 'Debug' if self.debug else 'Release'
        
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=qupled",
        ]
        
        build_args = ['--config', cfg]
        os.makedirs(self.build_temp, exist_ok=True)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

setup(
    # name="qupled",
    # version="0.1.0",
    # description="A package with quantum and classical utilities.",
    # packages=["qupled"],
    ext_modules=[CMakeExtension("qupled.native.native", sourcedir="qupled/native/src")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)
