from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import subprocess


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cfg = "Debug" if self.debug else "Release"
        use_mpi = os.environ.get("USE_MPI")
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            f"-DUSE_MPI={use_mpi}",
        ]
        # if use_mpi is not None:
        #     cmake_args = cmake_args + use_mpi
        build_args = ["--config", cfg]
        os.makedirs(self.build_temp, exist_ok=True)
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


setup(
    ext_modules=[CMakeExtension("qupled.native", sourcedir="src/qupled/native/src")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    package_dir={"": "src"},
)
