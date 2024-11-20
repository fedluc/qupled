import os
from cmake_build_extension import CMakeExtension, BuildExtension
from setuptools import setup

setup(
    ext_modules=[
        CMakeExtension(
            name="qupled.native",
            source_dir="src/qupled/native/src",
            cmake_configure_options=[
                f"-DUSE_MPI={os.environ.get('USE_MPI', 'ON')}",
            ],
        )
    ],
    cmdclass={"build_ext": BuildExtension},
    zip_safe=False,
    package_dir={"": "src"},
)
