import os
import platform
import shutil
from setuptools.command.build_ext import build_ext

class CustomBuildExt(build_ext):
    def run(self):
        super().run()
        # Determine the current OS and copy the corresponding shared library
        lib_dir = os.path.join(os.path.dirname(__file__), 'my_module_shared_libs')
        target_dir = os.path.join(self.build_lib, 'my_module')
        print(platform.system())
        # if platform.system() == 'Linux':
        #     lib_file = os.path.join(lib_dir, 'linux', 'my_lib.so')
        # elif platform.system() == 'Darwin':
        #     lib_file = os.path.join(lib_dir, 'macos', 'my_lib.dylib')
        # elif platform.system() == 'Windows':
        #     lib_file = os.path.join(lib_dir, 'windows', 'my_lib.dll')
        # else:
        #     raise RuntimeError("Unsupported OS")

        # # Copy the library to the target directory
        # shutil.copy(lib_file, target_dir)
