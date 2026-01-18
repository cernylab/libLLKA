import os
import sys
import subprocess
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # Required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        # Find Eigen
        eigen_dir = os.environ.get('EIGEN_INCLUDE_DIR')
        if not eigen_dir:
            # Try common locations
            for candidate in ['/usr/include/eigen3', '/usr/local/include/eigen3',
                            '/opt/homebrew/include/eigen3', '/opt/local/include/eigen3']:
                if os.path.exists(os.path.join(candidate, 'Eigen', 'Dense')):
                    eigen_dir = candidate
                    break

        if not eigen_dir:
            raise RuntimeError(
                "Eigen library not found. Please set EIGEN_INCLUDE_DIR environment variable "
                "or install Eigen (e.g., 'sudo zypper install eigen3-devel')"
            )

        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',
            '-DBUILD_PYTHON_BINDINGS=ON',
            f'-DEIGEN_INCLUDE_DIR={eigen_dir}',
        ]

        build_args = []

        # Set CMAKE_BUILD_TYPE
        cfg = 'Debug' if self.debug else 'Release'
        cmake_args += [f'-DCMAKE_BUILD_TYPE={cfg}']

        # Parallel build
        if hasattr(self, 'parallel') and self.parallel:
            build_args += [f'-j{self.parallel}']
        else:
            build_args += ['-j4']

        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        # CMake configuration
        subprocess.check_call(
            ['cmake', ext.sourcedir] + cmake_args,
            cwd=self.build_temp
        )

        # Build
        subprocess.check_call(
            ['cmake', '--build', '.', '--target', '_pyllka_core'] + build_args,
            cwd=self.build_temp
        )

        # Debug: Print paths
        import shutil
        import glob
        print(f"\n=== DEBUG INFO ===")
        print(f"build_temp: {self.build_temp}")
        print(f"extdir: {extdir}")
        print(f"Looking for module in: {self.build_temp}/python/pyllka/_pyllka_core*.so")

        built_module = glob.glob(f'{self.build_temp}/python/pyllka/_pyllka_core*.so')
        print(f"Found modules: {built_module}")

        if built_module:
            print(f"Copying {built_module[0]} -> {extdir}")
            # Ensure the destination directory exists
            os.makedirs(extdir, exist_ok=True)
            shutil.copy(built_module[0], extdir)
        else:
            print("ERROR: Built module not found!")
            print(f"Contents of {self.build_temp}/python:")
            if os.path.exists(f"{self.build_temp}/python"):
                for root, dirs, files in os.walk(f"{self.build_temp}/python"):
                    print(f"  {root}:")
                    for f in files:
                        print(f"    {f}")
        print(f"=== END DEBUG ===\n")

setup(
    ext_modules=[CMakeExtension('pyllka._pyllka_core', sourcedir='..')],
    cmdclass={'build_ext': CMakeBuild},
)
