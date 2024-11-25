#! /bin/sh

function print_usage() {
	echo "Usage:"
	echo "mkcc.sh path_to_eigen"
}

function clean_cmake_stuff() {
	rm Makefile
	rm -fr CMakeFiles
	rm cmake_install.cmake
	rm CTestTestfile.cmake
}

CC=gcc
CCX=g++;

IN_EIGEN_INCLUDE_DIR="$1"
if [ -z $IN_EIGEN_INCLUDE_DIR ]; then
	print_usage
	exit 1
fi

cmake . \
	-DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CCX \
	-DCMAKE_BUILD_TYPE=Debug \
	-DEIGEN_INCLUDE_DIR="$IN_EIGEN_INCLUDE_DIR" -DBUILD_EXAMPLES=ON \
	-DDNATCO_EXTRAS=ON \
	-DCMAKE_EXPORT_COMPILE_COMMANDS=ON

rm CMakeCache.txt
rm DartConfiguration.tcl
rm -fr Testing
clean_cmake_stuff

# Clean up tests
cd tests
clean_cmake_stuff

# Clean up examples (fun!)
cd ../examples
clean_cmake_stuff

cd simple_NtC_assignment
clean_cmake_stuff

cd ../gui_assigner
clean_cmake_stuff
