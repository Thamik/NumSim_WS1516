rm -rf CMakeFiles CMakeCache.txt cmake_install.cmake Makefile
cmake .
make clean
make -j 4

cd geom_gen
make clean
make -j 4
cd ..
