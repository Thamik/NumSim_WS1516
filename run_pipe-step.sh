cd geom_gen
make -j 2
./geomgen_release 7500 3
cd ..
make -j 2
mpirun -np 4 numsim data/assignment03.params data/complex_default.geom
