cd geom_gen
make -j 2
./geomgen_release 5000 4
cd ..
make -j 2
mpirun -np 4 numsim data/better.params data/complex_default.geom