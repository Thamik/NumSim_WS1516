cd geom_gen
make -j 2
./geomgen_release 1000 5 0
cd ..
make -j 2
mpirun -np 1 numsim data/better.params data/complex_default.geom
