cd geom_gen
make -j 2
./geomgen_release 250 3
cd ..
make -j 2
mpirun -np 1 numsim data/better.params data/complex_default.geom
