cd geom_gen
make -j 2
./geomgen_release 2500 2
cd ..
make -j 2
mpirun -np 4 numsim data/better.params data/complex_default.geom
