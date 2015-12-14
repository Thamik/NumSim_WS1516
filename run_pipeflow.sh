# clear VTK folder
cd VTK
rm *
cd ..


cd geom_gen
make -j 2
./geomgen_release 2500 2
cd ..
make -j 2
mpirun -np 1 numsim data/assignment03.params data/complex_default.geom
