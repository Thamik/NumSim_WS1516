# clear VTK folder
cd VTK
rm *
cd ..

# compile and run karman vortex street
cd geom_gen
make -j 2
./geomgen_release 7500 4
cd ..
sh compile.sh
mpirun -np 4 numsim data/assignment03.params data/complex_default.geom



