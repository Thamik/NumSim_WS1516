cd VTK
rm *
cd ..

mkdir VTK_DrivenCav
mkdir VTK_PipeFlow
mkdir VTK_PipeStep
mkdir VTK_Karmann

# delete 
cd VTK_DrivenCav
rm -r VTK
cd ..

cd VTK_PipeFlow
rm -r VTK
cd ..

cd VTK_PipeStep
rm -r VTK
cd ..

cd VTK_Karmann
rm -r VTK
cd ..

# compile main program and geometry generator
sh compile.sh
cd geom_gen
make clean
make -j 2
cd ..

# simulate Karmann Vortex Street
echo 'Simulate Stepflow\n\n'
echo 'Create Geometry...\n'
cd geom_gen
./geomgen_release 2500 4
cd ..
echo 'Geometry created!\n'

echo 'Simulate...\n'
mpirun -np 4 numsim data/assignment03.params data/complex_default.geom
echo 'done!\n'

echo 'Copy VTK data...\n'
cp -af VTK VTK_Karmann
cd VTK
rm *
cd ..
echo 'done!\n'

# simulate Driven Cavity
echo 'Simulate Driven Cavity\n\n'
echo 'Create Geometry...\n'
cd geom_gen
./geomgen_release 2500 1
cd ..
echo 'Geometry created!\n'

echo 'Simulate...\n'
mpirun -np 4 numsim data/assignment03.params data/complex_default.geom
echo 'done!\n'

echo 'Copy VTK data...\n'
cp -af VTK VTK_DrivenCav
cd VTK
rm *
cd ..
echo 'done!\n'

# simulate Pipe Flow
echo 'Simulate Pipeflow\n\n'
echo 'Create Geometry...\n'
cd geom_gen
./geomgen_release 2500 2
cd ..
echo 'Geometry created!\n'

echo 'Simulate...\n'
mpirun -np 4 numsim data/assignment03.params data/complex_default.geom
echo 'done!\n'

echo 'Copy VTK data...\n'
cp -af VTK VTK_PipeFlow
cd VTK
rm *
cd ..
echo 'done!\n'

# simulate Pipe Step FLow
echo 'Simulate Stepflow\n\n'
echo 'Create Geometry...\n'
cd geom_gen
./geomgen_release 2500 3
cd ..
echo 'Geometry created!\n'

echo 'Simulate...\n'
mpirun -np 4 numsim data/assignment03.params data/complex_default.geom
echo 'done!\n'

echo 'Copy VTK data...\n'
cp -af VTK VTK_PipeStep
cd VTK
rm *
cd ..
echo 'done!\n'
