for numproc in 9 16
do
	echo $numproc
	mpirun -np $numproc numsim data/my_default.params data/medium_grid.geom
done
