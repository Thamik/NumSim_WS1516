This readme should be an instruction on how to run the different geometry settings
properly and how to do the particle tracing.

1) Different geometry settings
==========================================================================================
1.1) Run the complete Script
The easiest, but most time consumpting way to run all four different geometry settings
is to run the script "run_assignment03.sh". This will execute the karman vortex street,
the driven cavity, the pipeflow and the pipeflow with step (in exactly that order) on
4 processes and grid with approximately 7500 cells.
In order to do so, the specific geometry and parameter files are automatically generated
by the script.
The VTK-files of each simulation will be stored in a different folder.

1.2) Run the single script
For each of the three general geometries (hence, not for the driven-cavity) there is
a script "run_[...].sh" which executes the specific geometry setting.

1.3) Run driven-cavity
The driven-cavity can be simulated by simply skipping the geometry file when calling the
program, e.g. execute
mpirun -np 4 numsim data/assignment03.params
for a simulation of the driven-cavity with the parameters listed in "assignment03.params".


2) Particle tracing and streaklines
==========================================================================================
In the course of each simulation, a file "timestep_all_particles.py" is written to the VTK
folder. This contains the position data of the particle tracing which is currently done
by the simulation tool. These data can be visualized by executing the python script
"plot_particles.py". Note, that the file hosting the particle positions has to be stored in
the VTK folder (otherwise the path, implemented in the python script has to be changed to
the new path)!
To change to streaklines can be done in the compute-class of the program.


3) Generation of geometry files
==========================================================================================
3.1 Building of the geometry generator
To build the geometry generator, one can either execute the makefile located in the folder
"geom_gen" or execute one of the scripts mentioned in 1) (all these will compile the
geometry generator).

3.2 Use of the geometry generator
After sucessfully compiling the generator, it can be executed on the terminal. It can handle
up to the arguments: The first one is the number of cells wich should be generated. Note that
this number cannot always be hit exaclty! The second argument referes to the type of geometry:
	1: Driven-cavity
	2: Pipeflow
	3: Pipeflow with step
	4: Karman vortex street
Finally, if the karman vortex street should be generated, one can also pass the angle of the
obstacle as a third argument.

3.3 Generation of parameter files
By default, the geometry generator will write also a parameter file according to the specifications
on the assignment sheet.

