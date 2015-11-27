#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "communicator.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* atof */

#include <mpi.h>

/* public methods */

/* constructor */
/**
Constructs a Geometry object with standart values<br>
- 128x128 grid
- physical x-length: 1
- physical y-length: 1 
*/
Geometry::Geometry()
: Geometry(NULL) // TODO: is there a nicer way?
{
	std::cout << "Warning: Geometry Constructor: no communicator given!\n" << std::flush;
}

Geometry::Geometry(Communicator *comm)
: _comm(comm), _size(128,128), _bsize(128,128), _length(1.0,1.0), _blength(1.0,1.0), _h(1.0,1.0), _velocity(0.0,0.0), _pressure(1.0) //standard values
//: _comm(comm), _size(32,32), _bsize(32,32), _length(1.0,1.0), _blength(1.0,1.0), _h(1.0,1.0), _velocity(0.0,0.0), _pressure(1.0)
{
	// handle total/partial size/length values
	set_meshwidth(); // set _h to the right values
	// TODO: something else to do?
}

/**
The geometry for the actual problem is given in the file in path file.
\param[in] file file path
\param[in] verbose printing debugging information
*/
void Geometry::Load(const char *file, bool verbose)
{
for (int i_rank=0; i_rank < _comm->getSize(); i_rank++){
if(i_rank == _comm->getRank()){ // read the files sequentially

	if (verbose && _comm->getRank()==0){
		std::cout << "Loading geometry file from path " << file << " ...\n";
	}
	std::string temp_string;
	std::ifstream infile;
	infile.open(file);
	if (!infile.is_open()){
		std::cout << "Warning: geometry file could not be read!\n";
		continue;
	}
	for(int i=1;i<=7;i++)
	{
		getline(infile,temp_string);
		switch(i)
		{
			case 1:
				_bsize[0] = atoi(temp_string.c_str());
				break;
			case 2:
				_bsize[1] = atoi(temp_string.c_str());
				break;
			case 3:
				_blength[0] = atof(temp_string.c_str());
				break;
			case 4:
				_blength[1] = atof(temp_string.c_str());
				break;
			case 5:
				_velocity[0] = atof(temp_string.c_str());
				break;
			case 6:
				_velocity[1] = atof(temp_string.c_str());
				break;
			case 7:
				_pressure = atof(temp_string.c_str());
				break;
		}
	}
	infile.close();
	//set_meshwidth(); // update _h
	//if (verbose){
	if (_comm->getRank() == 0){
		std::cout << "--------------------------------------------------\n";
		std::cout << "Geometry configuration read from file:\n";
		std::cout << "Total Size\t=\t(" << _bsize[0] << ", " << _bsize[1] << ")\n";
		std::cout << "Total Length\t=\t(" << _blength[0] << ", " << _blength[1] << ")\n";
		std::cout << "Velocity\t=\t(" << _velocity[0] << ", " << _velocity[1] << ")\n";
		std::cout << "Pressure\t=\t" << _pressure << "\n";
		std::cout << "--------------------------------------------------\n";
	}

}
MPI_Barrier(MPI_COMM_WORLD);
}

}

/**
\return number of cells in each direction
*/
const multi_index_t& Geometry::Size() const
{
	return _size;
}

const multi_index_t& Geometry::TotalSize() const
{
	return _bsize;
}

/**
\return physical length of the domain
*/
const multi_real_t& Geometry::Length() const
{
	return _length;
}

const multi_real_t &Geometry::TotalLength() const
{
	return _blength;
}

/**
\return width of a single grid cell
*/
const multi_real_t& Geometry::Mesh() const
{
	return _h;
}

/**
\return pointer to the communicator instance
*/
const Communicator* Geometry::getCommunicator() const
{
	return _comm;
}


//==============================================================================
// Methods for the update of boundary values (have to be done every timestep)
//==============================================================================

// TODO: something to edit here?

/**
Updates the boundary values for Grid u according to the pattern of u, i.e. Dirichlet boundary conditions
\param[out] u Grid in which the boundary conditions should be set
*/
void Geometry::Update_U(Grid *u) const
{
	//std::cout << "Geometry: " << _size[0] << ", " << _size[1] << ", " << _bsize[0] << ", " << _bsize[1] << ", " << _h[0] << ", " << _h[1] << "\n" << std::flush;

	//std::cout << "Geometry: Update_U\n" << std::flush; // only for debugging issues
	// see lecture, 3.1.2
	for (int i=1; i<=4; i++){
		// check if this is a domain boundary
		if (!is_global_boundary(i)) continue;

		//std::cout << "Process " << _comm->getRank() << " updates boundary " << i << ".\n" << std::flush;

		BoundaryIterator it(this);
		it.SetBoundary(i);
		it.First();
		while (it.Valid()){
			//std::cout << "i = " << i << ", it = " << it.Value() << "\n" << std::flush; // only for debugging issues
			//if (i==4){
			if (i==BoundaryIterator::boundaryTop){
				// upper boundary
				u->Cell(it) = 2*1.0 - u->Cell(it.Down());
				//std::cout << _comm->getRank() << ": setting top boundary values\n";
			//} else if (i==1){
			} else if (i==BoundaryIterator::boundaryLeft) {
				// left boundary
				u->Cell(it) = 0.0;
				//std::cout << _comm->getRank() << ": setting left boundary values\n";
			//} else if (i==2) {      
			} else if (i==BoundaryIterator::boundaryRight) {
				// right boundary
				u->Cell(it) = 0.0;
				u->Cell(it.Left()) = 0.0;
				//std::cout << _comm->getRank() << ": setting right boundary values\n";
			} else {
				// lower boundary
				u->Cell(it) = - u->Cell(it.Top());
				//std::cout << _comm->getRank() << ": setting bottom boundary values\n";
			}
			it.Next();
		}
	}
}

/**
Updates the boundary values for Grid v according to the pattern of v, i.e. Dirichlet boundary conditions
\param[out] v Grid in which the boundary conditions should be set
*/
void Geometry::Update_V(Grid *v) const
{
	// see lecture, 3.1.2
	for (int i=1; i<=4; i++){
		// check if this is a domain boundary
		if (!is_global_boundary(i)) continue;

		BoundaryIterator it(this);
		it.SetBoundary(i);
		it.First();
		while (it.Valid()){
			//if (i==4){
			if (i==BoundaryIterator::boundaryTop){
				// upper boundary
				v->Cell(it) = 0.0;
				v->Cell(it.Down()) = 0.0;
			//} else if (i==3){
			} else if (i==BoundaryIterator::boundaryBottom) {
				// lower boundary
				v->Cell(it) = 0.0;
			//} else if (i == 1) {
			} else if (i==BoundaryIterator::boundaryLeft) {
				v->Cell(it) = -v->Cell(it.Right());
			//} else if (i == 2) {
			} else if (i==BoundaryIterator::boundaryRight) {
				v->Cell(it) = -v->Cell(it.Left());
			}
			it.Next();
		}
	}
}

/**
Updates the boundary values for Grid p according to the pattern of p, i.e. Neumann boundary conditions
\param[out] p Grid in which the boundary conditions should be set
*/
void Geometry::Update_P(Grid *p) const
{
	// see lecture, 3.2.3
	for (int i=1; i<=4; i++){
		// check if this is a domain boundary
		if (!is_global_boundary(i)) continue;

		BoundaryIterator it(this);
		it.SetBoundary(i);
		it.First();
		while (it.Valid()){
			//if (i==4){
			if (i==BoundaryIterator::boundaryTop){
				// upper boundary
				p->Cell(it) = p->Cell(it.Down());
			//} else if (i==1){
			} else if (i==BoundaryIterator::boundaryLeft) {
				// left boundary
				p->Cell(it) = p->Cell(it.Right());
			//} else if (i==2){
			} else if (i==BoundaryIterator::boundaryRight) {
				// right boundary
				p->Cell(it) = p->Cell(it.Left());
			//} else if (i==3){
			} else if (i==BoundaryIterator::boundaryBottom) {
				// lower boundary
				p->Cell(it) = p->Cell(it.Top());
			}
			it.Next();
		}
	}
}

// own method
/**
Calculates the meshwidth such that it satisfies h = length/size
*/
void Geometry::set_meshwidth()
{
	// update meshwidth
	_h[0] = _blength[0]/(_bsize[0]-2);
	_h[1] = _blength[1]/(_bsize[1]-2);

	// update local size and length
	_size[0] = _comm->getLocalSize()[0];
	_size[1] = _comm->getLocalSize()[1];
	_length[0] = (_size[0]-2) * _h[0];
	_length[1] = (_size[1]-2) * _h[1];

	//std::cout << "set_meshwidth: " << _comm->getLocalSize()[0] << ", " << _comm->getLocalSize()[1] << "\n" << std::flush;
	//std::cout << "set_meshwidth: " << _size[0] << ", " << _size[1] << ", " << _bsize[0] << ", " << _bsize[1] << ", " << _h[0] << ", " << _h[1] << "\n" << std::flush;
}

/**
Decides whether the given boundary is a global boundary or is not
*/
bool Geometry::is_global_boundary(int boundary_index) const
{
	switch (boundary_index)
	{
		case BoundaryIterator::boundaryBottom :
			return _comm->isBottom();
			break;
		case BoundaryIterator::boundaryLeft :
			return _comm->isLeft();
			break;
		case BoundaryIterator::boundaryTop :
			return _comm->isTop();
			break;
		case BoundaryIterator::boundaryRight :
			return _comm->isRight();
			break;
		default:
			std::cout << "Warning: Geometry is_global_boundary: invalid boundary_index given!\n" << std::flush;
			return false;
			break;
	}
}

/**
This function updates the local geometry data stored in this object, using the latest information provided by the related communicator.
*/
void Geometry::update_values()
{
	//std::cout << "Geometry: update values!\n" << std::flush;
	set_meshwidth();
	//std::cout << "Geometry: " << _size[0] << ", " << _size[1] << ", " << _bsize[0] << ", " << _bsize[1] << ", " << _h[0] << ", " << _h[1] << "\n" << std::flush;
}

/**
Does the domain decomposition on the master process and broadcasts the resulting information, concerning the process distribution, to all other processes.
*/
void Geometry::do_domain_decomposition()
{
	multi_index_t tdim;
	int** rankDistri;
	multi_index_t** localSizes;
	if (_comm->getRank() == 0) {
		//std::cout << "Domain decomposition...\n" << std::flush;

		// horizontal decomposition
		//horizontal_domain_decomposition(tdim, rankDistri, localSizes);
		// vertical decomposition
		//vertical_domain_decomposition(tdim, rankDistri, localSizes);
		// 2d decomposition
		rect_domain_decomposition(tdim, rankDistri, localSizes);

		// send information
		//std::cout << "Broadcasting domain decomposition information...\n" << std::flush;
		int sendBuf(tdim[0]);
		MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		sendBuf = tdim[1];
		MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//std::cout << "First broadcasting completed.\n" << std::flush;

		for (int i=0; i<tdim[0]; i++){
			MPI_Bcast(rankDistri[i], tdim[1], MPI_INT, 0, MPI_COMM_WORLD);
		}

		//std::cout << "Second broadcasting completed.\n" << std::flush;

		for (int i=0; i<tdim[0]; i++){
			for (int j=0; j<tdim[1]; j++) {
				sendBuf = localSizes[i][j][0];
				MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				sendBuf = localSizes[i][j][1];
				MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//std::cout << "All broadcasting completed.\n" << std::flush;
	} else {
		// receive information
		int recBuf(0);
		MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		tdim[0] = recBuf;
		MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		tdim[1] = recBuf;

		//std::cout << "Recieved tdim: " << tdim[0] << ", " << tdim[1] << "\n" << std::flush;

		// allocate
		rankDistri = new int*[tdim[0]];
		for (int i=0; i<tdim[0]; i++){
			rankDistri[i] = new int[tdim[1]];
		}
		localSizes = new multi_index_t*[tdim[0]];
		for (int i=0; i<tdim[0]; i++){
			localSizes[i] = new multi_index_t[tdim[1]];
		}

		for (int i=0; i<tdim[0]; i++){
			MPI_Bcast(rankDistri[i], tdim[1], MPI_INT, 0, MPI_COMM_WORLD);
		}
		for (int i=0; i<tdim[0]; i++){
			for (int j=0; j<tdim[1]; j++) {
				MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				localSizes[i][j][0] = recBuf;
				MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				localSizes[i][j][1] = recBuf;
				//std::cout << "Recieved localSizes[i][j]: " << localSizes[i][j][0] << ", " << localSizes[i][j][1] << "\n" << std::flush;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	//std::cout << "Process no. " << _comm->getRank() << " is still working (1)!!\n" << std::flush;

	_bsize[0] += 2;
	_bsize[1] += 2;

	// write to communicator
	_comm->setProcDistribution(rankDistri, tdim, localSizes);

	//std::cout << "Process no. " << _comm->getRank() << " is still working (2)!!\n" << std::flush;

	// update geometry values
	update_values();

	//std::cout << "Process " << _comm->getRank() << ": left(" << is_global_boundary(BoundaryIterator::boundaryLeft) << "), right(" << is_global_boundary(BoundaryIterator::boundaryRight) << "), top(" << is_global_boundary(BoundaryIterator::boundaryTop) << "), bottom(" << is_global_boundary(BoundaryIterator::boundaryBottom) << ")\n" << std::flush;

	if(_comm->getRank()==0) {
		std::cout << "--------------------------------------------------\n";
		std::cout << "Process distribution (" << _comm->getSize() << " processes in total):\n" << std::flush;
		for(int j=tdim[1]-1; j >= 0; j--) {
			for(int i=0; i < tdim[0]; i++) {
				if (i != 0) std::cout << "\t-\t";
				std::cout << rankDistri[i][j];
			}
			std::cout << "\n" << std::flush;
		}
		std::cout << "--------------------------------------------------\n";
	}
}

/**
Does a horizontal decomposition of the physical domain.
*/
void Geometry::horizontal_domain_decomposition(multi_index_t& tdim, int**& rankDistri, multi_index_t**& localSizes) const
{
	index_t np = _comm->getSize();
	tdim[0] = np;
	tdim[1] = 1;
	// allocate
	rankDistri = new int*[tdim[0]];
	for (int i=0; i<tdim[0]; i++){
		rankDistri[i] = new int[tdim[1]];
	}
	localSizes = new multi_index_t*[tdim[0]];
	for (int i=0; i<tdim[0]; i++){
		localSizes[i] = new multi_index_t[tdim[1]];
	}
	// write values
	for (int i=0; i<tdim[0]; i++){
		for (int j=0; j<tdim[1]; j++){
			rankDistri[i][j] = i + j*tdim[0];
			if (i == (tdim[0]-1)){
				localSizes[i][j][0] = TotalSize()[0] - int(TotalSize()[0] / tdim[0]) * (tdim[0]-1);
			} else {
				localSizes[i][j][0] = TotalSize()[0] / tdim[0];
			}
			if (j == (tdim[1]-1)){
				localSizes[i][j][1] = TotalSize()[1] - int(TotalSize()[1] / tdim[1]) * (tdim[1]-1);
			} else {
				localSizes[i][j][1] = TotalSize()[1] / tdim[1];
			}
		}
	}

	//add ghost cells to domains
	/*for (int i=0; i<tdim[0]; i++) {
		for (int j=0; j<tdim[1]; j++) {
			// add ghost cells in x-dim
			if (i == 0 && i == (tdim[0]-1)) {
				// do nothing
			} else if (i == 0 || i == (tdim[0]-1)) {
				localSizes[i][j][0] += 1;
			} else {
				localSizes[i][j][0] += 2;
			}

			// add ghost cells in y-dim
			if (j == 0 && j == (tdim[1]-1)) {
				// do nothing
			} else if (j == 0 || j == (tdim[1]-1)) {
				localSizes[i][j][1] += 1;
			} else {
				localSizes[i][j][1] += 2; 
			}
		}
	}*/

	for (int i=0; i < tdim[0]; i++) {
		for (int j = 0; j < tdim[1]; j++) {
			localSizes[i][j][0] += 2;
			localSizes[i][j][1] += 2;
		}
	}
}

/**
Does a vertical decomposition of the physical domain.
*/
void Geometry::vertical_domain_decomposition(multi_index_t& tdim, int**& rankDistri, multi_index_t**& localSizes) const
{
	index_t np = _comm->getSize();
	tdim[0] = 1;
	tdim[1] = np;
	// allocate
	rankDistri = new int*[tdim[0]];
	for (int i=0; i<tdim[0]; i++){
		rankDistri[i] = new int[tdim[1]];
	}
	localSizes = new multi_index_t*[tdim[0]];
	for (int i=0; i<tdim[0]; i++){
		localSizes[i] = new multi_index_t[tdim[1]];
	}
	// write values
	for (int i=0; i<tdim[0]; i++){
		for (int j=0; j<tdim[1]; j++){
			rankDistri[i][j] = i + j*tdim[0];
			if (i == (tdim[0]-1)){
				localSizes[i][j][0] = TotalSize()[0] - int(TotalSize()[0] / tdim[0]) * (tdim[0]-1);
			} else {
				localSizes[i][j][0] = TotalSize()[0] / tdim[0];
			}
			if (j == (tdim[1]-1)){
				localSizes[i][j][1] = TotalSize()[1] - int(TotalSize()[1] / tdim[1]) * (tdim[1]-1);
			} else {
				localSizes[i][j][1] = TotalSize()[1] / tdim[1];
			}
		}
	}

	//add ghost cells to domains
	/*for (int i=0; i<tdim[0]; i++) {
		for (int j=0; j<tdim[1]; j++) {
			// add ghost cells in x-dim
			if (i == 0 && i == (tdim[0]-1)) {
				// do nothing
			} else if (i == 0 || i == (tdim[0]-1)) {
				localSizes[i][j][0] += 1;
			} else {
				localSizes[i][j][0] += 2;
			}

			// add ghost cells in y-dim
			if (j == 0 && j == (tdim[1]-1)) {
				// do nothing
			} else if (j == 0 || j == (tdim[1]-1)) {
				localSizes[i][j][1] += 1;
			} else {
				localSizes[i][j][1] += 2; 
			}
		}
	}*/

	for (int i=0; i < tdim[0]; i++) {
		for (int j = 0; j < tdim[1]; j++) {
			localSizes[i][j][0] += 2;
			localSizes[i][j][1] += 2;
		}
	}
}

/**
Tries to decompose the domain in a two-dimensional fashion. To do so, the function performs a simple prime factorization of the total number of processes.
*/
void Geometry::rect_domain_decomposition(multi_index_t& tdim, int**& rankDistri, multi_index_t**& localSizes) const
{
	index_t np = _comm->getSize();
	// search for a good decomposition
	tdim[0] = np;
	tdim[1] = 1;
	while (tdim[0] > tdim[1]){
		bool divisor_found = false;
		for (int i=2; i<=(np/2); i++){
			if ((tdim[0] % i) == 0){
				index_t temp_x = tdim[0] / i;
				index_t temp_y = tdim[1] * i;
				if (temp_x >= temp_y){
					tdim[0] = temp_x;
					tdim[1] = temp_y;
					divisor_found = true;
				}
				break;
			}
		}
		if (!divisor_found) break;
	}
	// allocate
	rankDistri = new int*[tdim[0]];
	for (int i=0; i<tdim[0]; i++){
		rankDistri[i] = new int[tdim[1]];
	}
	localSizes = new multi_index_t*[tdim[0]];
	for (int i=0; i<tdim[0]; i++){
		localSizes[i] = new multi_index_t[tdim[1]];
	}
	// write the other values
	for (int i=0; i<tdim[0]; i++){
		for (int j=0; j<tdim[1]; j++){
			rankDistri[i][j] = i + j*tdim[0];
			if (i == (tdim[0]-1)){
				localSizes[i][j][0] = TotalSize()[0] - int(TotalSize()[0] / tdim[0]) * (tdim[0]-1);
			} else {
				localSizes[i][j][0] = TotalSize()[0] / tdim[0];
			}
			if (j == (tdim[1]-1)){
				localSizes[i][j][1] = TotalSize()[1] - int(TotalSize()[1] / tdim[1]) * (tdim[1]-1);
			} else {
				localSizes[i][j][1] = TotalSize()[1] / tdim[1];
			}
		}
	}

	//add ghost cells to domains
	/*for (int i=0; i<tdim[0]; i++) {
		for (int j=0; j<tdim[1]; j++) {
			// add ghost cells in x-dim
			if (i == 0 && i == (tdim[0]-1)) {
				// do nothing
			} else if (i == 0 || i == (tdim[0]-1)) {
				localSizes[i][j][0] += 1;
			} else {
				localSizes[i][j][0] += 2;
			}

			// add ghost cells in y-dim
			if (j == 0 && j == (tdim[1]-1)) {
				// do nothing
			} else if (j == 0 || j == (tdim[1]-1)) {
				localSizes[i][j][1] += 1;
			} else {
				localSizes[i][j][1] += 2; 
			}
		}
	}*/

	for (int i=0; i < tdim[0]; i++) {
		for (int j = 0; j < tdim[1]; j++) {
			localSizes[i][j][0] += 2;
			localSizes[i][j][1] += 2;
		}
	}
}
