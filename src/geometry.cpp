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
: _comm(comm), _size(128,128), _bsize(128,128), _length(1.0,1.0), _blength(1.0,1.0), _h(1.0,1.0), _flags(NULL), _bval_u(NULL), _bval_v(NULL), _bval_p(NULL) //standard values
{
	// handle total/partial size/length values
	set_meshwidth(); // set _h to the right values
	// TODO: something else to do?
}

/**
The destructor of the Geometry class.
*/
Geometry::~Geometry()
{
	if (_flags != nullptr) delete[] _flags;
	if (_bval_u != nullptr) delete[] _bval_u;
	if (_bval_v != nullptr) delete[] _bval_v;
	if (_bval_p != nullptr) delete[] _bval_p;
}

void Geometry::load_domain_partitioning(const char* file)
{
	/*
		expects file = nullptr, if file not assigned as command line argument
	*/

	// variables where the information is stored
	index_t bsizeX, bsizeY;
	real_t blengthX, blengthY;
	char* flags(nullptr);
	real_t* bvu(nullptr);
	real_t* bvv(nullptr);
	real_t* bvp(nullptr);

	if (_comm->getRank() == 0){
		// on the master, load file and do partitioning
		char* total_flags;
		real_t* total_bvu;
		real_t* total_bvv;
		real_t* total_bvp;

		bool read = (file != nullptr);

		std::ifstream infile;
		if (read){
			infile.open(file);
			std::cout << "Loading geometry file from path " << file << " ...\n";
			if (!infile.is_open()){
				std::cout << "Warning: geometry file could not be read!\n";
				read = false;
				infile.close();
			}
		}

		if (!read){
			// set default values
			bsizeX = _bsize[0];
			bsizeY = _bsize[1];
			blengthX = _blength[0];
			blengthY = _blength[1];
		} else {
			// read from file
			infile >> bsizeX;
			infile >> bsizeY;
			infile >> blengthX;
			infile >> blengthY;
		}

		// allocate fields
		total_flags = new char[(bsizeX+2) * (bsizeY+2)];
		total_bvu = new real_t[(bsizeX+2) * (bsizeY+2)];
		total_bvv = new real_t[(bsizeX+2) * (bsizeY+2)];
		total_bvp = new real_t[(bsizeX+2) * (bsizeY+2)];

		// TODO: better!
		infile.get(); // try to read one more character
		infile.peek(); // eof here?
		bool b_eof = infile.eof();
		infile.unget();

		if (!read || b_eof){
			// set default values: driven cavity
			index_t ival;
			for (index_t i=0; i<bsizeX+2; i++){
				for (index_t j=0; j<bsizeY+2; j++){
					ival = j * (bsizeX+2) + i;
					if (j==bsizeY+1) {
						// upper boundary
						total_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v
						total_bvu[ival] = 1.0;
						total_bvv[ival] = 0.0;
						total_bvp[ival] = 0.0;
					} else if (j==0 || i==0 || i==bsizeX+1){
						// lower, left or right boundary
						total_flags[ival] = 1 | 1<<3; // neumann condition for p, dirichlet for u and v
						total_bvu[ival] = 0.0;
						total_bvv[ival] = 0.0;
						total_bvp[ival] = 0.0;
					}
				}
			}
		} else {
			// read from file
			for (index_t i=0; i<(bsizeX+2)*(bsizeY+2); i++){
				int tempint(0);
				infile >> tempint;
				total_flags[i] = char(tempint);
			}
			for (index_t i=0; i<(bsizeX+2)*(bsizeY+2); i++){
				infile >> total_bvu[i];
			}
			for (index_t i=0; i<(bsizeX+2)*(bsizeY+2); i++){
				infile >> total_bvv[i];
			}
			for (index_t i=0; i<(bsizeX+2)*(bsizeY+2); i++){
				infile >> total_bvp[i];
			}
		}

		// close file if necessary
		if (read){
			infile.close();
		}

		//std::cout << "Geometry file reading completed!" << std::endl;

		// send and store simple geometry information
		MPI_Bcast(&bsizeX, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&bsizeY, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&blengthX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&blengthY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		_bsize[0] = bsizeX;
		_bsize[1] = bsizeY;
		_blength[0] = blengthX;
		_blength[1] = blengthY;

		// console output
		std::cout << "--------------------------------------------------\n";
		std::cout << "Geometry configuration:\n";
		std::cout << "Total Size\t=\t(" << _bsize[0] << ", " << _bsize[1] << ")\n";
		std::cout << "Total Length\t=\t(" << _blength[0] << ", " << _blength[1] << ")\n";
		std::cout << "--------------------------------------------------\n";

		// domain partitioning
		multi_index_t tdim;
		int** rankDistri;
		multi_index_t** localSizes;
		do_domain_decomposition(tdim, rankDistri, localSizes);

		// compute indices of corner points
		multi_index_t** cornerPoints;
		cornerPoints = new multi_index_t*[tdim[0]];
		for (index_t i=0; i<tdim[0]; i++){
			cornerPoints[i] = new multi_index_t[tdim[1]];
		}

		for (index_t i=0; i<tdim[0]; i++){
			for (index_t j=0; j<tdim[1]; j++){
				if (i==0 && j==0){
					cornerPoints[i][j] = multi_index_t(0,0);
				} else if (j==0){
					cornerPoints[i][j][0] = cornerPoints[i-1][j][0] + localSizes[i-1][j][0] - 2;
					cornerPoints[i][j][1] = cornerPoints[i-1][j][1];
				} else {
					cornerPoints[i][j][0] = cornerPoints[i][j-1][0];
					cornerPoints[i][j][1] = cornerPoints[i][j-1][1] + localSizes[i][j-1][1] - 2;
				}
			}
		}

		// send information
		char* sendBuf_flags;
		real_t* sendBuf_u;
		real_t* sendBuf_v;
		real_t* sendBuf_p;
		for (index_t i=0; i<tdim[0]; i++){
			for (index_t j=0; j<tdim[1]; j++){

				index_t buflen = localSizes[i][j][0] * localSizes[i][j][1];

				if (rankDistri[i][j] == 0){
					flags = new char[buflen];
					bvu = new real_t[buflen];
					bvv = new real_t[buflen];
					bvp = new real_t[buflen];
					// on master, write to own local storage
					for (index_t n=0; n<localSizes[i][j][0]; n++){
						for (index_t m=0; m<localSizes[i][j][1]; m++){
							index_t totind = (m+cornerPoints[i][j][1])*bsizeX + (n+cornerPoints[i][j][0]);
							flags[m*localSizes[i][j][0]+n] = total_flags[totind];
							bvu[m*localSizes[i][j][0]+n] = total_bvu[totind];
							bvv[m*localSizes[i][j][0]+n] = total_bvv[totind];
							bvp[m*localSizes[i][j][0]+n] = total_bvp[totind];
						}
					}
				} else {
					sendBuf_flags = new char[buflen];
					sendBuf_u = new real_t[buflen];
					sendBuf_v = new real_t[buflen];
					sendBuf_p = new real_t[buflen];
					for (index_t n=0; n<localSizes[i][j][0]; n++){
						for (index_t m=0; m<localSizes[i][j][1]; m++){
							index_t totind = (m+cornerPoints[i][j][1])*bsizeX + (n+cornerPoints[i][j][0]);
							sendBuf_flags[m*localSizes[i][j][0]+n] = total_flags[totind];
							sendBuf_u[m*localSizes[i][j][0]+n] = total_bvu[totind];
							sendBuf_v[m*localSizes[i][j][0]+n] = total_bvv[totind];
							sendBuf_p[m*localSizes[i][j][0]+n] = total_bvp[totind];
						}
					}
					
					MPI_Send(sendBuf_flags, buflen, MPI_CHAR, rankDistri[i][j], 0, MPI_COMM_WORLD);
					MPI_Send(sendBuf_u, buflen, MPI_DOUBLE, rankDistri[i][j], 0, MPI_COMM_WORLD);
					MPI_Send(sendBuf_v, buflen, MPI_DOUBLE, rankDistri[i][j], 0, MPI_COMM_WORLD);
					MPI_Send(sendBuf_p, buflen, MPI_DOUBLE, rankDistri[i][j], 0, MPI_COMM_WORLD);

					delete[] sendBuf_flags;
					delete[] sendBuf_u;
					delete[] sendBuf_v;
					delete[] sendBuf_p;
				}

			}
		}

		// delete domain data
		for (index_t i=0; i<tdim[0]; i++){
			delete[] rankDistri[i];
		}
		delete[] rankDistri;
		for (index_t i=0; i<tdim[0]; i++){
			delete[] localSizes[i];
		}
		delete[] localSizes;
		for (index_t i=0; i<tdim[0]; i++){
			delete[] cornerPoints[i];
		}
		delete[] cornerPoints;

	} else {
		// not on the master, receive information

		// simple geometry information
		MPI_Bcast(&bsizeX, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&bsizeY, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&blengthX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&blengthY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		_bsize[0] = bsizeX;
		_bsize[1] = bsizeY;
		_blength[0] = blengthX;
		_blength[1] = blengthY;

		// domain partitioning (receive data)
		// does the same as "do_domain_decomposition(multi_index_t(0,0), nullptr, nullptr);"
		int** dummya;
		multi_index_t** dummyb;
		multi_index_t dummyc;
		do_domain_decomposition(dummyc, dummya, dummyb);

		// allocate
		index_t buflen = _size[0] * _size[1];
		flags = new char[buflen];
		bvu = new real_t[buflen];
		bvv = new real_t[buflen];
		bvp = new real_t[buflen];

		// receive complex geometry data
		MPI_Recv(flags, buflen, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(bvu, buflen, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(bvv, buflen, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(bvp, buflen, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}

	// dump the information in this Geometry object (both on master and other processes)
	if (_flags != nullptr) delete[] _flags;
	if (_bval_u != nullptr) delete[] _bval_u;
	if (_bval_v != nullptr) delete[] _bval_v;
	if (_bval_p != nullptr) delete[] _bval_p;

	if (flags == nullptr || bvu == nullptr || bvv == nullptr || bvp == nullptr){
		// safety check, if all fields are initialized
		std::cout << "Warning: An error did occur while loading the geometry data. Uninitialized field!\n" << std::flush;
	} else {
		_flags = flags;
		_bval_u = bvu;
		_bval_v = bvv;
		_bval_p = bvp;
	}

	// dont deleting local storage, because the member pointer not point to this data

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
	for (index_t i=1; i<=4; i++){
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
	for (index_t i=1; i<=4; i++){
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
	for (index_t i=1; i<=4; i++){
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


/*==================================================================
New Update Methods for General Geometry (GG)
==================================================================*/
void Geometry::UpdateGG_U(Grid *u) const
{
	BoundaryIteratorGG it(this);
	it.First();

	while(it.Valid()) {
 		// check where the fluid is
		bool topdown = !isObstacle(it.Top()) || !isObstacle(it.Down());
		bool leftright = !isObstacle(it.Left()) || !isObstacle(it.Right());
		if(isNeumannBoundaryU(it)) {
			//TODO: Implement Neumann boundary conditions
		} else {
			if(!isObstacle(it.Left()) && !isObstacle(it.Right())) {
				std::cout << "Warning: The obstacle is too thin (in x-direction)!!!\n" << std::flush;
			} else if (!isObstacle(it.Top()) && !isObstacle(it.Down()) {
				std::cout << "Warning: The obstacle is too thin (in y-direction)!!!\n" << std::flush;
			} else if (topdown) {
				//either at the top or the bottom cell is fluid
				if (!isObstacle(it.Left())) {
					u->Cell(it) = bvalU(it);
				} else if (!isObstacle(it.Right()) {
					u->Cell(it) = bvalU(it);
					u->Cell(it.Right()) = bvalU(it);
				} else {
					if (!isObstacle(it.Top())) {
						u->Cell(it) = 2.0*bvalU(it) - u->Cell(it.Top());
					} else if (!isObstacle(it.Down())) {
						u->Cell(it) = 2.0*bvalU(it) - u->Cell(it.Down());
					}
				}
			} else if (leftright) {
				// the corner cases don't have to be considered as they were already consider above!
				if (!isObstacle(it.Left())) {
					u->Cell(it) = bvalU(it);
				} else if (!isObstacle(it.Right())) {
					u->Cell(it) = bvalU(it);
					u->Cell(it.Right()) = bvalU(it);
				}
			}
		}
	}
}

void Geometry::UpdateGG_V(Grid *v) const
{

}

void Geometry::UpdateGG_P(Grid *p) const
{

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
	set_meshwidth();
}

/**
Does the domain decomposition on the master process and broadcasts the resulting information, concerning the process distribution, to all other processes.
*/
void Geometry::do_domain_decomposition(multi_index_t& tdim, int**& rankDistri, multi_index_t**& localSizes)
{
	if (_comm->getRank() == 0) {
		// horizontal decomposition
		//horizontal_domain_decomposition(tdim, rankDistri, localSizes);
		// vertical decomposition
		//vertical_domain_decomposition(tdim, rankDistri, localSizes);
		// 2d decomposition
		rect_domain_decomposition(tdim, rankDistri, localSizes);

		// send information
		int sendBuf(tdim[0]);
		MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		sendBuf = tdim[1];
		MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

		for (index_t i=0; i<tdim[0]; i++){
			MPI_Bcast(rankDistri[i], tdim[1], MPI_INT, 0, MPI_COMM_WORLD);
		}

		for (index_t i=0; i<tdim[0]; i++){
			for (index_t j=0; j<tdim[1]; j++) {
				sendBuf = localSizes[i][j][0];
				MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				sendBuf = localSizes[i][j][1];
				MPI_Bcast(&sendBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	} else {
		// receive information
		int recBuf(0);
		MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		tdim[0] = recBuf;
		MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		tdim[1] = recBuf;

		// allocate
		rankDistri = new int*[tdim[0]];
		for (index_t i=0; i<tdim[0]; i++){
			rankDistri[i] = new int[tdim[1]];
		}
		localSizes = new multi_index_t*[tdim[0]];
		for (index_t i=0; i<tdim[0]; i++){
			localSizes[i] = new multi_index_t[tdim[1]];
		}

		for (index_t i=0; i<tdim[0]; i++){
			MPI_Bcast(rankDistri[i], tdim[1], MPI_INT, 0, MPI_COMM_WORLD);
		}
		for (index_t i=0; i<tdim[0]; i++){
			for (index_t j=0; j<tdim[1]; j++) {
				MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				localSizes[i][j][0] = recBuf;
				MPI_Bcast(&recBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
				localSizes[i][j][1] = recBuf;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	_bsize[0] += 2;
	_bsize[1] += 2;

	// write to communicator
	_comm->setProcDistribution(rankDistri, tdim, localSizes);

	// update geometry values
	update_values();

	if(_comm->getRank()==0) {
		std::cout << "--------------------------------------------------\n";
		std::cout << "Process distribution (" << _comm->getSize() << " processes in total):\n" << std::flush;
		for(int j=int(tdim[1])-1; j >= 0; j--) {
			for(index_t i=0; i < tdim[0]; i++) {
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
	for (index_t i=0; i<tdim[0]; i++){
		rankDistri[i] = new int[tdim[1]];
	}
	localSizes = new multi_index_t*[tdim[0]];
	for (index_t i=0; i<tdim[0]; i++){
		localSizes[i] = new multi_index_t[tdim[1]];
	}
	// write values
	for (index_t i=0; i<tdim[0]; i++){
		for (index_t j=0; j<tdim[1]; j++){
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
	for (index_t i=0; i<tdim[0]; i++) {
		for (index_t j=0; j<tdim[1]; j++) {
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
	for (index_t i=0; i<tdim[0]; i++){
		rankDistri[i] = new int[tdim[1]];
	}
	localSizes = new multi_index_t*[tdim[0]];
	for (index_t i=0; i<tdim[0]; i++){
		localSizes[i] = new multi_index_t[tdim[1]];
	}
	// write values
	for (index_t i=0; i<tdim[0]; i++){
		for (index_t j=0; j<tdim[1]; j++){
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
	for (index_t i=0; i<tdim[0]; i++) {
		for (index_t j=0; j<tdim[1]; j++) {
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
	int np = _comm->getSize();
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
	for (index_t i=0; i<tdim[0]; i++){
		rankDistri[i] = new int[tdim[1]];
	}
	localSizes = new multi_index_t*[tdim[0]];
	for (index_t i=0; i<tdim[0]; i++){
		localSizes[i] = new multi_index_t[tdim[1]];
	}
	// write the other values
	for (index_t i=0; i<tdim[0]; i++){
		for (index_t j=0; j<tdim[1]; j++){
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
	for (index_t i=0; i<tdim[0]; i++) {
		for (index_t j=0; j<tdim[1]; j++) {
			localSizes[i][j][0] += 2;
			localSizes[i][j][1] += 2;
		}
	}
}

// Getter functions for the complex geometry data
bool Geometry::isObstacle(const Iterator& it) const
{
	return (_flags[it.Value()] >> 0) & 1;
}

bool Geometry::isNeumannBoundaryU(const Iterator& it) const
{
	return (_flags[it.Value()] >> 1) & 1;
}

bool Geometry::isNeumannBoundaryV(const Iterator& it) const
{
	return (_flags[it.Value()] >> 2) & 1;
}

bool Geometry::isNeumannBoundaryP(const Iterator& it) const
{
	return (_flags[it.Value()] >> 3) & 1;
}

// Getter functions for the boundary data
const real_t& Geometry::bvalU(const Iterator& it) const
{
	return _bval_u[it.Value()];
}

const real_t& Geometry::bvalV(const Iterator& it) const
{
	return _bval_v[it.Value()];
}

const real_t& Geometry::bvalP(const Iterator& it) const
{
	return _bval_p[it.Value()];
}
