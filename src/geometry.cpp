#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "communicator.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* atof */

#include <math.h>

#include <mpi.h>

// this is only for exercise sheet 3, the two-cell criterion is not satisfied
//#define twocell_criterion_check

/* public methods */

/* constructor */
/**
Constructs a Geometry object with standart values<br>
- 128x128 grid
- physical x-length: 1
- physical y-length: 1 
*/
Geometry::Geometry()
: Geometry(nullptr) // TODO: is there a nicer way?
{
	std::cout << "Warning: Geometry Constructor: no communicator given!\n" << std::flush;
}

Geometry::Geometry(Communicator *comm)
: _comm(comm), _size(128,128), _bsize(128,128), _total_offset(0,0), _offsets(nullptr), _localSizes(nullptr), _length(1.0,1.0), _blength(1.0,1.0), _h(1.0,1.0), _flags(nullptr), _bval_u(nullptr), _bval_v(nullptr), _bval_p(nullptr) //standard values
{
	// handle total/partial size/length values
	set_meshwidth(); // set _h to the right values
	// TODO: something else to do?
}

// copy constructor
Geometry::Geometry(const Geometry* geom)
: _comm(geom->_comm), _size(geom->_size), _bsize(geom->_bsize), _total_offset(geom->_total_offset), 
//_offsets(geom->_offsets), _localSizes(geom->_localSizes), 
_length(geom->_length), _blength(geom->_blength), _h(geom->_h)
//, _flags(geom->_flags), _bval_u(geom->_bval_u), _bval_v(geom->_bval_v), _bval_p(geom->_bval_p)
{
	set_meshwidth(); // set _h to the right values
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

	if (_offsets != nullptr){
		for (index_t i=0; i<_comm->ThreadDim()[0]; i++){
			delete[] _offsets[i];
		}
		delete[] _offsets;
	}

	if (_localSizes != nullptr){
		for (index_t i=0; i<_comm->ThreadDim()[0]; i++){
			delete[] _localSizes[i];
		}
		delete[] _localSizes;
	}
}

void Geometry::load_domain_partitioning(const char* file)
{
	/*
		expects file = nullptr, if file not assigned as the command line argument
	*/

	// variables where the information is stored
	index_t bsizeX, bsizeY;
	real_t blengthX, blengthY;
	char* flags(nullptr);
	real_t* bvu(nullptr);
	real_t* bvv(nullptr);
	real_t* bvp(nullptr);
	index_t totalOffsetX(0);
	index_t totalOffsetY(0);

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
#ifdef OUTPUT_GEOMETRY
			std::cout << "Loading geometry file from path " << file << " ...\n";
#endif
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
					} else {
						total_flags[ival] = 0; // fluid
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
#ifdef OUTPUT_GEOMETRY
		std::cout << "--------------------------------------------------\n";
		std::cout << "Geometry configuration:\n";
		std::cout << "Total Size\t=\t(" << _bsize[0] << ", " << _bsize[1] << ")\n";
		std::cout << "Total Length\t=\t(" << _blength[0] << ", " << _blength[1] << ")\n";
		std::cout << "--------------------------------------------------\n";
#endif

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
					totalOffsetX = cornerPoints[i][j][0];
					totalOffsetY = cornerPoints[i][j][1];

					flags = new char[buflen];
					bvu = new real_t[buflen];
					bvv = new real_t[buflen];
					bvp = new real_t[buflen];
					// on master, write to own local storage
					for (index_t n=0; n<localSizes[i][j][0]; n++){
						for (index_t m=0; m<localSizes[i][j][1]; m++){
							index_t totind = (m+cornerPoints[i][j][1])*_bsize[0] + (n+cornerPoints[i][j][0]);
							flags[m*localSizes[i][j][0]+n] = total_flags[totind];
							bvu[m*localSizes[i][j][0]+n] = total_bvu[totind];
							bvv[m*localSizes[i][j][0]+n] = total_bvv[totind];
							bvp[m*localSizes[i][j][0]+n] = total_bvp[totind];
						}
					}
				} else {
					MPI_Send(&(cornerPoints[i][j][0]), 1, MPI_INT, rankDistri[i][j], 0, MPI_COMM_WORLD);
					MPI_Send(&(cornerPoints[i][j][1]), 1, MPI_INT, rankDistri[i][j], 0, MPI_COMM_WORLD);

					sendBuf_flags = new char[buflen];
					sendBuf_u = new real_t[buflen];
					sendBuf_v = new real_t[buflen];
					sendBuf_p = new real_t[buflen];
					for (index_t n=0; n<localSizes[i][j][0]; n++){
						for (index_t m=0; m<localSizes[i][j][1]; m++){
							index_t totind = (m+cornerPoints[i][j][1])*_bsize[0] + (n+cornerPoints[i][j][0]);
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

		// dont delete localSizes data, it is stored in the rank 0 geometry object
		/*for (index_t i=0; i<tdim[0]; i++){
			delete[] localSizes[i];
		}
		delete[] localSizes;*/
		_localSizes = localSizes;

		// dont delete cornerPoints data, it is stored in the rank 0 geometry object
		/*for (index_t i=0; i<tdim[0]; i++){
			delete[] cornerPoints[i];
		}
		delete[] cornerPoints;*/
		// store it
		_offsets = cornerPoints;

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
		// does something like "do_domain_decomposition(multi_index_t(0,0), nullptr, nullptr);"
		int** dummya;
		multi_index_t** dummyb;
		multi_index_t dummyc;
		do_domain_decomposition(dummyc, dummya, dummyb);

		MPI_Recv(&totalOffsetX, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&totalOffsetY, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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

	_total_offset[0] = totalOffsetX;
	_total_offset[1] = totalOffsetY;

	// dont delete local storage, because the member pointer points to this data
	
	// output for debugging purpurse	
	//if(_comm->getRank() == 0) testIterator();

	//output_flags();
	//std::cout << "Rank: " << _comm->getRank() << ", total size: (" << _bsize[0] << ", " << _bsize[1] << "), " << "local size: (" << _size[0] << ", " << _size[1] << "), total offset: (" << _total_offset[0] << ", " << _total_offset[1] << ")\n" << std::flush;
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

const multi_index_t& Geometry::TotalOffset() const
{
	return _total_offset;
}

const multi_index_t& Geometry::Offset(index_t i, index_t j) const
{
	if (_comm->getRank() != 0){
		// the data is undefined!
		std::cout << "Fatal error! Geometry: Offset data is undefined on non-master process!\n" << std::flush;
		throw std::runtime_error("Geometry: Offset data is undefined on non-master process.");
	} else {
		return _offsets[i][j];
	}
}

const multi_index_t& Geometry::LocalSize(index_t i, index_t j) const
{
	if (_comm->getRank() != 0){
		// the data is undefined!
		std::cout << "Warning: Geometry: LocalSize data is undefined on non-master process!\n" << std::flush;
		throw std::runtime_error("Geometry: LocalSize data is undefined on non-master process.");
	} else {
		return _localSizes[i][j];
	}
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
#ifdef twocell_criterion_check
			if(!isObstacle(it.Left()) && !isObstacle(it.Right())) {
				std::cout << "Warning in u(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in x-direction)!!!\n" << std::flush;
			} else if (!isObstacle(it.Top()) && !isObstacle(it.Down())) {
				std::cout << "Warning in u(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in y-direction)!!!\n" << std::flush;
			} else 
#endif
			if (leftright) {
				// the corner cases don't have to be considered as they were already consider above!
				if (!isObstacle(it.Left())) {
					u->Cell(it) = u->Cell(it.Left()) + _h[0]*bvalU(it);
				} else if (!isObstacle(it.Right())) {
					u->Cell(it) = u->Cell(it.Right()) - _h[0]*bvalU(it);
				}
			} else if (topdown) {
				//either at the top or the bottom cell is fluid
				//if (!isObstacle(it.Left())) {
				//	u->Cell(it) = u->Cell(it.Left()) + _h[0]*bvalU(it);
				//} else if (!isObstacle(it.Right())) {
				//	u->Cell(it) = u->Cell(it.Right()) - _h[0]*bvalU(it);
				//} else {
					if (!isObstacle(it.Top())) {
						u->Cell(it) = u->Cell(it.Top()) - _h[1]*bvalU(it);
					} else if (!isObstacle(it.Down())) {
						u->Cell(it) = u->Cell(it.Down()) + _h[1]*bvalU(it);
					}
				//}
			}
		} else {
#ifdef twocell_criterion_check
			if(!isObstacle(it.Left()) && !isObstacle(it.Right())) {
				std::cout << "Warning in u(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in x-direction)!!!\n" << std::flush;
			} else if (!isObstacle(it.Top()) && !isObstacle(it.Down())) {
				std::cout << "Warning in u(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in y-direction)!!!\n" << std::flush;
			} else 
#endif
			if (topdown) {
				//either at the top or the bottom cell is fluid
				if (!isObstacle(it.Left())) {
					u->Cell(it) = bvalU(it);
					u->Cell(it.Left()) = bvalU(it);
				} else if (!isObstacle(it.Right())) {
					u->Cell(it) = bvalU(it);
					//u->Cell(it) = 0.0; // TODO: Remove
				} else {
					if (!isObstacle(it.Top())) {
						u->Cell(it) = 2.0*bvalU(it) - u->Cell(it.Top());
						u->Cell(it.Left()) = 2.0*bvalU(it.Left()) - u->Cell(it.Left().Top());
					} else if (!isObstacle(it.Down())) {
						u->Cell(it) = 2.0*bvalU(it) - u->Cell(it.Down());
						u->Cell(it.Left()) = 2.0*bvalU(it.Left()) - u->Cell(it.Left().Down());
					}
				}
			} else if (leftright) {
				// the corner cases don't have to be considered as they were already consider above!
				if (!isObstacle(it.Left())) {
					u->Cell(it) = bvalU(it);
					u->Cell(it.Left()) = bvalU(it);
				} else if (!isObstacle(it.Right())) {
					u->Cell(it) = bvalU(it);
				}
			}
		}
		it.Next();
	}
}

void Geometry::UpdateGG_V(Grid *v) const
{
	BoundaryIteratorGG it(this);
	it.First();

	while(it.Valid()) {
 		// check where the fluid is
		bool topdown = !isObstacle(it.Top()) || !isObstacle(it.Down());
		bool leftright = !isObstacle(it.Left()) || !isObstacle(it.Right());
		if(isNeumannBoundaryV(it)) {
#ifdef twocell_criterion_check
			if(!isObstacle(it.Left()) && !isObstacle(it.Right())) {
				std::cout << "Warning in v(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in x-direction)!!!\n" << std::flush;
			} else if (!isObstacle(it.Top()) && !isObstacle(it.Down())) {
				std::cout << "Warning in v(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in y-direction)!!!\n" << std::flush;
			} else 
#endif
			if (topdown) {
				//either at the top or the bottom cell is fluid
				//if (!isObstacle(it.Left())) {
				//	v->Cell(it) = v->Cell(it.Left()) + _h[0]*bvalV(it);
				//} else if (!isObstacle(it.Right())) {
				//	v->Cell(it) = v->Cell(it.Right()) - _h[0]*bvalV(it);
				//} else {
					if (!isObstacle(it.Top())) {
						v->Cell(it) = v->Cell(it.Top()) - _h[1]*bvalV(it);
					} else if (!isObstacle(it.Down())) {
						v->Cell(it) = v->Cell(it.Down()) + _h[1]*bvalV(it);
					}
				//}
			} else if (leftright) {
				// the corner cases don't have to be considered as they were already consider above!
				if (!isObstacle(it.Left())) {
					v->Cell(it) = v->Cell(it.Left()) + _h[0]*bvalV(it);
				} else if (!isObstacle(it.Right())) {
					v->Cell(it) = v->Cell(it.Right()) - _h[0]*bvalV(it);
				}
			}
		} else {
#ifdef twocell_criterion_check
			if(!isObstacle(it.Left()) && !isObstacle(it.Right())) {
				std::cout << "Warning in v(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in x-direction)!!!\n" << std::flush;
			} else if (!isObstacle(it.Top()) && !isObstacle(it.Down())) {
				std::cout << "Warning in v(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in y-direction)!!!\n" << std::flush;
			} else 
#endif
			if (leftright) {
				//either at the top or the bottom cell is fluid
				if (!isObstacle(it.Top())) {
					v->Cell(it) = bvalV(it);
				} else if (!isObstacle(it.Down())) {
					v->Cell(it) = bvalV(it);
					v->Cell(it.Down()) = bvalV(it);
				} else {
					if (!isObstacle(it.Left())) {
						v->Cell(it) = 2.0*bvalV(it) - v->Cell(it.Left());
						v->Cell(it.Down()) = 2.0*bvalV(it.Down()) - v->Cell(it.Down().Left());
					} else if (!isObstacle(it.Right())) {
						v->Cell(it) = 2.0*bvalV(it) - v->Cell(it.Right());
						v->Cell(it.Down()) = 2.0*bvalV(it.Down()) - v->Cell(it.Down().Right());
					}
				}
			} else if (topdown) {
				// the corner cases don't have to be considered as they were already considered above!
				if (!isObstacle(it.Top())) {
					v->Cell(it) = bvalV(it);
				} else if (!isObstacle(it.Down())) {
					v->Cell(it) = bvalV(it);
					v->Cell(it.Down()) = bvalV(it);
				}
			}
		}
		it.Next();
	}
}

void Geometry::UpdateGG_P(Grid *p) const
{
	BoundaryIteratorGG it(this);
	it.First();

	while(it.Valid()) {
 		// check where the fluid is
		bool topdown = !isObstacle(it.Top()) || !isObstacle(it.Down());
		bool leftright = !isObstacle(it.Left()) || !isObstacle(it.Right());
		if(isNeumannBoundaryP(it)) {
#ifdef twocell_criterion_check
			if(!isObstacle(it.Left()) && !isObstacle(it.Right())) {
				std::cout << "Warning in pN(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in x-direction)!!!\n" << std::flush;
				std::cout << it.Left().Pos()[0] << ", " << it.Left().Pos()[1] << "| " << it.Right().Pos()[0] << ", " << it.Right().Pos()[1] << std::endl;
			} else if (!isObstacle(it.Top()) && !isObstacle(it.Down())) {
				std::cout << "Warning in pN(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in y-direction)!!!\n" << std::flush;
			} else 
#endif
			if (topdown) {
				//either at the top or the bottom cell is fluid
				if (!isObstacle(it.Left())) {
					if (!isObstacle(it.Top())) {

						//std::cout << "P Neumann: LeftTop: " << p->Cell(it) << ", " << 0.5*( p->Cell(it.Left()) + p->Cell(it.Top()) ) << ", " << p->Cell(it.Left()) << ", " << p->Cell(it.Top()) << "\n";

						//std::cout << "P Neumann: left top\n";

						p->Cell(it) = 0.5*(p->Cell(it.Top()) - _h[1]*bvalP(it) + p->Cell(it.Left()) + _h[0]*bvalP(it));
					} else if (!isObstacle(it.Down())) {
						p->Cell(it) = 0.5*(p->Cell(it.Down()) + _h[1]*bvalP(it) + p->Cell(it.Left()) + _h[0]*bvalP(it));

						//std::cout << "P Neumann: left down\n";

					} /*else {
						while(true) std::cout << "Warning!\n";
					}*/
				} else if (!isObstacle(it.Right())) {
					if (!isObstacle(it.Top())) {

						//std::cout << "P Neumann: RightTop: " << p->Cell(it) << ", " << 0.5*( p->Cell(it.Right()) + p->Cell(it.Top()) ) << ", " << p->Cell(it.Right()) << ", " << p->Cell(it.Top()) << "\n";
						//if (abs(p->Cell(it)) >= 1.0) std::cout << "Attention, value of p (right top corner), it position: " << it.Pos()[0] << ", " << it.Pos()[1] << ", val p: " << p->Cell(it) << "\n";
						//if (abs(p->Cell(it.Right())) >= 1.0) std::cout << "Attention, value of p right, it position: " << it.Right().Pos()[0] << ", " << it.Right().Pos()[1] << ", val p: " << p->Cell(it.Right()) << "\n";
						//if (abs(p->Cell(it.Top())) >= 1.0) std::cout << "Attention, value of p top, it position: " << it.Top().Pos()[0] << ", " << it.Top().Pos()[1] << ", val p: " << p->Cell(it.Top()) << "\n";

						//std::cout << "P Neumann: right top\n";

						p->Cell(it) = 0.5*(p->Cell(it.Top()) - _h[1]*bvalP(it) + p->Cell(it.Right()) - _h[0]*bvalP(it));

					} else if (!isObstacle(it.Down())) {
						p->Cell(it) = 0.5*(p->Cell(it.Down()) + _h[1]*bvalP(it) + p->Cell(it.Right()) - _h[0]*bvalP(it));

						//std::cout << "P Neumann: right down\n";

					} else {
						std::cout << "Warning! Wrong boundary condition interpretation!\n";
					}
				} else {
					if (!isObstacle(it.Top())) {
						//std::cout << "P Neumann: top\n";
						p->Cell(it) = p->Cell(it.Top()) - _h[1]*bvalP(it);
					} else if (!isObstacle(it.Down())) {
						//std::cout << "P Neumann: down\n";
						p->Cell(it) = p->Cell(it.Down()) + _h[1]*bvalP(it);
					} else {
						std::cout << "Warning! Wrong boundary condition interpretation!\n";
					}
				}
			} else if (leftright) {
				// the corner cases don't have to be considered as they were already consider above!
				if (!isObstacle(it.Left())) {
					//std::cout << "P Neumann: left\n";
					p->Cell(it) = p->Cell(it.Left()) + _h[0]*bvalP(it);
				} else if (!isObstacle(it.Right())) {
					//std::cout << "P Neumann: right\n";
					p->Cell(it) = p->Cell(it.Right()) - _h[0]*bvalP(it);
				}
			}
		} else {
#ifdef twocell_criterion_check
			if(!isObstacle(it.Left()) && !isObstacle(it.Right())) {
				std::cout << "Warning in pD(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in x-direction)!!!\n" << std::flush;
				std::cout << it.Left().Pos()[0] << ", " << it.Left().Pos()[1] << "| " << it.Right().Pos()[0] << ", " << it.Right().Pos()[1] << std::endl;
			} else if (!isObstacle(it.Top()) && !isObstacle(it.Down())) {
				std::cout << "Warning in pD(" << it.Pos()[0] << ", " << it.Pos()[1] << "): The obstacle is too thin (in y-direction)!!!\n" << std::flush;
			} else 
#endif
			if (topdown) {
				//either at the top or the bottom cell is fluid
				if (!isObstacle(it.Left())) {
					if (!isObstacle(it.Top())) {
						//std::cout << "P Dirichlet: left top\n";
						p->Cell(it) = 2.0*bvalP(it) - 0.5*(p->Cell(it.Left()) + p->Cell(it.Top()));
					} else if (!isObstacle(it.Down())) {
						//std::cout << "P Dirichlet: left down\n";
						p->Cell(it) = 2.0*bvalP(it) - 0.5*(p->Cell(it.Left()) + p->Cell(it.Down()));
					}
				} else if (!isObstacle(it.Right())) {
					if (!isObstacle(it.Top())) {
						//std::cout << "P Dirichlet: right top\n";
						p->Cell(it) = 2.0*bvalP(it) - 0.5*(p->Cell(it.Right()) + p->Cell(it.Top()));
					} else if (!isObstacle(it.Down())) {
						//std::cout << "P Dirichlet: right down\n";
						p->Cell(it) = 2.0*bvalP(it) - 0.5*(p->Cell(it.Right()) + p->Cell(it.Down()));
					}
				} else {
					if (!isObstacle(it.Top())) {
						//std::cout << "P Dirichlet: top\n";
						p->Cell(it) = 2.0*bvalP(it) - p->Cell(it.Top());
					} else if (!isObstacle(it.Down())) {
						//std::cout << "P Dirichlet: down\n";
						p->Cell(it) = 2.0*bvalP(it) - p->Cell(it.Down());
					}
				}
			} else if (leftright) {
				// the corner cases don't have to be considered as they were already consider above!
				if (!isObstacle(it.Left())) {
					//std::cout << "P Dirichlet: left\n";
					p->Cell(it) = 2.0*bvalP(it) - p->Cell(it.Left());
				} else if (!isObstacle(it.Right())) {
					//std::cout << "P Dirichlet: right\n";
					p->Cell(it) = 2.0*bvalP(it) - p->Cell(it.Right());
				}
			}
		}
		it.Next();
	}
}

/*void Geometry::UpdateGG_F(Grid* F, Grid* u) const
{
	BoundaryIteratorGG it(this);
	it.First();

	while(it.Valid()) {
		if (!isObstacle(it.Left())) {
			F->Cell(it.Left()) = u->Cell(it.Left());
			F->Cell(it) = u->Cell(it);
		} else {
			F->Cell(it) = u->Cell(it);
		}
		it.Next();
	}
}

void Geometry::UpdateGG_G(Grid* G, Grid* v) const
{
	BoundaryIteratorGG it(this);
	it.First();

	while(it.Valid()) {
		if (!isObstacle(it.Down())) {
			G->Cell(it.Down()) = v->Cell(it.Down());
			G->Cell(it) = v->Cell(it);
		} else {
			G->Cell(it) = v->Cell(it);
		}
		it.Next();
	}
}*/

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

#ifdef OUTPUT_GEOMETRY
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
#endif
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

	// add ghost cells to domains
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
	//return (_flags[it.Value()] >> 0) & 1;
	return _flags[it.Value()] & 1;

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
	//bool boundary((_flags[it.Value()] >> 3) & 1);
	//if (_comm->getRank() == 1) std::cout << "Boundary Value for P(" << _comm->getRank() << "): " << boundary << " Position: " << it.Pos()[0] << ", " << it.Pos()[1] << std::endl;
	return (_flags[it.Value()] >> 3) & 1;
}

// Getter functions for the boundary data
const real_t& Geometry::bvalU(const Iterator& it) const
{
	// std::cout << "Position: " << it.Pos()[0] << ", " << it.Pos()[1] << "| bval_u: " << _bval_u[it.Value()] << std::endl;
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


void Geometry::output_flags() const
{
	fprintf(stderr,"====================Output of Flags====================\n");
	Iterator it(this);
	it.First();
	while (it.Valid()) {
		fprintf(stderr,"%i ", int(_flags[it.Value()]));
		if(it.Right().Value() == it.Value())
			fprintf(stderr,"\n");
		it.Next();
	}
	fprintf(stderr,"====================End of Output!====================\n\n");

	fprintf(stderr,"====================Output of U====================\n");
	Iterator it1(this);
	it1.First();
	while (it1.Valid()) {
		fprintf(stderr,"%.2f ", _bval_u[it1.Value()]);
		if(it1.Right().Value() == it1.Value())
			fprintf(stderr,"\n");
		it1.Next();
	}
	fprintf(stderr,"====================End of Output!====================\n\n");

	fprintf(stderr,"====================Output of V====================\n");
	Iterator it2(this);
	it2.First();
	while (it2.Valid()) {
		fprintf(stderr,"%.2f ", _bval_v[it2.Value()]);
		if(it2.Right().Value() == it2.Value())
			fprintf(stderr,"\n");
		it2.Next();
	}
	fprintf(stderr,"====================End of Output!====================\n\n");

	fprintf(stderr,"====================Output of P====================\n");
	Iterator it3(this);
	it3.First();
	while (it3.Valid()) {
		fprintf(stderr,"%.2f ", _bval_p[it3.Value()]);
		if(it3.Right().Value() == it3.Value())
			fprintf(stderr,"\n");
		it3.Next();
	}
	fprintf(stderr,"====================End of Output!====================\n\n");
}

/*void Geometry::output_flags() const
{
	for (int i = 0; i < _bsize[0]*_bsize[1]; i++)
	{
		std::cout << int(_flags[i]) << " " << std::flush;
	}
	std::cout << std::endl;
}*/

void Geometry::testIterator() const
{
	Grid* tmp = new Grid(this);
	tmp->Initialize(0.0);

	/*InteriorIteratorGG it1(this);
	it1.First();

	while(it1.Valid()) {
		tmp->Cell(it1) += 1;
		it1.Next();
	}

	BoundaryIteratorGG it2(this);
	it2.First();

	while(it2.Valid()) {
		tmp->Cell(it2) += 2;
		it2.Next();
	}*/

	JumpingInteriorIteratorGG it1(this, true);
	it1.First();

	while(it1.Valid()) {
		tmp->Cell(it1) += 1;
		it1.Next();
	}

	JumpingInteriorIteratorGG it2(this, false);
	it2.First();

	while(it2.Valid()) {
		tmp->Cell(it2) += 2;
		it2.Next();
	}

	tmp->Out();
	delete tmp;	
}

bool Geometry::isInsideObstacle(const multi_real_t& pos, const multi_real_t& offset) const
{
	real_t ix = ( (Size()[0] - 2.0)*(pos[0]/Length()[0]) ) - offset[0];
	real_t iy = ( (Size()[1] - 2.0)*(pos[1]/Length()[1]) ) - offset[1];
	
//	if (_offset[0] > 0 || _offset[1] > 0)
//		std::cout << "Warning: Positive Offset \n";

	index_t x = round(ix);
	index_t y = round(iy);

	index_t val = x + y * Size()[0];

	return isObstacle(Iterator(this, val));
}

bool Geometry::isInsideThisSubdomain(const multi_real_t& pos) const
{
	real_t minX = _blength[0] * real_t(_total_offset[0]) / real_t(_bsize[0]-2);
	real_t maxX = minX + _blength[0] * real_t(_size[0]-2) / real_t(_bsize[0]-2);
	real_t minY = _blength[1] * real_t(_total_offset[1]) / real_t(_bsize[1]-2);
	real_t maxY = minY + _blength[1] * real_t(_size[1]-2) / real_t(_bsize[1]-2);

	return pos[0] >= minX && pos[0] <= maxX && pos[1] >= minY && pos[1] <= maxY;
}

/*void Geometry::homogeneousBoundary()
{
	
}*/

void Geometry::fitToGeom(const Geometry* geom)
{
	_size = geom->_size;
	_bsize = geom->_bsize;
	_length = geom->_length;
	_blength = geom->_blength;

	set_meshwidth();
}

void Geometry::setSize(multi_index_t size)
{
	_size = size;

	set_meshwidth();
}

