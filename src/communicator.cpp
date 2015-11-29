#include "communicator.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "solver.hpp"

#include <cmath>	// sin, M_PI
#include <algorithm>    // std::min
#include <iostream>	// std::cout
#include <mpi.h>	// MPI

/* Public methods */

/** Communicator constructor; initializes MPI Environment
   *
   * \param [in] argc Number of arguments program was started with
   * \param [in] argv Arguments passed to the program on start
   * \param [in] verbose Debug output on terminal
   */
Communicator::Communicator(int *argc, char ***argv, bool verbose)
: _tidx(0,0), _tdim(0,0), _rank(0), _size(0), _evenodd(false), _rankDistribution(NULL), _localSize(0,0)
{
	if (verbose && _rank == 0) {
		std::cout << "MPI initializing..." << std::flush;
	}
	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank); // determine rank of process
	MPI_Comm_size(MPI_COMM_WORLD, &_size); // determine number of processes
	if (verbose && _rank == 0) {
		std::cout << "MPI initialized with " << _size << " processes.\n" << std::flush;
	}
}

/** Communicator destructor; finalizes MPI Environment
   */
Communicator::~Communicator()
{
	MPI_Finalize();
	if( _rankDistribution != NULL ) {
		for(index_t ii = 0; ii < _tdim[0]; ii++) {
			delete[] _rankDistribution[ii];
		}
		delete[] _rankDistribution;
	}
}

/** Returns the position of the current process with respect to the
   *  fields lower left corner
   *  \return position of current process within the process grid
   */
const multi_index_t& Communicator::ThreadIdx() const
{
	return _tidx;	
}

/**
	\return x- and y- sizes of process grid
*/
const multi_index_t& Communicator::ThreadDim() const
{
	return _tdim;
}

/**
*/
const bool &Communicator::EvenOdd() const
{
	return _evenodd;
}

/**  This is done by using the MPI-function MPI_Allreduce
   *
   * 	\param [in] val The data over which the sum is to be calculated
   *	\return gathered summation over all processes
   */
real_t Communicator::gatherSum(const real_t &val) const
{
	real_t res(0.0);
	real_t sendBuff(val);
	MPI_Allreduce(&sendBuff, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return res;
}

/** This is done by using the MPI-function MPI_Allreduce
   *
   * \param [in] val The data over which to find the minimum
   * \return gathered minimum over all processes
   */
real_t Communicator::gatherMin(const real_t &val) const
{
	real_t res(0.0);
	real_t sendBuff(val);
	MPI_Allreduce(&sendBuff, &res, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return res;
}

/** This is done by using the MPI-function MPI_Allreduce
	\param[in] val The data over which to find the maximum
	\return gathered maximum over all processes
*/
real_t Communicator::gatherMax(const real_t &val) const
{
	real_t res(0.0);
	real_t sendBuff(val);
	MPI_Allreduce(&sendBuff, &res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	return res;
}

void Communicator::copyBoundary(Grid *grid) const
{
	bool success(true);
	bool temp(true);

	// left boundary
	temp = copyLeftBoundary(grid);
	success = success && temp;

	// right boundary
	temp = copyRightBoundary(grid);
	success = success && temp;

	// top boundary
	temp = copyTopBoundary(grid);
	success = success && temp;
	
	// bottom boundary
	temp = copyBottomBoundary(grid);
	success = success && temp;
	
	/*// output for debugging
	if( success ) {
		std::cout << "Process no. " << _rank << " copied successful!\n";
	}
	else {
		std::cout << "Copy failed!!! (Process no. " << _rank << ")\n";
	}*/
}

bool Communicator::isLeft() const
{
	return _tidx[0] == 0;
}

bool Communicator::isRight() const
{
	return _tidx[0] == _tdim[0] - 1;
}

bool Communicator::isTop() const
{
	return _tidx[1] == _tdim[1] - 1;
}


bool Communicator::isBottom() const
{
	return _tidx[1] == 0;
}

const int &Communicator::getRank() const
{
	return _rank;
}

const int &Communicator::getSize() const
{
	return _size;
}

multi_index_t Communicator::getLocalSize() const
{
	return _localSize;
}

void Communicator::setProcDistribution(int** rankDistri, multi_index_t tdim, multi_index_t** localSizes)
{
	// std::cout << _rank << ": " << tdim[0] << ", " << tdim[1] << "\n" << std::flush;

	// delete _rankDistribution if not NULL
	if( _rankDistribution != NULL ) {
		for(index_t ii = 0; ii < _tdim[0]; ii++) {
			delete[] _rankDistribution[ii];
		}
		delete[] _rankDistribution;
	}

	_tdim = tdim;

	//std::cout << _rank << ": " << "start allocating..." << "\n" << std::flush;
	// allocate _rankDistribution
	_rankDistribution = new int*[_tdim[0]];
	for(index_t ii = 0; ii < _tdim[0]; ii++) {
		_rankDistribution[ii] = new int[_tdim[1]];
	}

	//std::cout << _rank << ": " << "done allocating!" << "\n" << std::flush;

	// write new values in _rankDistribution
	for(index_t ii = 0; ii < _tdim[0]; ii++) {
		for(index_t jj = 0; jj < _tdim[1]; jj++) {
			_rankDistribution[ii][jj] = int(rankDistri[ii][jj]);
		}
	}

	//std::cout << _rank << ": " << "start searching for position..." << "\n" << std::flush;
	for(index_t ii = 0; ii < _tdim[0]; ii++) {
		for(index_t jj = 0; jj < _tdim[1]; jj++) {
			if( _rankDistribution[ii][jj] == _rank ) {
				_tidx[0] = ii;
				_tidx[1] = jj;
			}
		}
	}
	//std::cout << _rank << ": " << "position found!" << "\n" << std::flush;

	_localSize = localSizes[_tidx[0]][_tidx[1]];

	//std::cout << _rank << ": " << "local sizes written!" << "\n" << std::flush;

	// set _evenodd
	int temp(0);
	for (index_t i=0; i<_tidx[0]; i++){
		temp += localSizes[i][0][0];
	}
	for (index_t i=0; i<_tidx[1]; i++){
		temp += localSizes[0][i][1];
	}
	if ((temp % 2) == 0){
		_evenodd = true;
	} else {
		_evenodd = false;
	}

	//std::cout << _rank << ": " << "evenodd set!" << "\n" << std::flush;
}


/* Private methods */

bool Communicator::copyLeftBoundary(Grid *grid) const
{
	int tag = 1013;

	//TODO Debug
	const Geometry* tempGeom = grid->getGeometry();
	real_t* sendBuff(NULL);
	real_t* recBuff(NULL);

	//initialize boundary iterator
	BoundaryIterator itSend(tempGeom);
	itSend.SetBoundary(BoundaryIterator::boundaryLeft);
	itSend.First();

	//TODO create buffer as local variables!!!!
	// copy boundary values to buffer
	if(!isLeft()) {
		sendBuff = new real_t[_localSize[1]];
		int index(0);
		while(itSend.Valid()) {
			sendBuff[index] = grid->Cell(itSend.Right());
			index++;
			itSend.Next();
		}
	}

	//allocate recieve buffer
	if(!isRight()) recBuff = new real_t[_localSize[1]];
	
	//Blocking SendReceive
	if(!isRight() && !isLeft()) {
		int leftRank = _rankDistribution[_tidx[0]-1][_tidx[1]];
		int rightRank = _rankDistribution[_tidx[0]+1][_tidx[1]];

		if( leftRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << leftRank << std::flush;
		if( rightRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << rightRank << std::flush;

		MPI_Sendrecv(sendBuff, _localSize[1], MPI_DOUBLE, leftRank, tag, recBuff, _localSize[1], MPI_DOUBLE, rightRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if(!isLeft()) {
		int leftRank = _rankDistribution[_tidx[0]-1][_tidx[1]];

		if( leftRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << leftRank << std::flush;
		
		MPI_Send(sendBuff, _localSize[1], MPI_DOUBLE, leftRank, tag, MPI_COMM_WORLD);
	}
	else if(!isRight()) {
		int rightRank = _rankDistribution[_tidx[0]+1][_tidx[1]];
		
		if( rightRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << rightRank << std::flush;

		MPI_Recv(recBuff, _localSize[1], MPI_DOUBLE, rightRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else {
		// the process is at the left as well as at the right boundary
		//TODO maybe output for debugging perpurse?
	}

	// Copy received data to grid
	BoundaryIterator itRecv(tempGeom);
	itRecv.SetBoundary(BoundaryIterator::boundaryRight);
	itRecv.First();
	
	if(!isRight()) {
		int index(0);
		while(itRecv.Valid()) {
			grid->Cell(itRecv) = recBuff[index];
			index++;
			itRecv.Next();
		}
	}

	if(!isLeft() && sendBuff!=NULL) delete[] sendBuff;
	if(!isRight() && recBuff!=NULL) delete[] recBuff;

	return true;
}

bool Communicator::copyRightBoundary(Grid *grid) const
{
	int tag = 1014;

	//TODO Debug
	const Geometry* tempGeom = grid->getGeometry();
	real_t* sendBuff(NULL);
	real_t* recBuff(NULL);

	//initialize boundary iterator
	BoundaryIterator itSend(tempGeom);
	itSend.SetBoundary(BoundaryIterator::boundaryRight);
	itSend.First();

	//TODO create buffer as local variables!!!!
	// copy boundary values to buffer
	if(!isRight()) {
		sendBuff = new real_t[_localSize[1]];
		int index(0);
		while(itSend.Valid()) {
			sendBuff[index] = grid->Cell(itSend.Left());
			index++;
			itSend.Next();
		}
	}

	//allocate recieve buffer
	if(!isLeft()) recBuff = new real_t[_localSize[1]];
	
	//Blocking SendReceive
	if(!isRight() && !isLeft()) {
		int leftRank = _rankDistribution[_tidx[0]-1][_tidx[1]];
		int rightRank = _rankDistribution[_tidx[0]+1][_tidx[1]];

		if( leftRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << leftRank << std::flush;
		if( rightRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << rightRank << std::flush;

		MPI_Sendrecv(sendBuff, _localSize[1], MPI_DOUBLE, rightRank, tag, recBuff, _localSize[1], MPI_DOUBLE, leftRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if(!isRight()) {
		int rightRank = _rankDistribution[_tidx[0]+1][_tidx[1]];

		if( rightRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << rightRank << std::flush;

		MPI_Send(sendBuff, _localSize[1], MPI_DOUBLE, rightRank, tag, MPI_COMM_WORLD);
	}
	else if(!isLeft()) {
		int leftRank = _rankDistribution[_tidx[0]-1][_tidx[1]];

		if( leftRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << leftRank << std::flush;

		MPI_Recv(recBuff, _localSize[1], MPI_DOUBLE, leftRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else {
		// the process is at the left as well as at the right boundary
		//TODO maybe output for debugging perpurse?
	}

	// Copy received data to grid
	BoundaryIterator itRecv(tempGeom);
	itRecv.SetBoundary(BoundaryIterator::boundaryLeft);
	itRecv.First();
	
	if(!isLeft()) {
		int index(0);
		while(itRecv.Valid()) {
			grid->Cell(itRecv) = recBuff[index];
			index++;
			itRecv.Next();
		}
	}

	if(!isRight() && sendBuff!=NULL) delete[] sendBuff;
	if(!isLeft() && recBuff!=NULL) delete[] recBuff;

	return true;
}

bool Communicator::copyTopBoundary(Grid *grid) const
{
	int tag = 1015;

	//TODO Debug
	const Geometry* tempGeom = grid->getGeometry();
	real_t* sendBuff(NULL);
	real_t* recBuff(NULL);

	//initialize boundary iterator
	BoundaryIterator itSend(tempGeom);
	itSend.SetBoundary(BoundaryIterator::boundaryTop);
	itSend.First();

	//TODO create buffer as local variables!!!!
	// copy boundary values to buffer
	if(!isTop()) {
		sendBuff = new real_t[_localSize[0]];
		int index(0);
		while(itSend.Valid()) {
			sendBuff[index] = grid->Cell(itSend.Down());
			index++;
			itSend.Next();
		}
	}

	//allocate recieve buffer
	if(!isBottom()) recBuff = new real_t[_localSize[0]];
	
	//Blocking SendReceive
	if(!isBottom() && !isTop()) {
		int topRank = _rankDistribution[_tidx[0]][_tidx[1]+1];
		int bottomRank = _rankDistribution[_tidx[0]][_tidx[1]-1];

		if( topRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << topRank << std::flush;
		if( bottomRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << bottomRank << std::flush;

		MPI_Sendrecv(sendBuff, _localSize[0], MPI_DOUBLE, topRank, tag, recBuff, _localSize[0], MPI_DOUBLE, bottomRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if(!isTop()) {
		int topRank = _rankDistribution[_tidx[0]][_tidx[1]+1];

		if( topRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << topRank << std::flush;

		MPI_Send(sendBuff, _localSize[0], MPI_DOUBLE, topRank, tag, MPI_COMM_WORLD);
	}
	else if(!isBottom()) {
		int bottomRank = _rankDistribution[_tidx[0]][_tidx[1]-1];
		
		if( bottomRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << bottomRank << std::flush;

		MPI_Recv(recBuff, _localSize[0], MPI_DOUBLE, bottomRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else {
		// the process is at the left as well as at the right boundary
		//TODO maybe output for debugging perpurse?
	}

	// Copy received data to grid
	BoundaryIterator itRecv(tempGeom);
	itRecv.SetBoundary(BoundaryIterator::boundaryBottom);
	itRecv.First();
	
	if(!isBottom()) {
		int index(0);
		while(itRecv.Valid()) {
			grid->Cell(itRecv) = recBuff[index];
			index++;
			itRecv.Next();
		}
	}

	if(!isTop() && sendBuff!=NULL) delete[] sendBuff;
	if(!isBottom() && recBuff!=NULL) delete[] recBuff;

	return true;
}

bool Communicator::copyBottomBoundary(Grid *grid) const
{
	int tag = 1016;

	//TODO Debug
	const Geometry* tempGeom = grid->getGeometry();
	real_t* sendBuff(NULL);
	real_t* recBuff(NULL);

	//initialize boundary iterator
	BoundaryIterator itSend(tempGeom);
	itSend.SetBoundary(BoundaryIterator::boundaryBottom);
	itSend.First();

	//TODO create buffer as local variables!!!!
	// copy boundary values to buffer
	if(!isBottom()) {
		sendBuff = new real_t[_localSize[0]];
		int index(0);
		while(itSend.Valid()) {
			sendBuff[index] = grid->Cell(itSend.Top());
			index++;
			itSend.Next();
		}
	}

	//allocate recieve buffer
	if(!isTop()) recBuff = new real_t[_localSize[0]];
	
	//Blocking SendReceive
	if(!isBottom() && !isTop()) {
		int topRank = _rankDistribution[_tidx[0]][_tidx[1]+1];
		int bottomRank = _rankDistribution[_tidx[0]][_tidx[1]-1];
		
		if( topRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << topRank << std::flush;
		if( bottomRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << bottomRank << std::flush;
		
		MPI_Sendrecv(sendBuff, _localSize[0], MPI_DOUBLE, bottomRank, tag, recBuff, _localSize[0], MPI_DOUBLE, topRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if(!isBottom()) {
		int bottomRank = _rankDistribution[_tidx[0]][_tidx[1]-1];

		if( bottomRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << bottomRank << std::flush;

		MPI_Send(sendBuff, _localSize[0], MPI_DOUBLE, bottomRank, tag, MPI_COMM_WORLD);
	}
	else if(!isTop()) {
		int topRank = _rankDistribution[_tidx[0]][_tidx[1]+1];

		if( topRank >= _size ) std::cout << "process " << _rank << " tried to call rank " << topRank << std::flush;

		MPI_Recv(recBuff, _localSize[0], MPI_DOUBLE, topRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else {
		// the process is at the left as well as at the right boundary
		//TODO maybe output for debugging perpurse?
		//std::cout << "process " << _rank << "is bottom and top" << std::flush;
	}

	// Copy received data to grid
	BoundaryIterator itRecv(tempGeom);
	itRecv.SetBoundary(BoundaryIterator::boundaryTop);
	itRecv.First();
	
	if(!isTop()) {
		int index(0);
		while(itRecv.Valid()) {
			grid->Cell(itRecv) = recBuff[index];
			index++;
			itRecv.Next();
		}
	}

	if(!isBottom() && sendBuff!=NULL) delete[] sendBuff;
	if(!isTop() && recBuff!=NULL) delete[] recBuff;

	return true;
}

int Communicator::getRankDistribution(int i, int j) const
{
	return _rankDistribution[i][j];
}
