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

Communicator::Communicator(int *argc, char ***argv)
: _tidx(0,0), _tdim(0,0), _rank(0), _size(0), _evenodd(false), _rankDistribution(NULL)
{
	//for (int i = 0; i<*argc; i++) std::cout << (*argv)[i] << "\n" << std::flush;
	MPI_Init(argc, argv); //TODO
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank); // determine rank of process
	MPI_Comm_size(MPI_COMM_WORLD, &_size); // determine number of processes
	//if (_rank == 0) std::cout << "MPI initialized with " << _size << " processes.\n" << std::flush;
	//TODO Noch was?
}

Communicator::~Communicator()
{
	MPI_Finalize();
	if( _rankDistribution != NULL ) {
		for(int ii = 0; ii < _tdim[0]; ii++) {
			delete[] _rankDistribution[ii];
		}
		delete[] _rankDistribution;
	}
	//TODO More to do?
}

const multi_index_t& Communicator::ThreadIdx() const
{
	//TODO Correct?
	return _tidx;	
}

const multi_index_t& Communicator::ThreadDim() const
{
	//TODO Correct?
	return _tdim;
}

const bool &Communicator::EvenOdd() const
{
	return _evenodd;
}

real_t Communicator::gatherSum(const real_t &val) const
{
	real_t res(0.0);
	real_t sendBuff(val);
	MPI_Allreduce(&sendBuff, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return res;
}

real_t Communicator::gatherMin(const real_t &val) const
{
	real_t res(0.0);
	real_t sendBuff(val);
	MPI_Allreduce(&sendBuff, &res, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return res;
}

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

const bool Communicator::isLeft() const
{
	return _tidx[0] == 0;
}

const bool Communicator::isRight() const
{
	return _tidx[0] == _tdim[0] - 1;
}

const bool Communicator::isTop() const
{
	return _tidx[1] == _tdim[1] - 1;
}


const bool Communicator::isBottom() const
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
	// delete _rankDistribution if not NULL
	if( _rankDistribution != NULL ) {
		for(int ii = 0; ii < _tdim[0]; ii++) {
			delete[] _rankDistribution[ii];
		}
		delete[] _rankDistribution;
	}

	_tdim = tdim;

	// allocate _rankDistribution
	_rankDistribution = new int*[_tdim[0]];
	for(int ii = 0; ii < _tdim[0]; ii++) {
		_rankDistribution[ii] = new int[_tdim[1]];
	}

	// write new values in _rankDistribution
	for(int ii = 0; ii < _tdim[0]; ii++) {
		for(int jj = 0; jj < _tdim[1]; jj++) {
			_rankDistribution[ii][jj] = int(rankDistri[ii][jj]);
		}
	}

	for(int ii = 0; ii < _tdim[0]; ii++) {
		for(int jj = 0; jj < _tdim[1]; jj++) {
			if( _rankDistribution[ii][jj] == _rank ) {
				_tidx[0] = ii;
				_tidx[1] = jj;
			}
		}
	}

	_localSize = localSizes[_tidx[0]][_tidx[1]];

	// set _evenodd
	int temp(0);
	for (int i=0; i<_tidx[0]; i++){
		temp += localSizes[i][0][0];
	}
	for (int i=0; i<_tidx[1]; i++){
		temp += localSizes[i][0][1];
	}
	if ((temp % 2) == 0){
		_evenodd = true;
	} else {
		_evenodd = false;
	}
}


/* Private methods */

bool Communicator::copyLeftBoundary(Grid *grid) const
{
	int tag = 1013;

	//TODO Debug
	const Geometry* tempGeom = grid->getGeometry();
	real_t* sendBuff;
	real_t* recBuff;

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
		MPI_Sendrecv(sendBuff, _localSize[1], MPI_DOUBLE, leftRank, tag, recBuff, _localSize[1], MPI_DOUBLE, rightRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if(!isLeft()) {
		int leftRank = _rankDistribution[_tidx[0]-1][_tidx[1]];
		MPI_Send(sendBuff, _localSize[1], MPI_DOUBLE, leftRank, tag, MPI_COMM_WORLD);
	}
	else if(!isRight()) {
		int rightRank = _rankDistribution[_tidx[0]+1][_tidx[1]];
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

	if(!isLeft()) delete[] sendBuff;
	if(!isRight()) delete[] recBuff;

	return true;
}

bool Communicator::copyRightBoundary(Grid *grid) const
{
	int tag = 1013;

	//TODO Debug
	const Geometry* tempGeom = grid->getGeometry();
	real_t* sendBuff;
	real_t* recBuff;

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
		MPI_Sendrecv(sendBuff, _localSize[1], MPI_DOUBLE, rightRank, tag, recBuff, _localSize[1], MPI_DOUBLE, leftRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if(!isRight()) {
		int rightRank = _rankDistribution[_tidx[0]+1][_tidx[1]];
		MPI_Send(sendBuff, _localSize[1], MPI_DOUBLE, rightRank, tag, MPI_COMM_WORLD);
	}
	else if(!isLeft()) {
		int leftRank = _rankDistribution[_tidx[0]-1][_tidx[1]];
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

	if(!isRight()) delete[] sendBuff;
	if(!isLeft()) delete[] recBuff;

	return true;
}

bool Communicator::copyTopBoundary(Grid *grid) const
{
	int tag = 1013;

	//TODO Debug
	const Geometry* tempGeom = grid->getGeometry();
	real_t* sendBuff;
	real_t* recBuff;

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
		MPI_Sendrecv(sendBuff, _localSize[0], MPI_DOUBLE, topRank, tag, recBuff, _localSize[0], MPI_DOUBLE, bottomRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if(!isTop()) {
		int topRank = _rankDistribution[_tidx[0]][_tidx[1]+1];
		MPI_Send(sendBuff, _localSize[0], MPI_DOUBLE, topRank, tag, MPI_COMM_WORLD);
	}
	else if(!isBottom()) {
		int bottomRank = _rankDistribution[_tidx[0]][_tidx[1]-1];
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

	if(!isTop()) delete[] sendBuff;
	if(!isBottom()) delete[] recBuff;

	return true;
}

bool Communicator::copyBottomBoundary(Grid *grid) const
{
	int tag = 1013;

	//TODO Debug
	const Geometry* tempGeom = grid->getGeometry();
	real_t* sendBuff;
	real_t* recBuff;

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
		MPI_Sendrecv(sendBuff, _localSize[0], MPI_DOUBLE, bottomRank, tag, recBuff, _localSize[0], MPI_DOUBLE, topRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if(!isBottom()) {
		int bottomRank = _rankDistribution[_tidx[0]][_tidx[1]-1];
		MPI_Send(sendBuff, _localSize[0], MPI_DOUBLE, bottomRank, tag, MPI_COMM_WORLD);
	}
	else if(!isTop()) {
		int topRank = _rankDistribution[_tidx[0]][_tidx[1]+1];
		MPI_Recv(recBuff, _localSize[0], MPI_DOUBLE, topRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else {
		// the process is at the left as well as at the right boundary
		//TODO maybe output for debugging perpurse?
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

	if(!isBottom()) delete[] sendBuff;
	if(!isTop()) delete[] recBuff;

	return true;
}


