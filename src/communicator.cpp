#include "communicator.hpp"
#include "compute.hpp"
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
	MPI_Init( argc, argv); //TODO
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank); // determine rank of process
	MPI_Comm_size(MPI_COMM_WORLD, &_size); // determine number of processes
	//TODO
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
	//TODO Debug
	return _evenodd;
}

real_t Communicator::gatherSum(const real_t &val) const
{
	//TODO Debug
	real_t res(0.0);
	real_t sendBuff(val);
	MPI_Allreduce(&sendBuff, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return res;
}

real_t Communicator::gatherMin(const real_t &val) const
{
	//TODO Debug
	real_t res(0.0);
	real_t sendBuff(val);
	MPI_Allreduce(&sendBuff, &res, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return res;
}

real_t Communicator::gatherMax(const real_t &val) const
{
	//TODO Debug
	real_t res(0.0);
	real_t sendBuff(val);
	MPI_Allreduce(&sendBuff, &res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	return res;
}

void Communicator::copyBoundary(Grid *grid) const
{
	//TODO
}

const bool Communicator::isLeft() const
{
	//TODO
}

const bool Communicator::isRight() const
{

}

const bool Communicator::isTop() const
{
	//TODO
}


const bool Communicator::isBottom() const
{
	//TODO
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

void Communicator::setProcDistribution(const int** rankDistri, const multi_index_t tdim, const multi_index_t** localSizes)
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
}


/* Private methods */

bool Communicator::copyLeftBoundary(Grid *grid) const
{
	//TODO
}

bool Communicator::copyRightBoundary(Grid *grid) const
{
	//TODO
}

bool Communicator::copyTopBoundary(Grid *grid) const
{
	//TODO
}

bool Communicator::copyBottomBoundary(Grid *grid) const
{
	//TODO
}


