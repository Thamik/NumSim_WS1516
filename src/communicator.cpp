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
{
	MPI_Init( argc, argv); //TODO
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank); // determine rank of process
	//TODO
}

Communicator::~Communicator()
{
	MPI_Finalize();
	//TODO More to do?
}

const multi_index_t& Communicator::ThreadIdx() const
{
	//TODO
}

const multi_index_t& Communicator::ThreadDim() const
{
	//TODO
}

const bool &Communicator::EvenOdd() const
{
	return _evenodd;
	//TODO Correct?
}

real_t Communicator::gatherSum(const real_t &val) const
{
	//TODO With Reduction MPI_SUM
}

real_t Communicator::gatherMin(const real_t &val) const
{
	//TODO With Reduction and MPI_MIN
}

real_t Communicator::gatherMax(const real_t &val) const
{
	//TODO With Reduction and MPI_MAX
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
	//TODO
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


