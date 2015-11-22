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

/* Public methods */

Communicator::Communicator(int *argc, char ***argv)
{
	//TODO
}

Communicator::~Communicator()
{
	//TODO
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
	//TODO
}

real_t Communicator::gatherSum(const real_t &val) const
{
	//TODO
}

real_t Communicator::gatherMin(const real_t &val) const
{
	//TODO
}

real_t Communicator::gatherMax(const real_t &val) const
{
	//TODO
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
	//TODO
}

const int &Communicator::getSize() const
{
	//TODO
}


/* Private methods */

bool copyLeftBoundary(Grid *grid) const
{
	//TODO
}

bool copyRightBoundary(Grid *grid) const
{
	//TODO
}

bool copyTopBoundary(Grid *grid) const
{
	//TODO
}

bool copyBottomBoundary(Grid *grid) const
{
	//TODO
}


