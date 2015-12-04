#include "iterator.hpp"
#include "geometry.hpp"

#include <iostream>

// Iterate through the grid like this:
/*
y
^
|...
|-------------------------->
|-------------------------->
|-------------------------->
|___________________________> x

Hence, from left to right, and then from the bottom to the top.
*/

/* Iterator */

/**
The constructor with the Iterator class.
\param[in]	geom	Geometry object, that provides all geometry information
\param[in]	value	Starting value of the iterator
*/
Iterator::Iterator(const Geometry* geom, const index_t& value)
: _geom(geom), _value(value), _valid(true)
{
}

/**
Constructor for the Iterator class.
\param[in]	geom	Geometry object, that provides all geometry information
*/
Iterator::Iterator(const Geometry* geom)
: Iterator(geom, 0) // default initial position is zero
{
}

/**
\return The value of the iterator, that is its current position in the grid.
*/
const index_t& Iterator::Value() const
{
	return _value;
}

/**
This is a cast opertor to convert an iterator to an index type. Returns the value of the iterator.
/return The value of the iterator
*/
Iterator::operator const index_t&() const
{
	return _value;
}

/** This method calculates the current position
of the Iterator in the cartesian x-y-coordinate
system.
\return The x and y coordinates
*/
multi_index_t Iterator::Pos() const
{
	index_t x, y;
	x = _value % _geom->Size()[0];
	y = _value / _geom->Size()[0];
	return multi_index_t(x,y);
}

/** The Iterator is set to the first element.
Here, the first element is the element located
in the bottom left corner
*/
void Iterator::First()
{
	_value = 0;
	//_valid = true; // maybe do this?
}

/**
The Iterator is set to the next element.
Here, the next element is the element to the right or the
most left element of the above line, respectively (cf. Iterationpath).
If the Iterator exits the field, it is disabled!
*/
void Iterator::Next()
{
	_value++;
	if (_value >= _geom->Size()[0]*_geom->Size()[1]) _valid = false;
}

/**
\return status of the Iterator
*/
bool Iterator::Valid() const
{
	return _valid;
}

/**
If the Iterator is currently at the left boundary,
it returns itself
\return Iterator to the Left
*/
Iterator Iterator::Left() const
{
	if (Pos()[0] == 0){
		// at the left boundary, return a copy of self
		return Iterator(_geom, _value);
	}else{
		// somewhere else, return the cell to the left
		return Iterator(_geom, _value - 1);
	}
}

/**
If the Iterator is currently at the right boundary,
it returns itself.
\return Iterator to the right
*/
Iterator Iterator::Right() const
{
	if (Pos()[0] == _geom->Size()[0]-1){
		// at the right boundary, return a copy of self
		return Iterator(_geom, _value);
	}else{
		// somewhere else, return the cell to the right
		return Iterator(_geom, _value + 1);
	}
}

/**
If the Iterator is currently at the top boundary,
it returns itself.
\return above Iterator
*/
Iterator Iterator::Top() const
{
	if (Pos()[1] == _geom->Size()[1]-1){
		// at the upper boundary, return a copy of self
		return Iterator(_geom, _value);
	}else{
		// somewhere else, return the cell above
		return Iterator(_geom, _value + _geom->Size()[0]);
	}
}

/**
If the Iterator is currently at the bottom boundary,
it returns itself.
\return lower Iterator
*/
Iterator Iterator::Down() const
{
	if (Pos()[1] == 0){
		// at the lower boundary, return a copy of self
		return Iterator(_geom, _value);
	}else{
		// somewhere else, return the cell below
		return Iterator(_geom, _value - _geom->Size()[0]);
	}
}

/* InteriorIterator */

/**
The constructor of the InteriorIterator class.
\param[in]	geom	Geometry object containing all geometrical data
*/
InteriorIterator::InteriorIterator(const Geometry* geom)
: Iterator(geom)
{
}

/**
Sets the position of the InteriorIterator to the first interior cell, which is the second cell from below and the second from left.
*/
void InteriorIterator::First()
{
	_value = _geom->Size()[0] + 1; // the second cell from below and the second from left
	// maybe set _valid as true?
}

/**
Moves the InteriorIterator to the next interior cell. If there are no further cells, the InteriorIterator is set to be invalid.
*/
void InteriorIterator::Next()
{
	if (Pos()[0] >= _geom->Size()[0]-2){
		_value += 3; // jump over the right and then the left boundary cells
	} else {
		_value++;
	}
	if (Pos()[1] >= _geom->Size()[1]-1) _valid = false; // maybe do this earlier?

	//std::cout << Pos()[0] << ", " << Pos()[1] << "\n";
}

//------------------------------------------------------------------------------

/* JumpingInteriorIterator */

/**
This function is the constructor of the JumpingInteriorIterator class.
\param[in]	geom	a Geometry object of the grids for which the iterator may be used
\param[in]	shifted	a boolean which specifies if the lower left (interior) cell of the grid, over which will be iterated, is a red or black cell (see RedOrBlackSOR)
*/
JumpingInteriorIterator::JumpingInteriorIterator(const Geometry* geom, bool shifted)
: InteriorIterator(geom), _shifted(shifted)
{
}

/**
This function sets the position of the iterator to the first valid cell.
*/
void JumpingInteriorIterator::First()
{
	InteriorIterator::First();
	if (_shifted) InteriorIterator::Next();
}

/**
This function moves the iterator to the next valid cell.
*/
void JumpingInteriorIterator::Next()
{
	if ((_geom->Size()[0] % 2) == 1){
		InteriorIterator::Next();
		InteriorIterator::Next();
	} else {
		// handle row change
		index_t temp_row = Pos()[1];
		InteriorIterator::Next();
		if (Pos()[1] == temp_row){
			InteriorIterator::Next();
			if (Pos()[1] != temp_row){
				InteriorIterator::Next();
			}
		}
	}
}

//------------------------------------------------------------------------------

/* BoundaryIterator */

/**
The constructor of the BoundaryIterator class. The boundary which shall be considered is undefined. To set the boundary, call the SetBoundary() method.
\param[in]	geom	Geometry object with all geometrical information
*/
BoundaryIterator::BoundaryIterator(const Geometry *geom)
: Iterator(geom), _boundary(0)
{
}

/**
This method sets the index of the boundary, which shall be considered, e.i. over which shall be iterated.
\param[in]	boundary	An index type value between 1 and 4. 1 for the left, 2 for the right, 3 for the bottom and 4 for the top boundary.
*/
void BoundaryIterator::SetBoundary(const index_t& boundary)
{
	_boundary = boundary;
}

/**
Moves the iterator to the first cell of the boundary specified with the method SetBoundary()
*/
void BoundaryIterator::First()
{
	if (_boundary == 0){
		// _boundary is unset, dont know what to do here
		std::cout << "Warning, BoundaryIterator: boundary index seems to not be set!\n";
	} else if (_boundary == 1){
		// the left boundary, set to the lower left corner cell
		_value = 0;
	} else if (_boundary == 2){
		// the right boundary, set to the lower right corner cell
		_value = _geom->Size()[0] - 1;
	} else if (_boundary == 3){
		// the lower boundary, set to the lower left corner cell
		_value = 0;
	} else if (_boundary == 4){
		// the upper boundary, set to the upper left corner cell
		_value = _geom->Size()[0] * (_geom->Size()[1] - 1);
	} else {
		// this should not happen, invalid value in _boundary
		std::cout << "Warning, BoundaryIterator: invalid boundary index!\n";
	}
	// maybe set _valid as true?
}

/**
Moves the BoundaryIterator to the next cell of the specified boundary. If there are no further cells, the BoundaryIterator is set to be invalid.
*/
void BoundaryIterator::Next()
{
	index_t temp(0);

	if (_boundary == 0){
		// _boundary is unset, dont know what to do here
	} else if (_boundary == 1){
		// the left boundary
		temp = Top().Value();
	} else if (_boundary == 2){
		// the right boundary
		temp = Top().Value();
	} else if (_boundary == 3){
		// the lower boundary	
		temp = Right().Value();
	} else if (_boundary == 4){
		// the upper boundary
		temp = Right().Value();
	} else {
		// this should not happen, invalid value in _boundary
	}

	if (temp == _value){
		// we were and are at the end (we stayed there)
		_valid = false;
	} else {
		_value = temp;
	}
}


//------------------------------------------------------------------------------
// InteriorIteratorGG

void InteriorIteratorGG::InteriorIteratorGG(const Geometry* geom) 
: InteriorIterator(geom)
{
}

void InteriorIteratorGG::First()
{
	InteriorIterator::First();
}

void InteriorIteratorGG::Next()
{
	InteriorIterator::Next();
	while (geom->isObstacle(this)) {
		InteriorIterator::Next();
	}
}

//------------------------------------------------------------------------------
// InteriorIteratorGG
void JumpingInteriorIteratorGG::JumpingInteriorIteratorGG(const Geometry* geom, bool shifted)
: JumpingInteriorIterator(geom, shifted)
{
}

void JumpingInteriorIteratorGG::First()
{
	JumpingInteriorIterator::First();
}

void JumpingInteriorIteratorGG::Next()
{
	JumpingInteriorIterator::Next();
	while (geom->isObstacle(this)) {
		JumpingInteriorIterator::Next();
	}
}

//------------------------------------------------------------------------------
// BoundaryIteratorGG
void BoundaryIteratorGG::BoundaryIteratorGG(const Geometry* geom)
: Iterator(geom)
{
}

void BoundaryIteratorGG::First()
{
	Iterator::First();
}

void BoundaryIteratorGG::Next()
{
	Iterator::Next();
	while (!geom->isObstacle(this)) {
		Iterator::Next();
	}

}
