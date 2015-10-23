#include "iterator.hpp"
#include "geometry.hpp"

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

/* Constructor */
Iterator::Iterator(const Geometry* geom, const index_t& value)
: _geom(geom), _value(value), _valid(true)
{
	// TODO: is this correct?
}

Iterator::Iterator(const Geometry* geom)
{
	Iterator(geom, 0); // Setting the default initial value (position) to zero
}

const index_t& Iterator::Value() const
{
	return _value;
}

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
	// TODO: test
	index_t x, y;
	x = _value / _geom->Size()[0];
	y = _value % _geom->Size()[0];
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
	if (_value >= _geom->Size()[0]*_geom->Size()[1]) _valid = false; // maybe one step earlier?
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
	// TODO: test
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
	// TODO: test
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
	// TODO: test
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
	// TODO: test
	if (Pos()[1] == 0){
		// at the lower boundary, return a copy of self
		return Iterator(_geom, _value);
	}else{
		// somewhere else, return the cell below
		return Iterator(_geom, _value - _geom->Size()[0]);
	}
}

/* not needed
// Protected Help Functions
index_t Iterator::Pos_To_Value(multi_index_t pos) const
{
	// TODO
}
*/

// InteriorIterator

// Constructor
InteriorIterator::InteriorIterator(const Geometry* geom)
: Iterator(geom)
{
	// TODO: constructor complete?
}

void InteriorIterator::First()
{
	// TODO: right?
	_value = _geom->Size()[0] + 1; // the second cell from below and the second from left
	// maybe set _valid as true?
}

void InteriorIterator::Next()
{
	// TODO: test
	if (Pos()[0] >= _geom->Size()[0]-1){
		_value += 3; // jump over the right and then the left boundary cells
	} else {
		_value++;
	}
	if (Pos()[1] >= _geom->Size()[1]-1) _valid = false; // maybe do this earlier?
}

// BoundaryIterator

// Constructor
BoundaryIterator::BoundaryIterator(const Geometry *geom)
: Iterator(geom), _boundary(0) // is this right?
{
	// TODO: constructor complete?
}

void BoundaryIterator::SetBoundary(const index_t& boundary)
{
	// TODO: test!
	_boundary = boundary;
	// handling boundary as counter only for boundary cells

	// setting _value
	if (_boundary <= _geom->Size()[0]-1){
		// at the lower boundary
		_value = _boundary;
	} else if (_boundary >= _geom->Size()[0] + 2*(_geom->Size()[1]-2)){
		// at the upper boundary
		_value = _geom->Size()[0]*(_geom->Size()[1]-1) + _boundary - _geom->Size()[0] + 2*(_geom->Size()[1]-2);
	} else {
		// at the left or right boundary, not at the lower, not at the upper
		if ( (_boundary-_geom->Size()[0])%2 == 0){
			_value = _geom->Size()[0] * ((_boundary-_geom->Size()[0])/2 + 1);
		} else {
			_value = _geom->Size()[0] * ((_boundary-_geom->Size()[0])/2 + 1) + _geom->Size()[0] - 1;
		}
	}
	// maybe set _valid as true?
}

void BoundaryIterator::First()
{
	// TODO: is this right?
	_value = 0;
	_boundary = 0;
	// maybe set _valid as true?
}

void BoundaryIterator::Next()
{
	// TODO: is this right?
	SetBoundary(_boundary+1);
	if (_boundary >= 2*_geom->Size()[0] + 2*_geom->Size()[1] - 4) _valid = false; // maybe do this earlier?
}
