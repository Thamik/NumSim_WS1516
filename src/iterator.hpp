/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __ITERATOR_HPP
#define __ITERATOR_HPP
//------------------------------------------------------------------------------
/// Iterator base class
/** Iterate through the grid like this:<br>
<code>
y<br>
^<br>
|...<br>
|--------------------------><br>
|--------------------------><br>
|--------------------------><br>
|___________________________> x<br>
</code>
Hence, from left to right, and then from the bottom to the top.
*/
class Iterator {
public:
  /// Constructs a new Iterator depending on a geometry
  Iterator(const Geometry *geom);
  /// Constructs a new Iterator on a geometry with a defined starting value
  Iterator(const Geometry *geom, const index_t &value);

	Iterator(const Geometry* geom, const multi_index_t& pos);

  /// Returns the current position value
  virtual const index_t &Value() const;
  /// Cast operator to convert Iterators to integers
  virtual operator const index_t &() const;
  /// Returns the (index) position coordinates
  virtual multi_index_t Pos() const;

	/// Returns the physical position coordinates
	virtual multi_real_t PhysPos(const multi_real_t& offset) const;

  /// Sets the iterator to the first element
  virtual void First();
  /// Goes to the next element of the iterator, disables it if position is end
  virtual void Next();

  /// Checks if the iterator still has a valid value
  virtual bool Valid() const;

  /// Returns an Iterator that is located left from this one.
  // If we are at the left boundary, the cell sees itself
  virtual Iterator Left() const;

  /// Returns an Iterator that is located right from this one
  // If we are at the right boundary, the cell sees itself
  virtual Iterator Right() const;

  /// Returns an Iterator that is located above this one
  // If we are at the upper domain boundary, the cell sees itself
  virtual Iterator Top() const;

  /// Returns an Iterator that is located below this one
  // If we are at the lower domain boundary, the cell sees itself
  virtual Iterator Down() const;

protected:
  const Geometry *_geom;
  index_t _value;
  bool _valid;

	multi_index_t _pos;

	Iterator(const Geometry* geom, const index_t &value, const multi_index_t& pos);

};

//------------------------------------------------------------------------------
/// Iterator for interior cells
class InteriorIterator : public Iterator {
public:
  /// Construct a new InteriorIterator
  InteriorIterator(const Geometry* geom);

  /// Sets the iterator to the first element
  virtual void First();
  /// Goes to the next element of the iterator, disables it if position is end
  virtual void Next();
};

//------------------------------------------------------------------------------
/// Iterator for interior cells in the red or black cycle
class JumpingInteriorIterator : public InteriorIterator {
public:
	/// Construct a new JumpingInteriorIterator
	JumpingInteriorIterator(const Geometry* geom, bool shifted);

	/// Sets the iterator to the first element
	void First();
	/// Goes to the next element of the iterator, disables it if position is end
	void Next();

private:
	bool _shifted;
};

//------------------------------------------------------------------------------
/// Iterator for domain boundary cells.
class BoundaryIterator : public Iterator {
public:

  // TODO: use this enum
  /*enum {
    boundaryBottom = 0,
    boundaryLeft = 1,
    boundaryTop = 2,
    boundaryRight = 3
  };*/
  enum { // our values
    boundaryBottom = 3,
    boundaryLeft = 1,
    boundaryTop = 4,
    boundaryRight = 2
  };

  /// Constructs a new BoundaryIterator
  BoundaryIterator(const Geometry *geom);

  /// Sets the boundary to iterate
  void SetBoundary(const index_t &boundary);

  /// Sets the iterator to the first element
  void First();
  /// Goes to the next element of the iterator, disables it if position is end
  void Next();

private:
  index_t _boundary;
	// _boundary enumerates the side of the area:
	// left = 1, right = 2, bottom = 3, top = 4
	// 0 = unset
};

//------------------------------------------------------------------------------
/// New InteriorIterator for General Geometries
class InteriorIteratorGG : public InteriorIterator {
public:
  /// Construct a new InteriorIterator
  InteriorIteratorGG(const Geometry* geom);

  /// Sets the iterator to the first element
  virtual void First();
  /// Goes to the next element of the iterator, disables it if position is end
  virtual void Next();
};


//------------------------------------------------------------------------------
/// New JumpingInteriorIterator for General Geometries
class JumpingInteriorIteratorGG : public JumpingInteriorIterator {
public:
  /// Constructs a new JumpingInteriorIteratorGG
  JumpingInteriorIteratorGG(const Geometry* geom, bool shifted);

  /// Sets the iterator to the first element
  virtual void First();
  /// Goes to the next element of the iterator, disables it if position is end
  virtual void Next();

private:
  bool _shifted;
};

//------------------------------------------------------------------------------
/// New JumpingInteriorIterator for General Geometries
class BoundaryIteratorGG : public Iterator {
public:
  /// Construct a new InteriorIterator
  BoundaryIteratorGG(const Geometry* geom);

  /// Sets the iterator to the first element
  virtual void First();
  /// Goes to the next element of the iterator, disables it if position is end
  virtual void Next();
};
//------------------------------------------------------------------------------
#endif // __ITERATOR_HPP
