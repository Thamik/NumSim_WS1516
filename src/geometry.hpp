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
#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP
//------------------------------------------------------------------------------
/// Class providing the geometrical information
class Geometry {
public:
  /// Constructs a default geometry:
  // driven cavity with 128 x 128 grid, no-slip boundary conditions
  // as shown below
  //
  //      u=1, v=0
  //    -------------
  //    |           |
  // u=0|           |u=0
  // v=0|           |v=0
  //    |           |
  //    |           |
  //    -------------
  //      u=0, v=0
  Geometry();

  /// Loads a geometry from a file
  void Load(const char *file, bool verbose = false);

  /// Returns the number of cells in each dimension
  const multi_index_t &Size() const;
  /// Returns the length of the domain
  const multi_real_t &Length() const;
  /// Returns the meshwidth
  const multi_real_t &Mesh() const;

  /// Updates the velocity field u
  void Update_U(Grid *u) const;
  /// Updates the velocity field v
  void Update_V(Grid *v) const;
  /// Updates the pressure field p
  void Update_P(Grid *p) const;

private:
  /// size of the grid
  multi_index_t _size;
  /// physical length of the domain
  multi_real_t _length;
  /// meshwidth
  multi_real_t _h;

  /// velocity
  multi_real_t _velocity;
  /// pressure
  real_t _pressure;

	// own method
	/// sets the meshwidth
	void set_meshwidth();
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP
