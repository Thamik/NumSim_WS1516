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
  Geometry(const Communicator *comm);

  /// Loads a geometry from a file
  void Load(const char *file, bool verbose = false);

  /// Returns the number of cells in each dimension
  const multi_index_t &Size() const;
  /// Returns the total number of cells in each dimension
  const multi_index_t &TotalSize() const;
  /// Returns the length of the domain
  const multi_real_t &Length() const;
  /// Returns the total length of the domain
  const multi_real_t &TotalLength() const;
  /// Returns the meshwidth
  const multi_real_t &Mesh() const;

	/// Returns a pointer to the communicator
	const Communicator* getCommunicator() const;

  /// Updates the velocity field u
  void Update_U(Grid *u) const;
  /// Updates the velocity field v
  void Update_V(Grid *v) const;
  /// Updates the pressure field p
  void Update_P(Grid *p) const;

private:
  const Communicator *_comm;

  /// size of the grid
  multi_index_t _size;
  multi_index_t _bsize;
  /// physical length of the domain
  multi_real_t _length;
  multi_real_t _blength;
  /// meshwidth
  multi_real_t _h;

  /// velocity
  multi_real_t _velocity;
  /// pressure
  real_t _pressure;

	// own method
	/// sets the meshwidth
	void set_meshwidth();

	/// decides whether the given boundary is a global boundary or is not
	bool is_global_boundary(int boundary_index) const;
	
};
//------------------------------------------------------------------------------
#endif // __GEOMETRY_HPP
