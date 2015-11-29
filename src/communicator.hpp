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

//------------------------------------------------------------------------------
#include "grid.hpp"
#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __COMMUNICATOR_HPP
#define __COMMUNICATOR_HPP
//------------------------------------------------------------------------------
/// Class for the communication between the different processes
class Communicator {
public:
  /// Constructor
  Communicator(int *argc, char ***argv, bool verbose = false);

  /// Destructor
  ~Communicator();

  /// requests the position of the current process
  const multi_index_t &ThreadIdx() const;

  /// Returns the way the domain is partitioned among all processes
  const multi_index_t &ThreadDim() const;

  /// Returns whether this process is a red or a black field
  const bool &EvenOdd() const;

  /// Gets the sum of all values and distributes the result among all processes
  real_t gatherSum(const real_t &val) const;

  /// Finds the minimum of the values and distributes the result among all processes
  real_t gatherMin(const real_t &val) const;

  /// Finds the maximum of the values and distributes the result among all processes
  real_t gatherMax(const real_t &val) const;

  /// Synchronizes ghost layer
  void copyBoundary(Grid *grid, bool verbose = false) const;

  /// Decide whether our left boundary is a domain boundary
  bool isLeft() const;

  /// Decide whether our right boundary is a domain boundary
  bool isRight() const;

  /// Decide whether our top boundary is a domain boundary
  bool isTop() const;

  /// Decide whether our bottom boundary is a domain boundary
  bool isBottom() const;

  /// Get MPI rank of current process
  const int &getRank() const;

  /// Get number of MPI processes
  const int &getSize() const;


	// own methods
	/// Writes important values after domain decomposition
	void setProcDistribution(int** rankDistri, multi_index_t tdim, multi_index_t** localSizes, bool verbose = false);

	/// Returns the size of the grid handled by the current process
	multi_index_t getLocalSize() const;

	/// Returns the Rank of the process assigned to the (i,j)-cell in the process grid
	int getRankDistribution(int i, int j) const;


private:
  /// position of the current process
  multi_index_t _tidx;
  /// dimensions of the process grid
  multi_index_t _tdim;
  /// rank of this process
  int _rank; 
  /// number of processes
  int _size; 
  /// determines if the grids belonging to this process start with a red or black cell
  bool _evenodd; 
  /// the distribution of the processes to the decomposed domain
  int** _rankDistribution;
  /// the size of the (local) domain belonging to this process
  multi_index_t _localSize; 


  /* Function to sync ghost layer on left boundary:
   *  send values of own left boundary to left neighbor and
   *  and receive values from his right boundary
   *
   *   ------------ ------------
   *  |           x|y           |
   *  |           x|y           |
   *  |           x|y           |
   *  |           x|y           |
   *  |           x|y           |
   *   ------------ ------------
   *
   *   y: values that are sent
   *   x: values that are received
   *
   * 
   */
  /// Function to sync ghost layer on left boundary:
  bool copyLeftBoundary(Grid *grid) const;

  /// Function to sync ghost layer on right boundary
  bool copyRightBoundary(Grid *grid) const;

  /// Function to sync ghost layer on top boundary
  bool copyTopBoundary(Grid *grid) const;

  /// Function to sync ghost layer on bottom boundary
  bool copyBottomBoundary(Grid *grid) const;
};
//------------------------------------------------------------------------------
#endif // __COMMUNICATOR_HPP
//------------------------------------------------------------------------------
