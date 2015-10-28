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
#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __SOLVER_HPP
#define __SOLVER_HPP
//------------------------------------------------------------------------------

/** abstract base class for an iterative solver
*/
class Solver {
public:
  /// Constructor of the abstract Solver class
  Solver(const Geometry *geom);
  /// Destructor of the Solver Class
  virtual ~Solver();

  /// This function must be implemented in a child class
  // @param [in][out] grid current values
  // @param [in]      rhs  right hand side values
  // @returns accumulated residual
  virtual real_t Cycle(Grid *grid, const Grid *rhs) const = 0;

	// own method
	// deletes the average value
	virtual void delete_average(Grid* grid) const;

protected:
  const Geometry *_geom;

  /// Returns the residual at [it] for the pressure-Poisson equation
  real_t localRes(const Iterator &it, const Grid *grid, const Grid *rhs) const;

	// own method
	// return the total residual for the pressure-Poisson equation
	real_t totalRes(const Grid* grid, const Grid* rhs) const;
};

//------------------------------------------------------------------------------

/** concrete SOR solver
*/
class SOR : public Solver {
public:
  /// Constructs an actual SOR solver
  SOR(const Geometry *geom, const real_t &omega);
  /// Destructor
  ~SOR();

  /// Returns the total residual and executes a solver cycle
  // @param grid current pressure values
  // @param rhs right hand side
  real_t Cycle(Grid *grid, const Grid *rhs) const;

protected:
  real_t _omega;
};
//------------------------------------------------------------------------------

class JacobiSolver : public Solver {
public:
	JacobiSolver(const Geometry* geom);
	~JacobiSolver();

	real_t Cycle(Grid* grid, const Grid* rhs) const;
};


class HeatConductionSolver : public Solver {
public:
	HeatConductionSolver(const Geometry* geom);
	~HeatConductionSolver();

	real_t Cycle(Grid* grid, const Grid* rhs) const;

};


#endif // __SOLVER_HPP
