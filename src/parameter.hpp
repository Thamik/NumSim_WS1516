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
#ifndef __PARAMETER_HPP
#define __PARAMETER_HPP
//------------------------------------------------------------------------------
/// Class providing the parameters needed for the simulation
class Parameter {
public:
  /// Constructs a new Parameter set with default values
  // Driven Cavity parameters; see exercise sheet 1
  Parameter();

  /// Loads the parameter values from a file
  void Load(const char *file, bool verbose = false);

  /// Getter functions for all parameters
  const real_t &Re() const;
  const real_t &Omega() const;
  const real_t &Alpha() const;
  const real_t &Dt() const;
  const real_t &Tend() const;
  const index_t &IterMax() const;
  const real_t &Eps() const;
  const real_t &Tau() const;

private:
  real_t _re; // Reynolds number
  real_t _omega; // relaxation factor
  real_t _alpha; // upwind differencing factor
  real_t _dt; // time step size
  real_t _tend; // end time
  real_t _eps; // tolerance for pressure iteration
  real_t _tau; // safety factor time step size
  index_t _itermax; // maximum number of iterations
};
//------------------------------------------------------------------------------
#endif // __PARAMETER_HPP
