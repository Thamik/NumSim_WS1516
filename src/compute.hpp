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
#include "console_output.hpp"

#include <string>

#define GO_FAST
//#define USE_PARTICLES
//#define OUTPUT
#define RUN_SERIAL

//------------------------------------------------------------------------------
#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP
//------------------------------------------------------------------------------
/// Class where the simulation data is calculated
class Compute {
public:
  /// Creates a compute instance with given geometry and parameter
  Compute(const Geometry *geom, const Parameter *param,
          const Communicator *comm = 0, const char* uq_filename = "", int sim_id = 0);
  /// Deletes all grids
  ~Compute();

  /// Execute one time step of the fluid simulation (with or without debug info)
  // @ param printInfo print information about current solver state (residual
  // etc.)
  void TimeStep(bool verbose = false);

  /// Returns the simulated time in total
  const real_t &GetTime() const;

  /// Returns the pointer to U
  const Grid *GetU() const;
  /// Returns the pointer to V
  const Grid *GetV() const;
  /// Returns the pointer to P
  const Grid *GetP() const;
  /// Returns the pointer to RHS
  const Grid *GetRHS() const;

  /// Computes and returns the absolute velocity
  const Grid *GetVelocity();
  /// Computes and returns the vorticity
  const Grid *GetVorticity();
  /// Computes and returns the stream line values
  const Grid *GetStream();

	void writeUQFile() const;

private:
  /// current timestep
  real_t _t;

  /// donor-cell diffusion condition (p. 27)
  real_t _dtlimit;

  /// limit for residual
  real_t _epslimit;

  /// velocities
  Grid *_u;
  Grid *_v;

  /// pressure
  Grid *_p;

  /// prel. vel
  Grid *_F;
  Grid *_G;

  /// right-hand side
   Grid *_rhs;

  // container for interpolating whichever values
  Grid *_tmp_velocity;
  Grid *_tmp_vorticity;
  Grid *_tmp_stream;

	// grid for the last pressure field to detect incontinuities (and other grids for the same reason)
	Grid* _p_old;
	Grid* _rhs_old;
	Grid* _F_old;
	Grid* _G_old;

	/// specifies if the solver is currently converging
	bool _solver_converging;

	real_t _diff_p;
	real_t _diff_rhs;
	real_t _diff_F;
	real_t _diff_G;

	RedOrBlackSOR *_solver;

  const Geometry *_geom;
  const Parameter *_param;
  const Communicator *_comm;

#ifdef USE_PARTICLES
	// this is only in use on process rank 0, otherwise not in use
	Particles* _particles;
#endif

  /// Compute the new velocites u,v
  void NewVelocities(const real_t &dt);
  /// Compute the temporary velocites F,G
  void MomentumEqu(const real_t &dt);
  /// Compute the RHS of the poisson equation
  void RHS(const real_t &dt);

	// own methods
	/// computing dt for stability
	real_t compute_dt() const;
	/// updates the boundary values for u,v,F,G
	void update_boundary_values();

	/// synchronize the grids
	void sync_FG();
	void sync_uv();
	void sync_p();
	void sync_all();

	/// the console clock
	ConsoleClock* _clock;

	void check_for_incontinuities();

	real_t errPipe();

	real_t breakOffPoint();

	const std::string _uq_filename;

	const int _sim_id;

};
//------------------------------------------------------------------------------
#endif // __COMPUTE_HPP
