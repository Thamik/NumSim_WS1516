#include "grid.hpp"

/* public methods */

/* constructor */
Grid::Grid(const Geometry *geom)
{
	// TODO
}

Grid::Grid(const Geometry *geom, const multi_real_t& offset)
{
	// TODO
}

/* destructor */
Grid::~Grid()
{
	// TODO
}

void Grid::Initialize(const real_t& value)
{
	// TODO
}

real_t& Grid::Cell(const Iterator &it)
{
	// TODO
}

const real_t& Grid::Cell(const Iterator &it) const
{
	// TODO
}

real_t Grid::Interpolate(const multi_real_t& pos) const
{
	// TODO
}

real_t Grid::dx_l(const Iterator& it) const
{
	// TODO
}

real_t Grid::dx_r(const Iterator& it) const
{
	// TODO
}

real_t Grid::dy_l(const Iterator& it) const
{
	// TODO
}

real_t Grid::dy_r(const Iterator& it) const
{
	// TODO
}

real_t Grid::dxx(const Iterator &it) const
{
	// TODO
}

real_t Grid::dyy(const Iterator& it) const
{
	// TODO
}

real_t Grid::DC_udu_x(const Iterator& it, const real_t& alpha) const
{
	// TODO
}

real_t Grid::DC_vdu_y(const Iterator& it, const real_t& alpha, const Grid* v) const
{
	// TODO
}

real_t Grid::DC_udv_x(const Iterator& it, const real_t& alpha, const Grid* u) const
{
	// TODO
}

real_t Grid::DC_vdv_y(const Iterator& it, const real_t& alpha) const
{
	// TODO
}

real_t Grid::Max() const
{
	// TODO
}

real_t Grid::Min() const
{
	// TODO
}

real_t Grid::AbsMax() const
{
	// TODO
}

real_t* Grid::Data()
{
	// TODO
}

