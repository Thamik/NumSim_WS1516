#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"

#include <stdio.h>      /* NULL */
#include <stdlib.h>     /* malloc, free */

/* public methods */

/* constructor */
Grid::Grid(const Geometry* geom, const multi_real_t& offset)
{
	_geom = geom;
	_data = (real_t*) malloc(_geom->Size()[0] * _geom->Size()[1] * sizeof(real_t));
	if (_data==NULL) exit(-1);
	// TODO: what about the offset?
	_offset = offset; // is this right?
}

Grid::Grid(const Geometry* geom)
{
	Grid(geom, 0.0); // does this make sense (offset=0)?
}

/* destructor */
Grid::~Grid()
{
	free(_data);
	// TODO
}

void Grid::Initialize(const real_t& value)
{
	// TODO: test this method
	for(int i=0; i<_geom->Size()[0]*_geom->Size()[1]; i++)
	{
		_data[i] = value;
	}
}

real_t& Grid::Cell(const Iterator& it)
{
	// TODO: test
	return _data[it.Value()];
}

const real_t& Grid::Cell(const Iterator& it) const
{
	// TODO: test
	return _data[it.Value()];
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
