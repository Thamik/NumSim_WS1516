#include "geometry.hpp"

/* public methods */

/* constructor */
Geometry::Geometry()
: _size(128,128), _length(1.0,1.0), _velocity(1.0,0.0), _pressure(1.0)
{
	_h = multi_real_t(_length[0]/_size[0],_length[1]/_size[0]);
	// TODO: are the values right?
}

void Geometry::Load(const char *file)
{
	// TODO
}

const multi_index_t& Geometry::Size() const
{
	return _size;
}

const multi_real_t& Geometry::Length() const
{
	return _length;
}

const multi_real_t& Geometry::Mesh() const
{
	return _h;
}

void Geometry::Update_U(Grid *u) const
{
	// TODO
}

void Geometry::Update_V(Grid *v) const
{
	// TODO
}

void Geometry::Update_P(Grid *p) const
{
	// TODO
}
