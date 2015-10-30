#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"

#include <stdio.h>      /* NULL */
#include <stdlib.h>     /* malloc, free */
#include <cmath>        // std::abs, pow
#include <iostream>

/* public methods */

/* constructor */
Grid::Grid(const Geometry* geom, const multi_real_t& offset)
: _geom(geom), _offset(offset)
{
	//_geom = geom;

	/*
	// TODO: kein malloc/free mehr, benutze new und delete
	_data = (real_t*) malloc(_geom->Size()[0] * _geom->Size()[1] * sizeof(real_t));
	if (_data==NULL) exit(-1);
	*/
	
	//std::cout << "Allocating memory with size " << _geom->Size()[0] << " * " << _geom->Size()[1] << "... " << std::flush; // only for debugging issues
	_data = new real_t[_geom->Size()[0] * _geom->Size()[1]];
	//std::cout << "Done.\n" << std::flush; // only for debugging issues

	// TODO: what about the offset?
	//_offset = offset; // is this right?
}

Grid::Grid(const Geometry* geom)
: Grid(geom, multi_real_t(0.0,0.0))
{
	// does this make sense (offset=0)?
}

/* destructor */
Grid::~Grid()
{
	// TODO: somethings more to do here?

	/*
	free(_data);
	*/

	delete _data;
}

void Grid::Initialize(const real_t& value)
{
	// TODO: test this method
	//std::cout << "Initializing: " << _geom->Size()[0]*_geom->Size()[1] << "\n" << std::flush;
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
	//std::cout << "Interpolate!\n" << std::flush;

	/*if(pos[0] > _offset[0] + 1 || pos[0] < _offset[0] || pos[1] > _offset[1] + 1 || pos[1] < _offset[1])
	{
		std::cout << "Warning: Interpolation - physical coordinates out of range! (" << pos[0] << ", " << pos[1] << ")\n";
		return 0.0;
	}*/
	/*real_t ix = ( (_geom->Size()[0] - 2.0)*pos[0] + 1.0 ) - (_offset[0]/_geom->Mesh()[0]);
	real_t iy = ( (_geom->Size()[1] - 2.0)*pos[1] + 1.0 ) - (_offset[1]/_geom->Mesh()[1]);*/
	real_t ix = ( (_geom->Size()[0] - 2.0)*pos[0] ) - _offset[0];
	real_t iy = ( (_geom->Size()[1] - 2.0)*pos[1] ) - _offset[1];
	if (_offset[0] > 0 || _offset[1] > 0)
		std::cout << "Warning: Positive Offset \n";
	
	//std::cout << "ix: " << ix << std::flush;
	//std::cout << "iy: " << iy << std::flush;

	index_t x1 = floor(ix);
	index_t x2 = ceil(ix);
	index_t y1 = floor(iy);
	index_t y2 = ceil(iy);

	//calculate values
	index_t vallu = x1 + y1*_geom->Size()[0];
	index_t valru = x2 + y1*_geom->Size()[0];
	index_t vallo = x1 + y2*_geom->Size()[0];	
	index_t valro = x2 + y2*_geom->Size()[0];

	//return _data[valro];

	if(vallu < 0)
		std::cout << "Warnung: " << vallu << "\n" << std::flush;
	else if(valro >= _geom->Size()[0]*_geom->Size()[1])
		std::cout << "Warnung: " << ix << ", " << iy << "\n" << std::flush;

	real_t alpha = ix - x1;
	real_t beta = iy - y1;

	//interpolation
	real_t xinter1 = (1.0 - alpha) * _data[vallu] + alpha * _data[valru];
	real_t xinter2 = (1.0 - alpha) * _data[vallo] + alpha * _data[valro];

	//std::cout << "Interpolate finished!\n" << std::flush;
	return (1.0 - beta) * xinter1 + beta * xinter2;
}

real_t Grid::dx_l(const Iterator& it) const
{
	// TODO: test
	return (_data[it.Value()] - _data[it.Left().Value()])/((_geom->Mesh())[0]);
}

real_t Grid::dx_r(const Iterator& it) const
{
	// TODO: test
	return (_data[it.Right().Value()] - _data[it.Value()])/((_geom->Mesh())[0]);
}

real_t Grid::dy_l(const Iterator& it) const
{
	// TODO: test
	return (_data[it.Value()] - _data[it.Down().Value()])/((_geom->Mesh())[1]);
}

real_t Grid::dy_r(const Iterator& it) const
{
	// TODO: test
	return (_data[it.Top().Value()] - _data[it.Value()])/((_geom->Mesh())[1]);
}

// own methods, central difference quotients
real_t Grid::dx_c(const Iterator& it) const
{
	return 0.5 * (dx_r(it)+dx_l(it));
}

real_t Grid::dy_c(const Iterator& it) const
{
	return 0.5 * (dy_r(it)+dy_l(it));
}

real_t Grid::dxx(const Iterator &it) const
{
	// TODO: test
	//return (_data[it.Right().Value()] - 2.0 * _data[it.Value()] + _data[it.Left().Value()])/( (_geom->Mesh())[0] * (_geom->Mesh())[0] );
	return (_data[it.Right().Value()] - 2.0 * _data[it.Value()] + _data[it.Left().Value()])/pow(_geom->Mesh()[0],2.0);
}

real_t Grid::dyy(const Iterator& it) const
{
	// TODO: test
	//return (_data[it.Top().Value()] - 2.0 * _data[it.Value()] + _data[it.Down().Value()])/( (_geom->Mesh())[1] * (_geom->Mesh())[1] );
	return (_data[it.Top().Value()] - 2.0 * _data[it.Value()] + _data[it.Down().Value()])/pow(_geom->Mesh()[1],2.0);
}

real_t Grid::DC_udu_x(const Iterator& it, const real_t& alpha) const
{
	// TODO: test
	// we use here that u*du/dx = 0.5 * d(u^2)/dx
	real_t res;
	real_t uij = _data[it.Value()];
	real_t uip1j = _data[it.Right().Value()];
	real_t uim1j = _data[it.Left().Value()];
	res = (1.0/(_geom->Mesh()[0])) * ( pow((uij + uip1j)/2.0, 2.0) - pow((uim1j + uij)/2.0, 2.0)) + (alpha/(_geom->Mesh()[0])) * ( std::abs(uij + uip1j)/2.0 * (uij - uip1j)/2.0 - std::abs(uim1j + uij)/2.0 * (uim1j - uij)/2.0 );
	res *= 0.5;
	return res;
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
	// TODO: test
	// we use here that v*dv/dx = 0.5 * d(v^2)/dx
	real_t res;
	real_t vij = _data[it.Value()];
	real_t vijp1 = _data[it.Top().Value()];
	real_t vijm1 = _data[it.Down().Value()];
	res = (1.0/(_geom->Mesh()[1])) * ( pow((vij + vijp1)/2.0, 2.0) - pow((vijm1 + vij)/2.0, 2.0)) + (alpha/(_geom->Mesh()[1])) * ( std::abs(vij + vijp1)/2.0 * (vij - vijp1)/2.0 - std::abs(vijm1 + vij)/2.0 * (vijm1 - vij)/2.0 );
	res *= 0.5;
	return res;
}

// die funktionen vdu_y und udv_x sind unklar, da in der Vorlesung nur d(uv)/dx und d(uv)/dy behandelt wurden. diese werden jetzt hier implementiert

/*
real_t Grid::DC_duu_x(const Iterator &it, const real_t &alpha) const
{
	return 2.0 * DC_udu_x(it,alpha);
}

real_t Grid::DC_dvv_y(const Iterator &it, const real_t &alpha) const
{
	return 2.0 * DC_vdv_y(it,alpha);
}

real_t Grid::DC_duv_x(const Iterator &it, const real_t &alpha, const Grid* u) const
{
	// TODO: test
	real_t res;
	real_t vij = _data[it.Value()];
	real_t uij = u->Data()[it.Value()];
	real_t uijp1 = u->Data()[it.Top().Value()];
	real_t vip1j = _data[it.Right().Value()];
	real_t uim1j = u->Data()[it.Left().Value()];
	real_t uim1jp1 = u->Data()[it.Left().Top().Value()];
	real_t vim1j = _data[it.Left().Value()];
	res = (1.0/_geom->Mesh()[0]) * ( (uij+uijp1)/2.0 * (vij+vip1j)/2.0 - (uim1j+uim1jp1)/2.0 * (vim1j+vij)/2.0 ) + (alpha/_geom->Mesh()[0]) * ( std::abs(uij+uijp1)/2.0 * (vij-vip1j)/2.0 - std::abs(uim1j+uim1jp1)/2.0 * (vim1j-vij)/2.0 );
	return res;
}

real_t Grid::DC_duv_y(const Iterator &it, const real_t &alpha, const Grid* v) const
{
	// TODO: test
	real_t res;
	real_t uij = _data[it.Value()];
	real_t vij = v->Data()[it.Value()];
	real_t vip1j = v->Data()[it.Right().Value()];
	real_t uijp1 = _data[it.Top().Value()];
	real_t vijm1 = v->Data()[it.Down().Value()];
	real_t vip1jm1 = v->Data()[it.Right().Down().Value()];
	real_t uijm1 = _data[it.Down().Value()];
	res = (1.0/_geom->Mesh()[1]) * ( (vij+vip1j)/2.0 * (uij+uijp1)/2.0 - (vijm1+vip1jm1)/2.0 * (uijm1+uij)/2.0 ) + (alpha/_geom->Mesh()[1]) * ( std::abs(vij+vip1j)/2.0 * (uij-uijp1)/2.0 - std::abs(vijm1+vip1jm1)/2.0 * (uijm1-uij)/2.0 );
	return res;
}
*/

//============================================================================================
//discretization of the convective terms without Donor Cell Method (just for testing purporse)

real_t Grid::DC_duu_x(const Iterator &it, const real_t &alpha) const
{
	real_t res;
	real_t uij = _data[it.Value()];
	real_t uip1j = _data[it.Right().Value()];
	real_t uim1j = _data[it.Left().Value()];
	real_t uphalf = 0.5*(uip1j + uij);
	real_t umhalf = 0.5*(uij + uim1j);

	res = (1.0/_geom->Mesh()[0]) * (pow(uphalf,2.0) - pow(umhalf,2.0));
	return res;
}

real_t Grid::DC_dvv_y(const Iterator &it, const real_t &alpha) const
{
	real_t res;
	real_t vij = _data[it.Value()];
	real_t vijp1 = _data[it.Top().Value()];
	real_t vijm1 = _data[it.Down().Value()];
	real_t vphalf = 0.5*(vijp1 + vij);
	real_t vmhalf = 0.5*(vij + vijm1);

	res = (1.0/_geom->Mesh()[1]) * (pow(vphalf,2.0) - pow(vmhalf,2.0));
	return res;
}

real_t Grid::DC_duv_y(const Iterator &it, const real_t &alpha, const Grid* v) const
{
	real_t res;
	real_t viphj = 0.5*(v->Data()[it.Right().Value()] + v->Data()[it.Value()]);
	real_t viphjm1 = 0.5*(v->Data()[it.Right().Down().Value()] + v->Data()[it.Down().Value()]);
	real_t uijph = 0.5*(_data[it.Top().Value()] + _data[it.Value()]);
	real_t uijmh = 0.5*(_data[it.Value()] + _data[it.Down().Value()]);

	res = (1.0/_geom->Mesh()[1]) * (viphj * uijph - viphjm1 * uijmh);
	return res;
}

real_t Grid::DC_duv_x(const Iterator &it, const real_t &alpha, const Grid* u) const
{
	real_t res;
	real_t vip1jph = 0.5 * (_data[it.Right().Top().Value()] + _data[it.Right().Value()]);
	real_t vijph = 0.5 * (_data[it.Top().Value()] + _data[it.Value()]);
	real_t uiphj = 0.5 * (u->_data[it.Value()] + u->_data[it.Right().Value()]);
	real_t uimhj = 0.5 * (u->_data[it.Value()] + u->_data[it.Left().Value()]);

	res = (1.0/_geom->Mesh()[0]) * (vip1jph*uiphj - vijph*uiphj);
	return res;
}
//====================================================================================

real_t Grid::Max() const
{
	// TODO: test
	real_t res = _data[0];
	for(int i=0; i<_geom->Size()[0]*_geom->Size()[1]; i++)
	{
		if (_data[i]>res) res = _data[i];
	}
	return res;
}

real_t Grid::Min() const
{
	// TODO: test
	real_t res = _data[0];
	for(int i=0; i<_geom->Size()[0]*_geom->Size()[1]; i++)
	{
		if (_data[i]<res) res = _data[i];
	}
	return res;
}

real_t Grid::AbsMax() const
{
	// TODO: test
	real_t res = std::abs(_data[0]);
	real_t temp;
	for(int i=0; i<_geom->Size()[0]*_geom->Size()[1]; i++)
	{
		temp = std::abs(_data[i]);
		if (temp>res) res = temp;
	}
	return res;
}

real_t* Grid::Data()
{
	return _data;
}

// own method
const real_t* Grid::Data() const
{
	return _data;
}

Grid* Grid::copy() const
{
	Grid* res = new Grid(_geom, _offset);
	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		res->Cell(it) = Cell(it);
		it.Next();
	}
	return res;
}

void Grid::Out() const
{
	fprintf(stderr,"====================Output of Grid====================\n");
	Iterator it(_geom);
	it.First();
	while (it.Valid()) {
		fprintf(stderr,"%.5f   ",Cell(it));
		if(it.Right().Value() == it.Value())
			fprintf(stderr,"\n");
		it.Next();
	}
	fprintf(stderr,"====================End of Output!====================\n\n");
}

void Grid::Laplace(Grid* in)
{
	InteriorIterator it(_geom);
	it.First();
	while(it.Valid()) {
		_data[it.Value()] = in->dxx(it) + in->dyy(it);
		it.Next();
	}
}

bool Grid::CheckNaN() const
{
	bool nan(false);

	Iterator it(_geom);
	it.First();
	while (it.Valid()) {
		if(std::abs(_data[it.Value()]) > 1e5){
			nan = true;
			std::cout << "NaN found: " << it.Pos()[0] << ", " << it.Pos()[1] << "\n" << std::flush;
			std::cin.get(); // pause
		}
		it.Next();	
	}

	return nan;
}
