#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include "communicator.hpp"

#include <stdio.h>      /* NULL */
#include <stdlib.h>     /* malloc, free */
#include <cmath>        // std::abs, pow
#include <iostream>

/* public methods */

/* constructor */
/**
Constructs a grid with offset. Note the convention for the offset used here: The offset has to be between -1.0 and 0.0 and gives the relative distance, i.e. offset=-1.0 means that the physical offset length is equal to the length of one grid-cell.
\param[in] geom the data on the grid size etc.
\param[in] offset offset to consider the fact, that not all data point are at the same physical points (needs to be between -1.0 and 0.0!!)
\param[in] verbose true if debugging information should be displayed
*/
Grid::Grid(const Geometry* geom, const multi_real_t& offset, bool verbose)
: _data(NULL), _offset(offset), _geom(geom)
{	
	if (verbose) std::cout << "Allocating memory with size " << _geom->Size()[0] << " * " << _geom->Size()[1] << "... " << std::flush;
	_data = new real_t[_geom->Size()[0] * _geom->Size()[1]];
	if (verbose) std::cout << "Done.\n" << std::flush;
}

/**
Constrcts a basic grid without offset
\param[in] geom the data on the grid size etc.
*/
Grid::Grid(const Geometry* geom)
: Grid(geom, multi_real_t(0.0,0.0))
{
}

/* destructor */
Grid::~Grid()
{
	delete[] _data;
}

/**
\param[in] value value to be written to all grid points
\param[in] verbose debugging output (default: false)
*/
void Grid::Initialize(const real_t& value, bool verbose)
{
	if (verbose) std::cout << "Initializing: " << _geom->Size()[0]*_geom->Size()[1] << "\n" << std::flush;
	for(index_t i=0; i<_geom->Size()[0]*_geom->Size()[1]; i++)
	{
		_data[i] = value;
	}
	if (verbose) std::cout << "Done!" << std::flush;
}

/**
\param[in] it position at which the grid should be evaluated
*/
real_t& Grid::Cell(const Iterator& it)
{
	return _data[it.Value()];
}

/**
\param[in] it position at which the grid should be evaluated
*/
const real_t& Grid::Cell(const Iterator& it) const
{
	return _data[it.Value()];
}

/**
The methods first calculates the four nearest indices on the grid (with offset!) and then does a double-linear interpolation to calculate the interpolated value at the given phyical coordinates
\param[in] pos physical position
\return value at physical position
*/
real_t Grid::Interpolate(const multi_real_t& pos) const
{
	real_t ix = ( (_geom->Size()[0] - 2.0)*(pos[0]/_geom->Length()[0]) ) - _offset[0];
	real_t iy = ( (_geom->Size()[1] - 2.0)*(pos[1]/_geom->Length()[1]) ) - _offset[1];
	
	if (_offset[0] > 0 || _offset[1] > 0)
		std::cout << "Warning: Positive Offset \n";

	index_t x1 = floor(ix);
	index_t x2 = ceil(ix);
	index_t y1 = floor(iy);
	index_t y2 = ceil(iy);

	index_t x = round(ix);
	index_t y = round(iy);

	//calculate values
	index_t vallu = x1 + y1*_geom->Size()[0];
	index_t valru = x2 + y1*_geom->Size()[0];
	index_t vallo = x1 + y2*_geom->Size()[0];	
	index_t valro = x2 + y2*_geom->Size()[0];

	index_t val = x + y*_geom->Size()[0];

	//return _data[val]; //TODO Remove
	if (_geom->isObstacle(Iterator(_geom, val))) return 0.0;

	/*if(vallu < 0){
		std::cout << "Warning, negative index value in interpolation: " << vallu << "\n" << std::flush;
	} else */
	if(valro >= _geom->Size()[0]*_geom->Size()[1]){
		std::cout << "Warnng, too large index value in interpolation: " << ix << ", " << iy << "\n" << std::flush;
	}

	real_t alpha = ix - x1;
	real_t beta = iy - y1;

	//interpolation
	real_t xinter1 = (1.0 - alpha) * _data[vallu] + alpha * _data[valru];
	real_t xinter2 = (1.0 - alpha) * _data[vallo] + alpha * _data[valro];

	return (1.0 - beta) * xinter1 + beta * xinter2;
}

/**
\param[in] it position where the difference quotient should be evaluated
\return left-sided x-difference quotient
*/
real_t Grid::dx_l(const Iterator& it) const
{
	return (_data[it.Value()] - _data[it.Left().Value()])/((_geom->Mesh())[0]);
}

/**
\param[in] it position where the difference quotient should be evaluated
\return right-sided x-difference quotient
*/
real_t Grid::dx_r(const Iterator& it) const
{
	return (_data[it.Right().Value()] - _data[it.Value()])/((_geom->Mesh())[0]);
}

/**
\param[in] it position where the difference quotient should be evaluated
\return left-sided y-difference quotient
*/
real_t Grid::dy_l(const Iterator& it) const
{
	return (_data[it.Value()] - _data[it.Down().Value()])/((_geom->Mesh())[1]);
}

/**
\param[in] it position where the difference quotient should be evaluated
\return right-sided y-difference quotient
*/
real_t Grid::dy_r(const Iterator& it) const
{
	return (_data[it.Top().Value()] - _data[it.Value()])/((_geom->Mesh())[1]);
}

// own methods, central difference quotients
/**
\param[in] it position where the difference quotient should be evaluated
\return central x-difference quotient
*/
real_t Grid::dx_c(const Iterator& it) const
{
	return 0.5 * (dx_r(it)+dx_l(it));
}

/**
\param[in] it position where the difference quotient should be evaluated
\return central y-difference quotient
*/
real_t Grid::dy_c(const Iterator& it) const
{
	return 0.5 * (dy_r(it)+dy_l(it));
}

/**
\param[in] it position where the difference quotient should be evaluated
\return central x-difference quotient of 2nd order
*/
real_t Grid::dxx(const Iterator &it) const
{
	//return (_data[it.Right().Value()] - 2.0 * _data[it.Value()] + _data[it.Left().Value()])/( (_geom->Mesh())[0] * (_geom->Mesh())[0] );
	return (_data[it.Right().Value()] - 2.0 * _data[it.Value()] + _data[it.Left().Value()])/pow(_geom->Mesh()[0],2.0);
}

/**
\param[in] it position where the difference quotient should be evaluated
\return central y-difference quotient of 2nd order
*/
real_t Grid::dyy(const Iterator& it) const
{
	//return (_data[it.Top().Value()] - 2.0 * _data[it.Value()] + _data[it.Down().Value()])/( (_geom->Mesh())[1] * (_geom->Mesh())[1] );
	return (_data[it.Top().Value()] - 2.0 * _data[it.Value()] + _data[it.Down().Value()])/pow(_geom->Mesh()[1],2.0);
}

/**
\param[in] it position where the approximation quotient should be evaluated
\param[in] alpha control parameter
\return difference quotient approximation to u*du/dx
*/
real_t Grid::DC_udu_x(const Iterator& it, const real_t& alpha) const
{
	// we use here that u*du/dx = 0.5 * d(u^2)/dx
	real_t res;
	real_t uij = _data[it.Value()];
	real_t uip1j = _data[it.Right().Value()];
	real_t uim1j = _data[it.Left().Value()];
	res = (1.0/(_geom->Mesh()[0])) * ( pow((uij + uip1j)/2.0, 2.0) - pow((uim1j + uij)/2.0, 2.0)) + (alpha/(_geom->Mesh()[0])) * ( std::abs(uij + uip1j)/2.0 * (uij - uip1j)/2.0 - std::abs(uim1j + uij)/2.0 * (uim1j - uij)/2.0 );
	res *= 0.5;
	return res;
}

/*real_t Grid::DC_vdu_y(const Iterator& it, const real_t& alpha, const Grid* v) const
{
	// TODO
}

real_t Grid::DC_udv_x(const Iterator& it, const real_t& alpha, const Grid* u) const
{
	// TODO
}*/

/**
\param[in] it position where the approximation quotient should be evaluated
\param[in] alpha control parameter
\return difference quotient approximation to v*dv/dy
*/
real_t Grid::DC_vdv_y(const Iterator& it, const real_t& alpha) const
{
	// we use here that v*dv/dx = 0.5 * d(v^2)/dx
	real_t res;
	real_t vij = _data[it.Value()];
	real_t vijp1 = _data[it.Top().Value()];
	real_t vijm1 = _data[it.Down().Value()];
	res = (1.0/(_geom->Mesh()[1])) * ( pow((vij + vijp1)/2.0, 2.0) - pow((vijm1 + vij)/2.0, 2.0)) + (alpha/(_geom->Mesh()[1])) * ( std::abs(vij + vijp1)/2.0 * (vij - vijp1)/2.0 - std::abs(vijm1 + vij)/2.0 * (vijm1 - vij)/2.0 );
	res *= 0.5;
	return res;
}

// The use of the functions  vdu_y and udv_x are unclear, since only d(uv)/dx and d(uv)/dy have been presented in the lecture. The last mentioned functions are implemented below.

// the original donor cell methods
/**
\param[in] it position where the approximation quotient should be evaluated
\param[in] alpha control parameter
\return difference quotient approximation to d(u^2)/dx
*/
real_t Grid::DC_duu_x(const Iterator &it, const real_t &alpha) const
{
	return 2.0 * DC_udu_x(it,alpha);
}

/**
\param[in] it position where the approximation quotient should be evaluated
\param[in] alpha control parameter
\return difference quotient approximation to d(v^2)/dx
*/
real_t Grid::DC_dvv_y(const Iterator &it, const real_t &alpha) const
{
	return 2.0 * DC_vdv_y(it,alpha);
}

/**
\param[in] it position where the approximation quotient should be evaluated
\param[in] alpha control parameter
\return difference quotient approximation to d(uv)/dx
*/
real_t Grid::DC_duv_x(const Iterator &it, const real_t &alpha, const Grid* u) const
{
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

/**
\param[in] it position where the approximation quotient should be evaluated
\param[in] alpha control parameter
\return difference quotient approximation to d(uv)/dy
*/
real_t Grid::DC_duv_y(const Iterator &it, const real_t &alpha, const Grid* v) const
{
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


//============================================================================================
//discretization of the convective terms without Donor Cell Method (just for testing purporse)
/*
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
*/
//====================================================================================

/**
\return maximal value of the grid
*/
real_t Grid::Max() const
{
	real_t res = _data[0];
	for(index_t i=0; i<_geom->Size()[0]*_geom->Size()[1]; i++)
	{
		if (_data[i]>res) res = _data[i];
	}
	return res;
}

real_t Grid::InnerMax() const
{
	InteriorIterator it(_geom);
	it.First();
	real_t res = Cell(it);
	while (it.Valid()){
		if (Cell(it)>res) res = Cell(it);
		it.Next();
	}
	return res;
}

real_t Grid::InnerMin() const
{
	InteriorIterator it(_geom);
	it.First();
	real_t res = Cell(it);
	while (it.Valid()){
		if (Cell(it)<res) res = Cell(it);
		it.Next();
	}
	return res;
}

/**
\return minimal value of the grid
*/
real_t Grid::Min() const
{
	real_t res = _data[0];
	for(index_t i=0; i<_geom->Size()[0]*_geom->Size()[1]; i++)
	{
		if (_data[i]<res) res = _data[i];
	}
	return res;
}

/**
\return maximal value of the grid
*/
real_t Grid::AbsMax() const
{
	real_t res = std::abs(_data[0]);
	real_t temp;
	for(index_t i=0; i<_geom->Size()[0]*_geom->Size()[1]; i++)
	{
		temp = std::abs(_data[i]);
		if (temp>res) res = temp;
	}
	return res;
}

/**
\return the maximal value of all grids
*/
real_t Grid::TotalMax() const
{
	// this needs to be called by all processes!
	return _geom->getCommunicator()->gatherMax(Max());
}

/**
\return the maximal value of the interior cell points of all grids
*/
real_t Grid::TotalInnerMax() const
{
	// this needs to be called by all processes!
	return _geom->getCommunicator()->gatherMax(InnerMax());
}

/**
\return the minimal value of the interior cell points of all grids
*/
real_t Grid::TotalInnerMin() const
{
	// this needs to be called by all processes!
	return _geom->getCommunicator()->gatherMin(InnerMin());
}

/**
\return the minimal value of all grids
*/
real_t Grid::TotalMin() const
{
	// this needs to be called by all processes!
	return _geom->getCommunicator()->gatherMin(Min());
}

/**
\return the maximal absolute value of all grids
*/
real_t Grid::TotalAbsMax() const
{
	// this needs to be called by all processes!
	return _geom->getCommunicator()->gatherMax(AbsMax());
}

/**
Gives write acces directly to the data (better use Cell())
\return pointer to the raw data
*/
real_t* Grid::Data()
{
	return _data;
}

// own method
/**
Gives read access directly to the data (better use Cell())
\return constant pointer to the raw data
*/
const real_t* Grid::Data() const
{
	return _data;
}

/**
Note that the copy is a deep copy, i.e. a completely new grid is generated here
\return copied grid
*/
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

/**
Prints the grid to the terminal.
*/
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

/**
Calculates the laplacian of the data in the interior grid.
Note, that the data on the calling grid object are overwritten!!
\param[in] in Grid of which the laplacian method should be calculated
*/
void Grid::Laplace(Grid* in)
{
	InteriorIterator it(_geom);
	it.First();
	while(it.Valid()) {
		_data[it.Value()] = in->dxx(it) + in->dyy(it);
		it.Next();
	}
}

/**
Debugging method: Checks the interior grid for extremly large values
\return bool if bad entry found
*/
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

/**
\return mean value
*/
real_t Grid::average_value() const
{
	// compute the average value
	real_t avg(0.0);
	Iterator it(_geom);
	it.First();
	while (it.Valid()){
		avg += _data[it.Value()];
		it.Next();
	}
	return avg / real_t(_geom->Size()[0] * _geom->Size()[1]);
}

/**
\return the offset of this grid
*/
const multi_real_t& Grid::getOffset() const
{
	return _offset;
}

/**
\return a pointer to the underlying geometry object
*/
const Geometry* Grid::getGeometry() const
{
	return _geom;
}
