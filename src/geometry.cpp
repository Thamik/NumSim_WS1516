#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "communicator.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* atof */

/* public methods */

/* constructor */
/**
Constructs a Geometry object with standart values<br>
- 128x128 grid
- physical x-length: 1
- physical y-length: 1 
*/
Geometry::Geometry()
: _size(128,128), _length(1.0,1.0), _h(1.0,1.0), _velocity(0.0,0.0), _pressure(1.0) //standard values
{
	set_meshwidth(); // set _h to the right values
	// TODO: are the values right?
}

Geometry::Geometry(const Communicator *comm)
{
	_comm = comm;
	//TODO: something else to do?
}

/**
The geometry for the actual problem is given in the file in path file.
\param[in] file file path
\param[in] verbose printing debugging information
*/
void Geometry::Load(const char *file, bool verbose)
{
	if (verbose){
		std::cout << "Loading geometry file from path " << file << " ...\n";
	}
	std::string temp_string;
	std::ifstream infile;
	infile.open(file);
	if (!infile.is_open()){
		std::cout << "Warning: geometry file could not be read!\n";
		return;
	}
	for(int i=1;i<=7;i++)
	{
		getline(infile,temp_string);
		switch(i)
		{
			case 1:
				_size[0] = atoi(temp_string.c_str());
				break;
			case 2:
				_size[1] = atoi(temp_string.c_str());
				break;
			case 3:
				_length[0] = atof(temp_string.c_str());
				break;
			case 4:
				_length[1] = atof(temp_string.c_str());
				break;
			case 5:
				_velocity[0] = atof(temp_string.c_str());
				break;
			case 6:
				_velocity[1] = atof(temp_string.c_str());
				break;
			case 7:
				_pressure = atof(temp_string.c_str());
				break;
		}
	}
	infile.close();
	set_meshwidth(); // update _h
	//if (verbose){
		std::cout << "--------------------------------------------------\n";
		std::cout << "Geometry configuration read from file:\n";
		std::cout << "Size\t\t=\t(" << _size[0] << ", " << _size[1] << ")\n";
		std::cout << "Length\t\t=\t(" << _length[0] << ", " << _length[1] << ")\n";
		std::cout << "Velocity\t=\t(" << _velocity[0] << ", " << _velocity[1] << ")\n";
		std::cout << "Pressure\t=\t" << _pressure << "\n";
		std::cout << "--------------------------------------------------\n";
	//}
}

/**
\return number of cells in each direction
*/
const multi_index_t& Geometry::Size() const
{
	return _size;
}

const multi_index_t& Geometry::TotalSize() const
{
	//TODO
}

/**
\return physical length of the domain
*/
const multi_real_t& Geometry::Length() const
{
	return _length;
}

const multi_real_t &Geometry::TotalLength() const
{
	//TODO
}

/**
\return width of a single grid cell
*/
const multi_real_t& Geometry::Mesh() const
{
	return _h;
}

/**
\return pointer to the communicator instance
*/
const Communicator* Geometry::getCommunicator() const
{
	return _comm;
}


//==============================================================================
// Methods for the update of boundary values (have to be done every timestep)
//==============================================================================

// TODO: something to edit here?

/**
Updates the boundary values for Grid u according to the pattern of u, i.e. Dirichlet boundary conditions
\param[out] u Grid in which the boundary conditions should be set
*/
void Geometry::Update_U(Grid *u) const
{
	//std::cout << "Geometry: Update_U\n" << std::flush; // only for debugging issues
	// see lecture, 3.1.2
	for (int i=1; i<=4; i++){
		// check if this is a domain boundary
		if (!is_global_boundary(i)) continue;

		BoundaryIterator it(this);
		it.SetBoundary(i);
		it.First();
		while (it.Valid()){
			//std::cout << "i = " << i << ", it = " << it.Value() << "\n" << std::flush; // only for debugging issues
			if (i==4){
				// upper boundary
				u->Cell(it) = 2*1.0 - u->Cell(it.Down());
			} else if (i==1){
				// left boundary
				u->Cell(it) = 0.0;
			} else if (i==2) {
				// right boundary
				u->Cell(it) = 0.0;
				u->Cell(it.Left()) = 0.0;
			} else {
				// lower boundary
				u->Cell(it) = - u->Cell(it.Top());
			}
			it.Next();
		}
	}
}

/**
Updates the boundary values for Grid v according to the pattern of v, i.e. Dirichlet boundary conditions
\param[out] v Grid in which the boundary conditions should be set
*/
void Geometry::Update_V(Grid *v) const
{
	// see lecture, 3.1.2
	for (int i=1; i<=4; i++){
		// check if this is a domain boundary
		if (!is_global_boundary(i)) continue;

		BoundaryIterator it(this);
		it.SetBoundary(i);
		it.First();
		while (it.Valid()){
			if (i==4){
				// upper boundary
				v->Cell(it) = 0.0;
				v->Cell(it.Down()) = 0.0;
			} else if (i==3){
				// lower boundary
				v->Cell(it) = 0.0;
			} else if (i == 1) {
				v->Cell(it) = -v->Cell(it.Right());
			} else if (i == 2) {
				v->Cell(it) = -v->Cell(it.Left());
			}
			it.Next();
		}
	}
}

/**
Updates the boundary values for Grid p according to the pattern of p, i.e. Neumann boundary conditions
\param[out] p Grid in which the boundary conditions should be set
*/
void Geometry::Update_P(Grid *p) const
{
	// see lecture, 3.2.3
	for (int i=1; i<=4; i++){
		// check if this is a domain boundary
		if (!is_global_boundary(i)) continue;

		BoundaryIterator it(this);
		it.SetBoundary(i);
		it.First();
		while (it.Valid()){
			if (i==4){
				// upper boundary
				p->Cell(it) = p->Cell(it.Down());
			} else if (i==1){
				// left boundary
				p->Cell(it) = p->Cell(it.Right());
			} else if (i==2){
				// right boundary
				p->Cell(it) = p->Cell(it.Left());
			} else if (i==3){
				// lower boundary
				p->Cell(it) = p->Cell(it.Top());
			}
			it.Next();
		}
	}
}

// own method
/**
Calculates the meshwidth such that it satisfies h = length/size
*/
void Geometry::set_meshwidth()
{
	_h[0] = _length[0]/_size[0];
	_h[1] = _length[1]/_size[1];
}

/**
Decides whether the given boundary is a global boundary or is not
*/
bool Geometry::is_global_boundary(int boundary_index) const
{
	switch (boundary_index)
	{
		case BoundaryIterator::boundaryBottom :
			return _comm->isBottom();
			break;
		case BoundaryIterator::boundaryLeft :
			return _comm->isLeft();
			break;
		case BoundaryIterator::boundaryTop :
			return _comm->isTop();
			break;
		case BoundaryIterator::boundaryRight :
			return _comm->isRight();
			break;
	}
}
