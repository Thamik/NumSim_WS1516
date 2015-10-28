#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* atof */

/* public methods */

/* constructor */
Geometry::Geometry()
: _size(5,5), _length(1.0,1.0), _h(1.0,1.0), _velocity(0.0,0.0), _pressure(1.0)
//: _size(128,128), _length(1.0,1.0), _h(1.0,1.0), _velocity(0.0,0.0), _pressure(1.0) //standart values
{
	set_meshwidth(); // set _h to the right values
	// TODO: are the values right?
}

void Geometry::Load(const char *file)
{
	// TODO: Test this method
	std::string temp_string;
	std::ifstream infile;
	infile.open(file);
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

// update the boundary values (have to be done every timestep)
void Geometry::Update_U(Grid *u) const
{
	//std::cout << "Geometry: Update_U\n" << std::flush; // only for debugging issues
	// TODO: test
	// see lecture, 3.1.2
	for (int i=1; i<=4; i++){
		//std::cout << "!!!!!!!!!!!\n" << std::flush; // only for debugging issues
		BoundaryIterator it(this);
		it.SetBoundary(i);
		it.First();
		while (it.Valid()){
			//std::cout << "i = " << i << ", it = " << it.Value() << "\n" << std::flush; // only for debugging issues
			if (i==4){
				// upper boundary
				u->Cell(it) = 2*1.0 - u->Cell(it.Down());
			} else if (i==1 || i==2){
				// left or right boundary
				u->Cell(it) = 0.0;
			} else {
				// lower boundary
				u->Cell(it) = - u->Cell(it.Top());
			}
			it.Next();
		}
	}
}

// update the boundary values (have to be done every timestep)
void Geometry::Update_V(Grid *v) const
{
	// TODO: test
	// see lecture, 3.1.2
	for (int i=1; i<=4; i++){
		BoundaryIterator it(this);
		it.SetBoundary(i);
		it.First();
		while (it.Valid()){
			if (i==4 || i==3){
				// upper boundary
				//v->Cell(it) = - v->Cell(it.Down());
				v->Cell(it) = 0.0;
			/*} else if (i==1 || i==2){
				// left or right boundary
				v->Cell(it) = 0.0;*/
			} else if (i == 1) {
				v->Cell(it) = -v->Cell(it.Right());
			} else if (i == 2) {
				v->Cell(it) = -v->Cell(it.Left());
			} /*else {
				// lower boundary
				v->Cell(it) = - v->Cell(it.Top());
			}*/
			it.Next();
		}
	}
}

// update the boundary values for p
void Geometry::Update_P(Grid *p) const
{
	// TODO: test
	// see lecture, 3.2.3
	for (int i=1; i<=4; i++){
		BoundaryIterator it(this);
		it.SetBoundary(i);
		it.First();
		while (it.Valid()){
			if (i==4){
				// upper boundary
				//std::cout << it.Pos()[0] << ", " << it.Pos()[1] << "\n";
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
void Geometry::set_meshwidth()
{
	_h[0] = _length[0]/_size[0];
	_h[1] = _length[1]/_size[1];
}
