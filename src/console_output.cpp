#include "console_output.hpp"

#include <string>

ConsoleClock::ConsoleClock()
: _position_counter(0), _repr(NULL), _n_pos(0), _n_lines(0), _n_cols(0)
{
	//init_classic_clock();
	//init_drawn_through();
	init_classic_clock_broader();
}

ConsoleClock::~ConsoleClock()
{
	if (_repr != NULL){
		for (int i=0; i<_n_pos; i++) delete[] _repr[i];
		delete[] _repr;
	}
}

const std::string& ConsoleClock::repr(int row) const
{
	return _repr[_position_counter % _n_pos][row];
}

void ConsoleClock::nextPos()
{
	_position_counter++;
}

void ConsoleClock::init_classic_clock()
{
	_n_pos = 8;
	_n_lines = 3;
	_n_cols = 3;
	_repr = new std::string*[_n_pos];
	for (int i=0; i<_n_pos; i++) _repr[i] = new std::string[_n_lines];

	_repr[0][0] = " | ";
	_repr[0][1] = " o ";
	_repr[0][2] = "   ";

	_repr[1][0] = "  /";
	_repr[1][1] = " o ";
	_repr[1][2] = "   ";

	_repr[2][0] = "   ";
	_repr[2][1] = " o-";
	_repr[2][2] = "   ";

	_repr[3][0] = "   ";
	_repr[3][1] = " o ";
	_repr[3][2] = "  \\";

	_repr[4][0] = "   ";
	_repr[4][1] = " o ";
	_repr[4][2] = " | ";

	_repr[5][0] = "   ";
	_repr[5][1] = " o ";
	_repr[5][2] = "/  ";

	_repr[6][0] = "   ";
	_repr[6][1] = "-o ";
	_repr[6][2] = "   ";

	_repr[7][0] = "\\  ";
	_repr[7][1] = " o ";
	_repr[7][2] = "   ";
}

void ConsoleClock::init_classic_clock_broader()
{
	_n_pos = 8;
	_n_lines = 3;
	_n_cols = 5;
	_repr = new std::string*[_n_pos];
	for (int i=0; i<_n_pos; i++) _repr[i] = new std::string[_n_lines];

	_repr[0][0] = "  |  ";
	_repr[0][1] = "  o  ";
	_repr[0][2] = "     ";

	_repr[1][0] = "   / ";
	_repr[1][1] = "  o  ";
	_repr[1][2] = "     ";

	_repr[2][0] = "     ";
	_repr[2][1] = "  o--";
	_repr[2][2] = "     ";

	_repr[3][0] = "     ";
	_repr[3][1] = "  o  ";
	_repr[3][2] = "   \\ ";

	_repr[4][0] = "     ";
	_repr[4][1] = "  o  ";
	_repr[4][2] = "  |  ";

	_repr[5][0] = "     ";
	_repr[5][1] = "  o  ";
	_repr[5][2] = " /   ";

	_repr[6][0] = "     ";
	_repr[6][1] = "--o  ";
	_repr[6][2] = "     ";

	_repr[7][0] = " \\   ";
	_repr[7][1] = "  o  ";
	_repr[7][2] = "     ";
}

void ConsoleClock::init_drawn_through()
{
	_n_pos = 4;
	_n_lines = 3;
	_n_cols = 3;
	_repr = new std::string*[_n_pos];
	for (int i=0; i<_n_pos; i++) _repr[i] = new std::string[_n_lines];

	_repr[0][0] = " | ";
	_repr[0][1] = " | ";
	_repr[0][2] = " | ";

	_repr[1][0] = "  /";
	_repr[1][1] = " / ";
	_repr[1][2] = "/  ";

	_repr[2][0] = "   ";
	_repr[2][1] = "---";
	_repr[2][2] = "   ";

	_repr[3][0] = "\\  ";
	_repr[3][1] = " \\ ";
	_repr[3][2] = "  \\";
}
