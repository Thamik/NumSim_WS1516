#include "typedef.hpp"

#include <string>

#ifndef __CONSOLE_OUTPUT_HPP
#define __CONSOLE_OUTPUT_HPP

class ConsoleClock {
public:
	ConsoleClock();
	~ConsoleClock();

	const std::string& repr(int row) const;

	void nextPos();

private:
	int _position_counter;

	std::string** _repr;
	int _n_pos;
	int _n_lines;
	int _n_cols;

	void init_classic_clock();
	void init_classic_clock_broader();
	void init_drawn_through();
};

#endif // __CONSOLE_OUTPUT_HPP
