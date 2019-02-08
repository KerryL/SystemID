// File:  slice.h
// Date:  2/8/2019
// Auth:  K. Loux
// Desc:  Slice structure.

#ifndef SLICE_H_
#define SLICE_H_

struct Slice
{
	Slice(const double& time, const double& input, const double& response)
		: time(time), input(input), response(response) {}

	double time;// [sec]
	double input;
	double response;
};

#endif// SLICE_H_
