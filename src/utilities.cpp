// File:  utilities.cpp
// Date:  8/4/2016
// Auth:  K. Loux
// Desc:  Tire data modeling utilities.

// Standard C++ headers
#include <limits>
#include <sstream>
#include <iomanip>

// Local headers
#include "utilities.h"

//==========================================================================
// Namespace:		Utilities
// Function:		ComputeCoefficientOfDetermination
//
// Description:		Computes the r^2 value for the given data.
//
// Input Arguments:
//		observedValues	= const std::vector<double>&
//		fitValues		= const std::vector<double>&
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//==========================================================================
double Utilities::ComputeCoefficientOfDetermination(
	const std::vector<double>& observedValues,
	const std::vector<double>& fitValues)
{
	assert(observedValues.size() > 1);
	assert(observedValues.size() == fitValues.size());

	double totalSumSquares(0.0);
	double residualSumSquares(0.0);
	const double mean(ComputeMean(observedValues));

	unsigned int i;
	for (i = 0; i < observedValues.size(); i++)
	{
		totalSumSquares += (observedValues[i] - mean) * (observedValues[i] - mean);
		residualSumSquares += (observedValues[i] - fitValues[i]) * (observedValues[i] - fitValues[i]);
	}

	return 1.0 - residualSumSquares / totalSumSquares;
}

//==========================================================================
// Namespace:		Utilities
// Function:		ComputeStanardErrorOfRegression
//
// Description:		Computes the standard error of the regression for the
//					given data.
//
// Input Arguments:
//		observedValues	= const std::vector<double>&
//		fitValues		= const std::vector<double>&
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//==========================================================================
double Utilities::ComputeStanardErrorOfRegression(
	const std::vector<double>& observedValues,
	const std::vector<double>& fitValues)
{
	assert(observedValues.size() > 0);
	assert(observedValues.size() == fitValues.size());

	double residualSumSquares(0.0);
	unsigned int i;
	for (i = 0; i < observedValues.size(); i++)
		residualSumSquares += (observedValues[i] - fitValues[i]) * (observedValues[i] - fitValues[i]);

	return std::sqrt(residualSumSquares / observedValues.size());
}
