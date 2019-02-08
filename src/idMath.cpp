// File:  idMath.cpp
// Date:  2/8/2019
// Auth:  K. Loux
// Desc:  Collection of methods related to mathematical operations.

// Local headers
#include "idMath.h"

// Standard C++ headers
#include <cstdlib>
#include <cassert>
#include <limits>
#include <cstdarg>
#include <sstream>

//=============================================================================
// Namespace:		IDMath
// Function:		IsZero
//
// Description:		Returns true if a number is small enough to regard as zero.
//
// Input Arguments:
//		n	= const double& to be checked for being close to zero
//
// Output Arguments:
//		None
//
// Return Value:
//		bool, true if the number is less than NEARLY_ZERO
//
//=============================================================================
bool IDMath::IsZero(const double &n, const double &eps)
{
	if (fabs(n) < eps)
		return true;
	else
		return false;
}

//=============================================================================
// Namespace:		IDMath
// Function:		IsZero
//
// Description:		Returns true if a number is small enough to regard as zero.
//					This function checks the magnitude of the Vector.
//
// Input Arguments:
//		v	= const Eigen::VectorXd& to be checked for being close to zero
//
// Output Arguments:
//		None
//
// Return Value:
//		bool, true if the magnitude is less than NEARLY_ZERO
//
//=============================================================================
bool IDMath::IsZero(const Eigen::VectorXd &v, const double &eps)
{
	if (v.norm() < eps)
		return true;

	return false;
}

//=============================================================================
// Namespace:		IDMath
// Function:		Clamp
//
// Description:		Ensures the specified value is between the limits.  In the
//					event that the value is out of the specified bounds, the
//					value that is returned is equal to the limit that the value
//					has exceeded.
//
// Input Arguments:
//		value		= const double& reference to the value which we want to
//					  clamp
//		lowerLimit	= const double& lower bound of allowable values
//		upperLimit	= const double& upper bound of allowable values
//
// Output Arguments:
//		None
//
// Return Value:
//		double, equal to the clamped value
//
//=============================================================================
double IDMath::Clamp(const double &value, const double &lowerLimit,
	const double &upperLimit)
{
	// Make sure the arguments are valid
	assert(lowerLimit < upperLimit);

	if (value < lowerLimit)
		return lowerLimit;
	else if (value > upperLimit)
		return upperLimit;

	return value;
}

//=============================================================================
// Namespace:		IDMath
// Function:		RangeToPlusMinusPi
//
// Description:		Adds or subtracts 2 * pi to the specified angle until the
//					angle is between -pi and pi.
//
// Input Arguments:
//		angle		= const double& reference to the angle we want to bound
//
// Output Arguments:
//		None
//
// Return Value:
//		double, equal to the re-ranged angle
//
//=============================================================================
double IDMath::RangeToPlusMinusPi(const double &angle)
{
	return fmod(angle + M_PI, 2.0 * M_PI) - M_PI;
}

//=============================================================================
// Namespace:		IDMath
// Function:		RangeToPlusMinus180
//
// Description:		Adds or subtracts 180 to the specified angle until the
//					angle is between -180 and 180.
//
// Input Arguments:
//		angle		= const double& reference to the angle we want to bound
//
// Output Arguments:
//		None
//
// Return Value:
//		double, equal to the re-ranged angle
//
//=============================================================================
double IDMath::RangeToPlusMinus180(const double &angle)
{
	return fmod(angle + 180.0, 360.0) - 180.0;
}

//=============================================================================
// Namespace:		IDMath
// Function:		Sign
//
// Description:		Returns 1.0 for positive, -1.0 for negative and 0.0 for
//					zero.
//
// Input Arguments:
//		value		= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//=============================================================================
double IDMath::Sign(const double &value)
{
	if (value > 0.0)
		return 1.0;
	else if (value < 0.0)
		return -1.0;
	else
		return 0.0;
}

//=============================================================================
// Namespace:		IDMath
// Function:		GetPrecision
//
// Description:		Determines the best number of digits after the decimal place
//					for a string representation of the specified value (for
//					use with printf-style %0.*f formatting.
//
// Input Arguments:
//		value				= const double&
//		significantDigits	= const unsigned int&
//		dropTrailingZeros	= const bool&
//
// Output Arguments:
//		None
//
// Return Value:
//		unsigned int
//
//=============================================================================
unsigned int IDMath::GetPrecision(const double &value,
	const unsigned int &significantDigits, const bool &dropTrailingZeros)
{
	int precision(significantDigits - static_cast<unsigned int>(floor(log10(value)) - 1));
	if (precision < 0)
		precision = 0;
	if (!dropTrailingZeros)
		return precision;

	std::ostringstream ss;
	ss.precision(precision);
	ss << value;

	std::string number(ss.str());
	unsigned int i;
	for (i = number.size() - 1; i > 0; --i)
	{
		if (number[i] == '0')
			--precision;
		else
			break;
	}

	if (precision < 0)
		precision = 0;

	return precision;
}

//=============================================================================
// Namespace:		IDMath
// Function:		GetPrecision
//
// Description:		Returns the required precision (digits past zero) to
//					distinguish between adjacent graduations.
//
// Input Arguments:
//		minimum			= const double&
//		majorResolution	= const double&
//		isLogarithmic	= const bool&
//
// Output Arguments:
//		None
//
// Return Value:
//		unsigned int
//
//=============================================================================
unsigned int IDMath::GetPrecision(const double &minimum,
	const double &majorResolution, const bool &isLogarithmic)
{
	double baseValue;
	if (isLogarithmic)
		baseValue = minimum;
	else
		baseValue = majorResolution;

	if (log10(baseValue) >= 0.0)
		return 0;

	return static_cast<unsigned int>(-log10(baseValue) + 1);
}
