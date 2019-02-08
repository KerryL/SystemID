// File:  idMath.h
// Date:  2/8/2019
// Auth:  K. Loux
// Desc:  Collection of methods related to mathematical operations.

#ifndef ID_MATH_H_
#define ID_MATH_H_

// Eigen headers
#include <Eigen/Eigen>

// Standard C++ headers
#include <limits>// For QNaN

/// Namespace containing commonly used mathematical and number processing
/// methods.
namespace IDMath
{
	/// Constant for evauluating whether or not a number is "very small."
	const double NearlyZero = 1.0e-12;
	//const double QNAN = std::numeric_limits<double>::quiet_NaN();// Not currently used

	/// Determines if the specified value is close enough to zero to be
	/// considered zero.
	///
	/// \param n   Value to consider.
	/// \param eps Threshold value.
	///
	/// \returns True if the absolute value of \p n is less than \p eps.
	bool IsZero(const double &n, const double &eps = NearlyZero);

	/// Determines if the norm of the specified vector is close enough to zero
	/// to be considered zero.
	///
	/// \param v   Vector to consider.
	/// \param eps Threshold value.
	///
	/// \returns True if the norm value of \p v is less than \p eps.
	bool IsZero(const Eigen::VectorXd &v, const double &eps = NearlyZero);

	/// Checks to see if the specified value is not a number.
	///
	/// \param value Value to consider.
	///
	/// \returns True if \p value is not a number.
	template<typename T>
	bool IsNaN(const T &value);

	/// Checks to see if the specified value is infinite.
	///
	/// \param value Value to consider.
	///
	/// \returns True if \p value is infinite.
	template<typename T>
	bool IsInf(const T &value);

	/// Checks to see if the specified value is a number and is not infinite.
	///
	/// \param value Value to consider.
	///
	/// \returns True if \p value is a finite number.
	template<typename T>
	bool IsValid(const T &value);

	/// Forces the value to lie between the specified limits.
	///
	/// \param value      Value to clamp.
	/// \param lowerLimit Maximum allowed output value.
	/// \param upperLimit Minimum allowed outptu value.
	///
	/// \returns The \p upperLimit if \p value is greater than \p upperLimit,
	///          the \p lowerLimit if \p value is less than \p lowerLimit, or
	///          the unmodified \p value otherwise.
	double Clamp(const double &value, const double &lowerLimit,
		const double &upperLimit);

	/// Converts the argument to lie within the range <&plusmn><&pi>.
	///
	/// \param angle Angle in radians.
	///
	/// \returns An equivalent angle within the range <&plusmn><&pi>.
	double RangeToPlusMinusPi(const double &angle);

	/// Converts the argument to lie within the range <&plusmn>180 deg.
	///
	/// \param angle Angle in degrees.
	///
	/// \returns An equivalent angle within the range <&plusmn>180 deg.
	double RangeToPlusMinus180(const double &angle);

	/// Determines the sign of the argument.
	///
	/// \param value The value to consider.
	///
	/// \returns -1.0 if the \p value is negative, +1.0 if the \p value is
	///          positive, or 0.0 if the \p value is zero.
	double Sign(const double &value);

	/// Returns the required precision to represent the specified \p value with
	/// the specified number of significant digits.
	///
	/// \param value             The value to assess.
	/// \param significantDigits The desired number of significant digits.
	/// \param dropTrailingZeros Indicates whether or not trailing zeros should
	///                          be included in the count.
	///
	/// \returns The number of digits to the right of the decimal point
	///          required to show the specified value with the specified number
	///          of significant digits.
	unsigned int GetPrecision(const double &value,
		const unsigned int &significantDigits = 2,
		const bool &dropTrailingZeros = true);

	/// Returns the required precision to differentiate between adjacent
	/// graduations on a number line.
	///
	/// \param minimum         The minimum value to be considered.
	/// \param majorResolution The delta between the minimum value and the next
	///                        graduation.
	/// \param isLogarithmic   Flag indicating whether or not the graduations
	///                        are logarithmically spaced.
	///
	/// \returns The required precision to differentiate between adjacent
	///          graduations.
	unsigned int GetPrecision(const double &minimum,
		const double &majorResolution, const bool &isLogarithmic = false);
}

//=============================================================================
// Namespace:		IDMath
// Function:		IsNaN
//
// Description:		Determines if the specified number is or is not a number.
//
// Input Arguments:
//		value	= const T& to check
//
// Output Arguments:
//		None
//
// Return Value:
//		bool, true if the argument is NOT a number
//
//=============================================================================
template<typename T>
bool IDMath::IsNaN(const T &value)
{
	return value != value;
}

//=============================================================================
// Namespace:		IDMath
// Function:		IsInf
//
// Description:		Determines if the specified number is infinite.
//
// Input Arguments:
//		value	= const T&
//
// Output Arguments:
//		None
//
// Return Value:
//		bool, true if the argument is ininite
//
//=============================================================================
template<typename T>
bool IDMath::IsInf(const T &value)
{
	return std::numeric_limits<T>::has_infinity &&
		value == std::numeric_limits<T>::infinity();
}

//=============================================================================
// Namespace:		IDMath
// Function:		IsValid
//
// Description:		Determines if the specified value is a valid number.
//
// Input Arguments:
//		value	= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		bool, true if the argument is valid
//
//=============================================================================
template<typename T>
bool IDMath::IsValid(const T &value)
{
	return !IsNaN<T>(value) && !IsInf<T>(value);
}

#endif// PLOT_MATH_H_
