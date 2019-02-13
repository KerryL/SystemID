// File:  filter.cpp
// Date:  2/8/2019
// Auth:  K. Loux
// Desc:  Base class (abstract) for digital filters.

// Standard C++ headers
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <sstream>

// Local headers
#include "filter.h"
#include "expressionTree.h"
#include "idMath.h"

//=============================================================================
// Class:			Filter
// Function:		Filter
//
// Description:		Constructor for the Filter class.
//
// Input Arguments:
//		_sampleRate		= const double& specifying the sampling rate in Hz
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
Filter::Filter(const double &sampleRate) : sampleRate(sampleRate)
{
}

//=============================================================================
// Class:			Filter
// Function:		Filter
//
// Description:		Constructor for the Filter class for arbitrary filters.
//					Passed arguments are assumed to be for a continuous time
//					filter (s-domain), and will be translated into filter
//					coefficients according to the specified sample rate.
//
// Input Arguments:
//		sampleRate		= const double& specifying the sampling rate in Hz
//		numerator		= const std::vector<double> containing the numerator
//						  coefficients from highest power to zero power
//		denominaotr		= const std::vector<double> containing the denominator
//						  coefficients from highest power to zero power
//		initialValue	= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
Filter::Filter(const double &sampleRate, const std::vector<double> &numerator,
	const std::vector<double> &denominator, const double &initialValue)
	: sampleRate(sampleRate)
{
	GenerateCoefficients(numerator, denominator);
	Initialize(initialValue);
}

//=============================================================================
// Class:			Filter
// Function:		GenerateCoefficients
//
// Description:		Generates the discrete-time (z-domain) coefficients for
//					a filter equivalent to the continuous-time (s-domain)
//					arguments.  Uses bilinear transform:
//					s = 2 * (1 - z^-1) / (T * (1 + z^-1)).
//
// Input Arguments:
//		numerator	= const std::vector<double>& continuous time coefficients,
//					  highest power of s to lowest power of s
//		denominator	= const std::vector<double>& continuous time coefficients,
//					  highest power of s to lowest power of s
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
void Filter::GenerateCoefficients(const std::vector<double> &numerator,
	const std::vector<double> &denominator)
{
	const unsigned int highestPower(
		std::max(numerator.size(), denominator.size()) - 1);
	std::string numString = AssembleZExpression(numerator, highestPower);
	std::string denString = AssembleZExpression(denominator, highestPower);

	std::vector<double> zNum = CoefficientsFromString(numString);
	std::vector<double> zDen = CoefficientsFromString(denString);
	ResizeArrays(zNum.size(), zDen.size());

	unsigned int i;
	for (i = 0; i < zNum.size(); ++i)
		mA[i] = zNum[i] / zDen[0];
	for (i = 0; i < zDen.size() - 1; ++i)
		mB[i] = zDen[i + 1] / zDen[0];
}

//=============================================================================
// Class:			Filter
// Function:		AssembleZExpression
//
// Description:		Assembles the z-domain expression equivalent to the s-domain
//					coefficients provided.
//
// Input Arguments:
//		coefficients	= const std::vector<double>& continuous time coefficients,
//						  highest power of s to lowest power of s
//		highestPower	= const unsigned int&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string Filter::AssembleZExpression(
	const std::vector<double>& coefficients,
	const unsigned int &highestPower) const
{
	std::ostringstream ss;
	ss << "(" << std::fixed << 1.0 / sampleRate << "*(1+z^-1))";
	std::string posBilinTerm(ss.str()), negBilinTerm("(2*(1-z^-1))");
	std::string result;

	unsigned int i;
	for (i = 0; i < coefficients.size(); ++i)
	{
		if (IDMath::IsZero(coefficients[i]))
			continue;

		std::ostringstream term;
		term << std::fixed << coefficients[i];
		if (!result.empty() && coefficients[i] > 0.0)
			result.append("+");
		result.append(term.str());

		if (coefficients.size() - 1 > i)
			result.append("*" + negBilinTerm);
		if (coefficients.size() - 1 > 1 + i)
		{
			std::ostringstream power;
			power << "^" << coefficients.size() - i - 1;
			result.append(power.str());
		}
		if (highestPower + i > coefficients.size() - 1)
			result.append("*" + posBilinTerm);
		if (highestPower + i > coefficients.size())
		{
			std::ostringstream power;
			power << "^" << highestPower - coefficients.size() + 1 + i;
			result.append(power.str());
		}
	}

	if (result.empty())
		result = "0";

	return result;
}

//=============================================================================
// Class:			Filter
// Function:		Initialize
//
// Description:		Initializes the filter to the specified value.
//
// Input Arguments:
//		initialValue	= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
void Filter::Initialize(const double &initialValue)
{
	for (auto& v : mU)
		v = initialValue;

	double ssGain(ComputeSteadyStateGain());
	for (auto& v : mY)
		v = initialValue * ssGain;
}

//=============================================================================
// Class:			Filter
// Function:		Apply
//
// Description:		Applies the filter to the new input value.
//
// Input Arguments:
//		u0	= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		double containing the filtered value
//
//=============================================================================
double Filter::Apply(const double &u0)
{
	ShiftArray(mU);
	mU[0] = u0;

	ShiftArray(mY);

	mY[0] = 0.0;
	unsigned int i;
	for (i = 0; i < mA.size(); ++i)
		mY[0] += mA[i] * mU[i];
	for (i = 1; i < mY.size(); ++i)
		mY[0] -= mB[i - 1] * mY[i];

	return mY[0];
}

//=============================================================================
// Class:			Filter
// Function:		ShiftArray
//
// Description:		Shifts the array values by one index (value with highest
//					index is lost).
//
// Input Arguments:
//		s		= std::vector<double>&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
void Filter::ShiftArray(std::vector<double>& s) const
{
	unsigned int i;
	for (i = s.size() - 1; i > 0; --i)
		s[i] = s[i - 1];
}

//=============================================================================
// Class:			Filter
// Function:		ResizeArrays
//
// Description:		Allocates the coefficient and input/output storage arrays.
//
// Input Arguments:
//		_inSize		= const unsigned int&
//		_outSize	= const unsigned int&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
void Filter::ResizeArrays(const unsigned int &inSize,
	const unsigned int &outSize)
{
	mA.resize(inSize);
	mB.resize(outSize - 1);
	mU.resize(inSize);
	mY.resize(outSize);
}

//=============================================================================
// Class:			Filter
// Function:		CoefficientsFromString
//
// Description:		Creates a vector of coefficients from highest power to zero
//					power, based on a string representing the expression.
//
// Input Arguments:
//		s	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::vector<double>
//
//=============================================================================
std::vector<double> Filter::CoefficientsFromString(const std::string &s)
{
	ExpressionTree e;
	std::string expression;
	std::string errorString = e.Solve(s, expression);
	if (!errorString.empty())
	{
		// TODO:  Generate a warning here?
	}

	std::vector<std::pair<int, double>> terms =
		ExpressionTree::FindPowersAndCoefficients(
			ExpressionTree::BreakApartTerms(expression));

	terms = CollectLikeTerms(terms);
	terms = PadMissingTerms(terms);
	std::sort(terms.begin(), terms.end());

	std::vector<double> coefficients;
	unsigned int i;
	for (i = 0; i < terms.size(); ++i)
		coefficients.push_back(terms[terms.size() - 1 - i].second);

	return coefficients;
}

//=============================================================================
// Class:			Filter
// Function:		CollectLikeTerms
//
// Description:		Collects all terms with the same exponent and adds the
//					coefficients together.
//
// Input Arguments:
//		terms	= std::vector<std::pair<int, double>>
//
// Output Arguments:
//		None
//
// Return Value:
//		std::vector<std::pair<int, double>>
//
//=============================================================================
std::vector<std::pair<int, double>> Filter::CollectLikeTerms(
	std::vector<std::pair<int, double>> terms)
{
	unsigned int i, j;
	for (i = 0; i < terms.size(); ++i)
	{
		for (j = i + 1; j < terms.size(); ++j)
		{
			if (terms[i].first == terms[j].first)
			{
				terms[i].second += terms[j].second;
				terms.erase(terms.begin() + j);
				--j;
			}
		}
	}

	return terms;
}

//=============================================================================
// Class:			Filter
// Function:		PadMissingTerms
//
// Description:		If a power between the maximum power and zero is missing,
//					a zero-coefficient value for that power is inserted at
//					the appropriate location in the vector.
//
// Input Arguments:
//		terms	= std::vector<std::pair<int, double>>
//
// Output Arguments:
//		None
//
// Return Value:
//		std::vector<std::pair<int, double>>
//
//=============================================================================
std::vector<std::pair<int, double>> Filter::PadMissingTerms(
	std::vector<std::pair<int, double>> terms)
{
	std::sort(terms.begin(), terms.end(),
		std::greater<std::pair<int, double>>());// Sort in descending order of power

	int expectedPower(terms[0].first - 1);
	while (expectedPower < -1)
	{
		++expectedPower;
		terms.insert(terms.begin(), std::make_pair(expectedPower + 1, 0.0));
	}

	unsigned int i;
	for (i = 1; i < terms.size(); ++i)
	{
		if (terms[i].first != expectedPower)
			terms.insert(terms.begin() + i, std::make_pair(expectedPower, 0.0));
		--expectedPower;
	}

	while (expectedPower >= 0)
	{
		terms.insert(terms.end(), std::make_pair(expectedPower, 0.0));
		--expectedPower;
	}

	return terms;
}

//=============================================================================
// Class:			Filter
// Function:		ComputeSteadyStateGain
//
// Description:		Returns the steady-state value resluting from a unity step
//					input.  Returns zero if the results of this analysis are
//					not guaranteed valid - den must have these properties for
//					Final Value Theorem to yield valid results:
//					1.  Roots of den must have negative real parts
//					2.  den must have no more than one zero root
//					The final value is evaluated by multiplying by s and the
//					input (step input is 1/s, so we multiply by s/s, which
//					means we do nothing), canceling esses, then evaluating for
//					s=0.
//
// Input Arguments:
//		num	= const std::string&
//		den = const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//=============================================================================
double Filter::ComputeSteadyStateGain(const std::string &num,
	const std::string &den)
{
	std::vector<double> numeratorCoefficients = CoefficientsFromString(num);
	std::vector<double> denominatorCoefficients = CoefficientsFromString(den);

	unsigned int numEndZeros(0), denEndZeros(0);
	int i;
	for (i = numeratorCoefficients.size() - 1; i >= 0; --i)
	{
		if (IDMath::IsZero(numeratorCoefficients[i]))
			++numEndZeros;
		else
			break;
	}

	for (i = denominatorCoefficients.size() - 1; i >= 0; --i)
	{
		if (IDMath::IsZero(denominatorCoefficients[i]))
			++denEndZeros;
		else
			break;
	}

	// TODO:  Check condition 1

	// Check condition 2
	if (denEndZeros > 1)
		return 0.0;

	unsigned int essesToCancel(std::min<unsigned int>(numEndZeros, denEndZeros));

	// When evaluating for s=0, everything except the polynomial's constant term drops out
	return numeratorCoefficients[numeratorCoefficients.size() - 1 - essesToCancel] /
		denominatorCoefficients[denominatorCoefficients.size() - 1 - essesToCancel];
}

//=============================================================================
// Class:			Filter
// Function:		ComputeSteadyStateGain
//
// Description:		Returns the steady-state value resluting from a unity step
//					input.  Returns zero if the results of this analysis are
//					not guaranteed valid (i.e. if the steady-state gain is
//					undefined).
//
// Input Arguments:
//		None
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//=============================================================================
double Filter::ComputeSteadyStateGain() const
{
	double numeratorSum(0.0);
	double denominatorSum(1.0);
	for (const auto& v : mA)
		numeratorSum += v;
	for (const auto& v : mB)
		denominatorSum += v;

	return numeratorSum / denominatorSum;
}
