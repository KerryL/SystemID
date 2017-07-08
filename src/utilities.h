// File:  utilities.h
// Date:  8/4/2016
// Auth:  K. Loux
// Desc:  Tire data modeling utilities.

#ifndef UTILITIES_H_
#define UTILITIES_H_

// Standard C++ headers
#include <vector>
#include <algorithm>
#include <utility>
#include <fstream>
#include <cassert>

namespace Utilities
{
	template <typename T>
	double ComputeStandardDeviation(const std::vector<T>& v);
	template <typename T>
	double ComputeMean(const std::vector<T>& v);
	template <typename T>
	T ComputeMedian(const std::vector<T>& v);

	template <typename Iter>
	double ComputeStandardDeviation(Iter start, Iter end);
	template <typename Iter>
	double ComputeMean(Iter start, Iter end);
	template <typename Iter>
	double ComputeMedian(Iter start, Iter end);

	double ComputeCoefficientOfDetermination(
		const std::vector<double>& observedValues,
		const std::vector<double>& fitValues);
	double ComputeStanardErrorOfRegression(
		const std::vector<double>& observedValues,
		const std::vector<double>& fitValues);

	typedef std::vector<double>::size_type RowType;

	template<typename T>
	static std::vector<T> ExtractSubVector(const std::vector<T>& v,
		const RowType& start, const RowType& end);
	template<typename T>
	static std::vector<T> ExtractSubVector(const std::vector<T>& v,
		const RowType& start);

	template <typename T>
	int Sign(const T& value)
	{
		return (T(0) < value) - (value < T(0));
	}

	template <typename T1, typename T2>
	std::vector<std::pair<T1, T2> > Zip(const std::vector<T1>& a, const std::vector<T2>& b)
	{
		assert(a.size() == b.size());

		std::vector<std::pair<T1, T2> > z(a.size());
		unsigned int i;
		for (i = 0; i < z.size(); i++)
		{
			z[i].first = a[i];
			z[i].second = b[i];
		}

		return z;
	}

	template <typename T1, typename T2>
	std::vector<std::pair<T1, T2> > Zip(const std::vector<T1>& a, const T2* b)
	{
		std::vector<std::pair<T1, T2> > z(a.size());
		unsigned int i;
		for (i = 0; i < z.size(); i++)
		{
			z[i].first = a[i];
			z[i].second = b[i];
		}

		return z;
	}

	template <typename T1, typename T2>
	std::vector<std::pair<T1, T2> > Zip(const T1* a, const std::vector<T2>& b)
	{
		std::vector<std::pair<T1, T2> > z(b.size());
		unsigned int i;
		for (i = 0; i < z.size(); i++)
		{
			z[i].first = a[i];
			z[i].second = b[i];
		}

		return z;
	}

	template <typename T1, typename T2>
	void Unzip(const std::vector<std::pair<T1, T2> >& z, std::vector<T1>* a, std::vector<T2>* b)
	{
		unsigned int i;

		if (a)
		{
			a->resize(z.size());
			for (i = 0; i < z.size(); i++)
				a->operator[](i) = z[i].first;
		}

		if (b)
		{
			b->resize(z.size());
			for (i = 0; i < z.size(); i++)
				b->operator[](i) = z[i].second;
		}
	}

	struct Sequence
	{
	public:
		Sequence(const double& minimum, const double& step)
			: minimum(minimum), step(step) { i = 0; }

		double operator()() { return minimum + ++i * step; }

	private:
		const double minimum;
		const double step;
		unsigned int i;
	};
}

// Standard C++ headers
#include <numeric>

//==========================================================================
// Namespace:		Utilities
// Function:		ComputeStandardDeviation
//
// Description:		Computes the standard deviation of the specified vector.
//
// Input Arguments:
//		v	= const std::vector<T>&
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//==========================================================================
template <typename T>
double Utilities::ComputeStandardDeviation(const std::vector<T>& v)
{
	return ComputeStandardDeviation(v.begin(), v.end());
}

//==========================================================================
// Namespace:		Utilities
// Function:		ComputeMean
//
// Description:		Computes the mean of the specified vector.
//
// Input Arguments:
//		v	= const std::vector<T>&
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//==========================================================================
template <typename T>
double Utilities::ComputeMean(const std::vector<T>& v)
{
	return ComputeMean(v.begin(), v.end());
}

//==========================================================================
// Namespace:		Utilities
// Function:		ComputeMedian
//
// Description:		Computes the median of the specified vector.
//
// Input Arguments:
//		v	= const std::vector<T>&
//
// Output Arguments:
//		None
//
// Return Value:
//		T
//
//==========================================================================
template <typename T>
T Utilities::ComputeMedian(const std::vector<T>& v)
{
	return ComputeMedian(v.cbegin(), v.cend());
}

//==========================================================================
// Namespace:		Utilities
// Function:		ComputeStandardDeviation
//
// Description:		Computes the standard deviation of the specified vector.
//
// Input Arguments:
//		start	= Iter
//		end		= Iter
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//==========================================================================
template <typename Iter>
double Utilities::ComputeStandardDeviation(Iter start, Iter end)
{
	const double mean(ComputeMean(start, end));
	std::vector<double> diff(std::distance(start, end));
	std::transform(start, end, diff.begin(), [mean](double x) { return x - mean; });
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	return std::sqrt(sq_sum / diff.size());
}

//==========================================================================
// Namespace:		Utilities
// Function:		ComputeMean
//
// Description:		Computes the mean of the specified vector.
//
// Input Arguments:
//		start	= Iter
//		end		= Iter
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//==========================================================================
template <typename Iter>
double Utilities::ComputeMean(Iter start, Iter end)
{
	return std::accumulate(start, end, static_cast<typename std::iterator_traits<Iter>::value_type>(0))
		/ static_cast<double>(std::distance(start, end));
}

//==========================================================================
// Namespace:		Utilities
// Function:		ComputeMedian
//
// Description:		Computes the median of the specified vector.
//
// Input Arguments:
//		start	= Iter
//		end		= Iter
//
// Output Arguments:
//		None
//
// Return Value:
//		double
//
//==========================================================================
template <typename Iter>
double Utilities::ComputeMedian(Iter start, Iter end)
{
	std::vector<typename std::iterator_traits<Iter>::value_type> vCopy(start, end);
	std::sort(vCopy.begin(), vCopy.end());
	return static_cast<double>(vCopy[vCopy.size() / 2]);
}

//==========================================================================
// Class:			Utilities
// Function:		ExtractSubVector
//
// Description:		Splits a portion out of the specified vector.
//
// Input Arguments:
//		v		= const std::vector<T>&
//		start	= const RowType&
//		end		= const RowType&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::vector<T>
//
//==========================================================================
template<typename T>
std::vector<T> Utilities::ExtractSubVector(const std::vector<T>& v,
	const RowType& start, const RowType& end)
{
	typename std::vector<T>::const_iterator first(v.begin() + start);
	typename std::vector<T>::const_iterator last(v.begin() + end);
	return std::vector<T>(first, last);
}

//==========================================================================
// Class:			Utilities
// Function:		ExtractSubVector
//
// Description:		Splits a portion out of the specified vector.
//
// Input Arguments:
//		v		= const std::vector<T>&
//		start	= const RowType&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::vector<T>
//
//==========================================================================
template<typename T>
std::vector<T> Utilities::ExtractSubVector(const std::vector<T>& v,
	const RowType& start)
{
	typename std::vector<T>::const_iterator first(v.begin() + start);
	return std::vector<T>(first, v.end());
}

#endif// UTILITIES_H_
