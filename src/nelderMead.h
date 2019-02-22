// File:  nelderMead.h
// Date:  8/4/2016
// Auth:  K. Loux
// Desc:  Implementation of Nelder-Mead (simplex) non-linear optimization.
//        For minimizing an objective function.

#ifndef NELDER_MEAD_H_
#define NELDER_MEAD_H_

// Standard C++ headers
#include <type_traits>
#include <cassert>

// Local headers
#include "optimizer.h"

template <int paramCount>
class NelderMead : public Optimizer
{
public:
	NelderMead(ObjectiveFunction objectiveFunction, const unsigned int& iterationLimit);

	virtual ~NelderMead() {}

	// Initialization parameters
	void SetSimplexInitializerDelta(const double& delta) { simplexDelta = delta; }
	void SetSimplexInitializerDeltaNearZero(const double& delta) { simplexDeltaNearZero = delta; }
	void SetSimplexInitializerNearlyZero(const double& zero) { assert(nearlyZero > 0.0); nearlyZero = zero; }

	// Convergence criteria
	void SetConvergedDeltaStep(const double& delta) { assert(delta > 0.0); convergedDeltaStep = delta; }
	void SetConvergedDeltaFunction(const double& delta) { assert(delta > 0.0); convergedDeltaFunction = delta; }

	// Parameters for reflection, expansion, contraction and shrinking
	void SetReflectionFactor(const double& factor) { assert(factor > 0.0); reflectionFactor = factor; }
	void SetExpansionFactor(const double& factor) { assert(factor > 1.0); expansionFactor = factor; }
	void SetContractionFactor(const double& factor) { assert(factor > 0.0 && factor < 1.0); contractionFactor = factor; }
	void SetShrinkFactor(const double& factor) { assert(factor > 0.0 && factor < 1.0); shrinkFactor = factor; }

	// Parameters for assisting with stabilization of calculations that can result in NaN
	void SetAdjustmentFactor(const double& factor) { assert(factor > 0.0 && factor < 1.0); adjustmentFactor = factor; }
	void SetEvaluationIterationLimit(const unsigned int& limit) { evalIterationLimit = limit; }

	typedef Eigen::Matrix<double, paramCount, 1> PointVec;

	void SetInitialGuess(const PointVec& initialGuess) { guess = initialGuess; }

	Eigen::VectorXd Optimize() const override;

	inline unsigned int GetIterationCount() const { return iterationCount; }

private:
	double simplexDelta;
	double simplexDeltaNearZero;
	double nearlyZero;

	double convergedDeltaStep;
	double convergedDeltaFunction;

	double reflectionFactor;
	double expansionFactor;
	double contractionFactor;
	double shrinkFactor;

	double adjustmentFactor;
	unsigned int evalIterationLimit;

	mutable unsigned int iterationCount;

	// These typedefs are funky because we want to be able to use dynamically sized
	// matrices if requested
	typedef typename std::conditional<paramCount == Eigen::Dynamic,
		Eigen::MatrixXd,
		Eigen::Matrix<double, paramCount, paramCount + 1>>::type SimplexMat;
	typedef typename std::conditional<paramCount == Eigen::Dynamic,
		Eigen::VectorXd,
		Eigen::Matrix<double, paramCount + 1, 1>>::type ValueVec;
	typedef typename std::conditional<paramCount == Eigen::Dynamic,
		Eigen::MatrixXd,
		Eigen::Matrix<double, paramCount, paramCount>>::type NbyNMat;
	typedef typename std::conditional<paramCount == Eigen::Dynamic,
		Eigen::PermutationMatrix<paramCount, paramCount>,
		Eigen::PermutationMatrix<paramCount + 1, paramCount + 1>>::type PermutationMat;

	PointVec guess;
	void Initialize(SimplexMat& simplex, ValueVec& functionValue) const;
	bool IsConverged(const SimplexMat& simplex, const ValueVec& functionValue) const;

	void Shrink(SimplexMat& simplex, ValueVec& functionValue) const;

	static void SortByFunctionValue(SimplexMat& simplex, ValueVec& functionValue);
	static PointVec AverageColumns(const NbyNMat& m);

	struct DynamicResize
	{
		static void Resize(PointVec& v, const int& trueParamCount);
		static void ResizeValue(ValueVec& v, const int& trueParamCount);
		static void Resize(SimplexMat& m, const int& trueParamCount);
	};

	struct FixedResize
	{
		static void Resize(PointVec& /*v*/, const int& /*trueParamCount*/) {}
		static void ResizeValue(ValueVec& /*v*/, const int& /*trueParamCount*/) {}
		static void Resize(SimplexMat& /*m*/, const int& /*trueParamCount*/) {}
	};

	typedef typename std::conditional<paramCount == Eigen::Dynamic,
		DynamicResize, FixedResize>::type ConditionalResize;

	bool StableEvaluate(const PointVec& centroid, const PointVec& worst,
		double factor, PointVec& point, double& value) const;
};

// Standard C++ headers
#include <numeric>

// Local headers
#include "utilities.h"

//==========================================================================
// Class:			NelderMead
// Function:		NelderMead
//
// Description:		Constructor for NelderMead class.
//
// Input Arguments:
//		objectiveFunction	= ObjectiveFunction
//		iterationLimit		= const unsigned int&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//==========================================================================
template <int paramCount>
NelderMead<paramCount>::NelderMead(ObjectiveFunction objectiveFunction,
	const unsigned int& iterationLimit)
	: Optimizer(objectiveFunction, iterationLimit)
{
	simplexDelta = 0.05;
	simplexDeltaNearZero = 0.00025;
	nearlyZero = 1.0e-8;

	convergedDeltaStep = 1.0e-6;
	convergedDeltaFunction = 1.0e-6;

	reflectionFactor = 1.0;
	expansionFactor = 2.0;
	contractionFactor = 0.5;
	shrinkFactor = 0.5;

	adjustmentFactor = 0.8;
	evalIterationLimit = 100;
}

//==========================================================================
// Class:			NelderMead
// Function:		Optimize
//
// Description:		Sovles for and returns arguments for the objective
//					function that optimize the result.
//
// Input Arguments:
//		None
//
// Output Arguments:
//		None
//
// Return Value:
//		Eigen::VectorXd
//
//==========================================================================
template <int paramCount>
Eigen::VectorXd NelderMead<paramCount>::Optimize() const
{
	assert(guess.size() > 0);

	ValueVec functionValue;
	SimplexMat simplex;
	Initialize(simplex, functionValue);

	const double expansionSacle(reflectionFactor * expansionFactor);
	const double outsideContractionScale(reflectionFactor * contractionFactor);

	PointVec centroid;
	PointVec reflection;
	PointVec expansion;
	PointVec contraction;

	ConditionalResize::Resize(centroid, guess.size());
	ConditionalResize::Resize(reflection, guess.size());
	ConditionalResize::Resize(expansion, guess.size());
	ConditionalResize::Resize(contraction, guess.size());

	double reflectionValue;
	double expansionValue;
	double contractionValue;

	unsigned int i(0);
	while (!IsConverged(simplex, functionValue) && i < iterationLimit)
	{
		// Compute reflection of the worst point and its corresponding function value
		centroid = AverageColumns(simplex.leftCols(simplex.rows()));

		if (!StableEvaluate(centroid, simplex.template rightCols<1>(), reflectionFactor, reflection, reflectionValue))
			return centroid;// TODO:  Warn user?

		if (reflectionValue < functionValue(0))// Reflected point is a new minimum
		{
			// Compute the expansion of the worst point and its corresponding function value
			if (!StableEvaluate(centroid, simplex.template rightCols<1>(), expansionSacle, expansion, expansionValue))
				return centroid;// TODO:  Warn user?

			if (expansionValue < reflectionValue)
			{
				simplex.template rightCols<1>() = expansion;
				functionValue.template tail<1>()(0) = expansionValue;
			}
			else
			{
				simplex.template rightCols<1>() = reflection;
				functionValue.template tail<1>()(0) = reflectionValue;
			}
		}
		else if (reflectionValue < functionValue.template tail<2>()(0))// Reflected point is better than the 2nd-to-worst point
		{
			simplex.template rightCols<1>() = reflection;
			functionValue.template tail<1>()(0) = reflectionValue;
		}
		else// Reflected point is no better than any other point we've tested
		{
			if (reflectionValue < functionValue.template tail<1>()(0))// Reflected point is not the worst
			{
				// Perform an outside contraction
				if (!StableEvaluate(centroid, simplex.template rightCols<1>(), outsideContractionScale, contraction, contractionValue))
					return centroid;// TODO:  Warn user?

				if (contractionValue <= reflectionValue)
				{
					simplex.template rightCols<1>() = contraction;
					functionValue.template tail<1>()(0) = contractionValue;
				}
				else
					Shrink(simplex, functionValue);
			}
			else
			{
				// Perform an inside contraction
				if (!StableEvaluate(centroid, simplex.template rightCols<1>(), -contractionFactor, contraction, contractionValue))
					return centroid;// TODO:  Warn user?

				if (contractionValue < functionValue.template tail<1>()(0))
				{
					simplex.template rightCols<1>() = contraction;
					functionValue.template tail<1>()(0) = contractionValue;
				}
				else
					Shrink(simplex, functionValue);
			}
		}

		SortByFunctionValue(simplex, functionValue);
		i++;
	}

	iterationCount = i;

	return simplex.template leftCols<1>();
}

//==========================================================================
// Class:			NelderMead
// Function:		StableEvaluate
//
// Description:		Evaluates the objective function for the calculated point,
//					looping to adjust the factor if necessary to arrive at a
//					finite value.
//
// Input Arguments:
//		centroid			= const PointVec&
//		worst				= const PointVec&
//		factor				= double
//
// Output Arguments:
//		point				= PointVec&
//		value				= double&
//
// Return Value:
//		bool
//
//==========================================================================
template <int paramCount>
bool NelderMead<paramCount>::StableEvaluate(const PointVec& centroid,
	const PointVec& worst, double factor, PointVec& point, double& value) const
{
	for (unsigned int i = 0; i < iterationLimit; ++i)
	{
		point = (1.0 + factor) * centroid - factor * worst;
		value = objectiveFunction(point);
		if (std::isfinite(value))
			return true;

		factor *= adjustmentFactor;
	}
	
	return false;
}

//==========================================================================
// Class:			NelderMead
// Function:		Initialize
//
// Description:		Initializes the simplex and the function value vector.
//
// Input Arguments:
//		None
//
// Output Arguments:
//		simplex			= SimplexMat&
//		functionValue	= ValueVec&
//
// Return Value:
//		None
//
//==========================================================================
template <int paramCount>
void NelderMead<paramCount>::Initialize(SimplexMat& simplex, ValueVec& functionValue) const
{
	ConditionalResize::Resize(simplex, guess.size());
	ConditionalResize::ResizeValue(functionValue, guess.size());

	simplex.col(0) = guess;
	functionValue(0) = objectiveFunction(guess);

	int i;
	for (i = 1; i < simplex.cols(); i++)
	{
		simplex.col(i) = guess;
		if (fabs(simplex(i - 1, i)) < nearlyZero)
			simplex(i - 1, i) = simplexDeltaNearZero;
		else
			simplex(i - 1, i) *= (1.0 + simplexDelta);

		functionValue(i) = objectiveFunction(simplex.col(i));
	}

	SortByFunctionValue(simplex, functionValue);
}

//==========================================================================
// Class:			NelderMead
// Function:		SortByFunctionValue
//
// Description:		Sorts the simplex columns and functionValue according
//					to functionValue.
//
// Input Arguments:
//		simplex			= SimplexMat&
//		functionValue	= ValueVec&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//==========================================================================
template <int paramCount>
void NelderMead<paramCount>::SortByFunctionValue(SimplexMat& simplex, ValueVec& functionValue)
{
	// Determine the required order of the permutation matrix
	std::vector<int> order(functionValue.size());
	std::iota(order.begin(), order.end(), 0);
	std::vector<std::pair<double, int> > valueZip(Utilities::Zip(functionValue.data(), order));
	std::sort(valueZip.begin(), valueZip.end());
	Utilities::Unzip<double, int>(valueZip, nullptr, &order);

	PermutationMat p(functionValue.rows());
	int i;
	for (i = 0; i < p.size(); i++)
		p.indices()(i) = order[i];

	simplex *= p;
	functionValue.transpose() *= p;
}

//==========================================================================
// Class:			NelderMead
// Function:		IsConverged
//
// Description:		Checks convergence criteria.
//
// Input Arguments:
//		simplex			= const SimplexMat&
//		functionValue	= const ValueVec&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//==========================================================================
template <int paramCount>
bool NelderMead<paramCount>::IsConverged(const SimplexMat& simplex,
	const ValueVec& functionValue) const
{
	double maxDeltaFunction(0.0), maxDeltaStep(0.0);
	int i;
	for (i = 1; i < functionValue.size(); i++)
	{
		maxDeltaFunction = std::max(maxDeltaFunction,
			fabs(functionValue(i) - functionValue(0)));
		maxDeltaStep = std::max(maxDeltaStep,
			(simplex.col(i) - simplex.col(0)).norm());
	}

	return maxDeltaFunction < convergedDeltaFunction &&
		maxDeltaStep < convergedDeltaStep;
}

//==========================================================================
// Class:			NelderMead
// Function:		AverageColumns
//
// Description:		Averages the columns of the specified matrix to create an
//					average vector.
//
// Input Arguments:
//		m	= const NbyNMat&
//
// Output Arguments:
//		None
//
// Return Value:
//		Eigen::VectorXd
//
//==========================================================================
template <int paramCount>
typename NelderMead<paramCount>::PointVec NelderMead<paramCount>::AverageColumns(const NbyNMat& m)
{
	PointVec sum(m.col(0));
	int i;
	for (i = 1; i < m.cols(); i++)
		sum += m.col(i);

	return sum / m.cols();
}

//==========================================================================
// Class:			NelderMead
// Function:		Shrink
//
// Description:		Performs the "shrink" step of the Nelder-Mead algorithm.
//					(Sometimes referred to as "multiple contraction," i.e. all
//					points move toward the current minimum).
//
// Input Arguments:
//		simplex			= SimplexMat&
//		functionValue	= ValueVec&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//==========================================================================
template <int paramCount>
void NelderMead<paramCount>::Shrink(SimplexMat& simplex, ValueVec& functionValue) const
{
	int i;
	for (i = 1; i < simplex.cols(); i++)
	{
		simplex.col(i) = simplex.col(0) + shrinkFactor * (simplex.col(i) - simplex.col(0));
		functionValue(i) = objectiveFunction(simplex.col(i));
	}
}

//==========================================================================
// Class:			NelderMead::DynamicResize
// Function:		Resize
//
// Description:		Resizes the vector to the specified size.
//
// Input Arguments:
//		v				= PointVec&
//		trueParamCount	= const int&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//==========================================================================
template <int paramCount>
void NelderMead<paramCount>::DynamicResize::Resize(PointVec& v,
	const int& trueParamCount)
{
	v.resize(trueParamCount);
}

//==========================================================================
// Class:			NelderMead::DynamicResize
// Function:		ResizeValue
//
// Description:		Resizes the vector to the specified size.
//
// Input Arguments:
//		v				= ValueVec&
//		trueParamCount	= const int&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//==========================================================================
template <int paramCount>
void NelderMead<paramCount>::DynamicResize::ResizeValue(ValueVec& v,
	const int& trueParamCount)
{
	v.resize(trueParamCount + 1);
}

//==========================================================================
// Class:			NelderMead::DynamicResize
// Function:		Resize
//
// Description:		Resizes the matrix to the specified size.
//
// Input Arguments:
//		m				= SimplexMat&
//		trueParamCount	= const int&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//==========================================================================
template <int paramCount>
void NelderMead<paramCount>::DynamicResize::Resize(
	SimplexMat& m, const int& trueParamCount)
{
	m.resize(trueParamCount, trueParamCount + 1);
}

#endif// NELDER_MEAD_H_
