// File:  modelFitter.cpp
// Date:  7/10/2017
// Auth:  K. Loux
// Desc:  Model fitter class.

// Local headers
#include "modelFitter.h"

// Standard C++ headers
#include <algorithm>

bool ModelFitter::DetermineParameters(const std::vector<std::vector<Slice>>& data,
	double& bandwidthFrequency, double& sampleTime)
{
	sampleTime = ComputeSampleTime(data);

	ResponseModeller modeller(data, sampleTime, rolloverPoint, ResponseModeller::ModelType::FirstOrder);
	auto objectiveFunction = std::bind(&ResponseModeller::ComputeModelError, &modeller, std::placeholders::_1);
	NelderMead<1> optimization(objectiveFunction, iterationLimit);
	Eigen::VectorXd initialGuess(1);
	initialGuess(0) = bwGuess;
	optimization.SetInitialGuess(initialGuess);

	Eigen::VectorXd parameters(optimization.Optimize());
	bandwidthFrequency = parameters(0);

	iterationCount = optimization.GetIterationCount();
	maximumError = modeller.GetMaximumError();

	return iterationCount < iterationLimit;
}

bool ModelFitter::DetermineParameters(const std::vector<std::vector<Slice>>& data,
	double& bandwidthFrequency, double& dampingRatio, double& sampleTime)
{
	sampleTime = ComputeSampleTime(data);

	ResponseModeller modeller(data, sampleTime, rolloverPoint, ResponseModeller::ModelType::SecondOrder);
	auto objectiveFunction = std::bind(&ResponseModeller::ComputeModelError, &modeller, std::placeholders::_1);
	NelderMead<2> optimization(objectiveFunction, iterationLimit);
	Eigen::Vector2d initialGuess(bwGuess, 1.0);
	optimization.SetInitialGuess(initialGuess);

	Eigen::VectorXd parameters(optimization.Optimize());
	bandwidthFrequency = parameters(0);
	dampingRatio = parameters(1);

	iterationCount = optimization.GetIterationCount();
	maximumError = modeller.GetMaximumError();

	return iterationCount < iterationLimit;
}

std::vector<double> ModelFitter::BuildTriangularCoefficients(const unsigned int& order)
{
	std::vector<double> coefficients(order);
	std::iota(coefficients.begin(), coefficients.end(), 1.0);
	std::reverse(coefficients.begin(), coefficients.end());

	for (unsigned int i = 2; i < order; ++i)// for each row in the triangle
	{
		for (unsigned int j = 1; j < i; ++j)// for each element in the middle of this row
			coefficients[order - i + j - 1] += coefficients[order - i + j];
	}

	return coefficients;
}

bool ModelFitter::DetermineParameters(const std::vector<std::vector<Slice>>& data,
	std::vector<double>& numerator, std::vector<double>& denominator, const unsigned int& order, double& sampleTime)
{
	assert(order > 0);
	sampleTime = ComputeSampleTime(data);

	ResponseModeller modeller(data, sampleTime, rolloverPoint, ResponseModeller::ModelType::NthOrder);
	auto objectiveFunction = std::bind(&ResponseModeller::ComputeModelError, &modeller, std::placeholders::_1);
	NelderMead<Eigen::Dynamic> optimization(objectiveFunction, iterationLimit);
	Eigen::VectorXd initialGuess(2 * order + 1);
	const double somethingSmall(0.1);
	initialGuess.head(order) = Eigen::VectorXd::Ones(order) * somethingSmall;
	initialGuess(order) = pow(bwGuess, order);
	auto coefficients(BuildTriangularCoefficients(order));
	for (unsigned int i = 0; i < order; ++i)
		initialGuess(order + i + 1) = coefficients[i] * pow(bwGuess, i + 1);
	optimization.SetInitialGuess(initialGuess);

	Eigen::VectorXd parameters(optimization.Optimize());
	unsigned int i;
	numerator.resize(order + 1);
	denominator.resize(order);
	for (i = 0; i < numerator.size(); ++i)
		numerator[i] = parameters(i);
	for (auto& c : denominator)
		c = parameters(i++);

	iterationCount = optimization.GetIterationCount();
	maximumError = modeller.GetMaximumError();

	return iterationCount < iterationLimit;
}

bool ModelFitter::DetermineParameters(const std::vector<std::vector<Slice>>& data,
	const std::string& numerator, const std::string& denominator,
	std::map<std::string, double>& userParameters, double& sampleTime)
{
	assert(!numerator.empty());
	assert(!denominator.empty());
	assert(!userParameters.empty());

	sampleTime = ComputeSampleTime(data);

	std::vector<std::string> userParamVector;
	for (const auto& p : userParameters)
		userParamVector.push_back(p.first);

	ResponseModeller modeller(data, sampleTime, rolloverPoint, numerator, denominator, userParamVector);
	auto objectiveFunction = std::bind(&ResponseModeller::ComputeModelError, &modeller, std::placeholders::_1);
	NelderMead<Eigen::Dynamic> optimization(objectiveFunction, iterationLimit);
	Eigen::VectorXd initialGuess(Eigen::VectorXd::Ones(userParameters.size()));
	optimization.SetInitialGuess(initialGuess);

	Eigen::VectorXd parameters(optimization.Optimize());
	auto it(userParameters.begin());
	for (int i = 0; i < parameters.size(); ++i)
		it++->second = parameters(i);

	iterationCount = optimization.GetIterationCount();
	maximumError = modeller.GetMaximumError();

	return iterationCount < iterationLimit;
}

double ModelFitter::ComputeSampleTime(const std::vector<std::vector<Slice>>& data)
{
	std::vector<double> sampleTimes(data.size());
	for (unsigned int i = 0; i < data.size(); ++i)
		sampleTimes[i] = (data[i].back().time - data[i].front().time) / (data[i].size() - 1);

	return std::accumulate(sampleTimes.begin(), sampleTimes.end(), 0.0) / sampleTimes.size();
}
