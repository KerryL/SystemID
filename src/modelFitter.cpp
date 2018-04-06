// File:  modelFitter.cpp
// Date:  7/10/2017
// Auth:  K. Loux
// Desc:  Model fitter class.

// Local headers
#include "modelFitter.h"
#include "nelderMead.h"
#include "responseModeller.h"

bool ModelFitter::DetermineParameters(const std::vector<std::vector<Slice>>& data,
	double& bandwidthFrequency, double& sampleTime)
{
	sampleTime = ComputeSampleTime(data);

	ResponseModeller modeller(data, sampleTime, rolloverPoint, ResponseModeller::ModelType::FirstOrder);
	auto objectiveFunction = std::bind(&ResponseModeller::ComputeModelError, &modeller, std::placeholders::_1);
	NelderMead<1> optimization(objectiveFunction, iterationLimit);
	Eigen::VectorXd initialGuess(1);
	initialGuess(0) = 1.0;
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
	Eigen::Vector2d initialGuess(1.0, 1.0);
	optimization.SetInitialGuess(initialGuess);

	Eigen::VectorXd parameters(optimization.Optimize());
	bandwidthFrequency = parameters(0);
	dampingRatio = parameters(1);

	iterationCount = optimization.GetIterationCount();
	maximumError = modeller.GetMaximumError();

	return iterationCount < iterationLimit;
}

double ModelFitter::ComputeSampleTime(const std::vector<std::vector<Slice>>& data)
{
	std::vector<double> sampleTimes(data.size());
	unsigned int i;
	for (i = 0; i < data.size(); ++i)
		sampleTimes[i] = (data[i].back().time - data[i].front().time) / (data[i].size() - 1);

	return std::accumulate(sampleTimes.begin(), sampleTimes.end(), 0.0) / sampleTimes.size();
}
