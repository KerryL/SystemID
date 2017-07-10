// File:  responseModeller.cpp
// Date:  7/10/2017
// Auth:  K. Loux
// Desc:  Model fitter class.

// Local headers
#include "responseModeller.h"

// Standard C++ headers
#include <cassert>

double ResponseModeller::ComputeModelError(const Eigen::VectorXd& parameters)
{
	if (modelType == ModelType::FirstOrder)
		ComputeModelledResponse(parameters(0));
	else// if (modelType == ModelType::SecondOrder)
		ComputeModelledResponse(parameters(0), parameters(1));

	maximumError = 0.0;
	double totalError(0.0);
	unsigned int i;
	for (i = 0; i < data.size(); ++i)
	{
		double error(fabs(data[i].response - modelledResponse[i]));
		if (rolloverPoint > 0.0)
			error = Unwind(error);

		totalError += error;
		if (error > maximumError)
			maximumError = error;
	}

	return totalError;
}

void ResponseModeller::ComputeModelledResponse(const double& bandwidthFrequency,
	const double& dampingRatio)
{
	modelledResponse.resize(data.size());
	double a, b1, b2;
	ComputeCoefficients(bandwidthFrequency, dampingRatio, a, b1, b2);

	modelledResponse[0] = data[0].response;
	modelledResponse[1] = data[1].response;

	unsigned int i;
	for (i = 2; i < data.size(); ++i)
		modelledResponse[i] = a * (data[i].input + 2 * data[i - 1].input + data[i - 2].input)
			- b1 * data[i - 1].response - b2 * data[i - 2].response;
}

void ResponseModeller::ComputeModelledResponse(const double& bandwidthFrequency)
{
	modelledResponse.resize(data.size());
	double a, b;
	ComputeCoefficients(bandwidthFrequency, a, b);

	modelledResponse[0] = data[0].response;

	unsigned int i;
	for (i = 1; i < data.size(); ++i)
		modelledResponse[i] = a * (data[i].input + data[i - 1].input) - b * data[i - 1].response;
}

void ResponseModeller::ComputeCoefficients(const double& bandwidthFrequency,
	const double& dampingRatio, double& a, double& b1, double& b2)
{
	const double denominator(4 + 4 * dampingRatio * bandwidthFrequency * sampleTime
		+ bandwidthFrequency * bandwidthFrequency * sampleTime * sampleTime);
	a = (bandwidthFrequency * bandwidthFrequency * sampleTime * sampleTime) / denominator;
	b1 = (2 * bandwidthFrequency * bandwidthFrequency * sampleTime * sampleTime - 8) / denominator;
	b2 = (4 - 4 * dampingRatio * bandwidthFrequency * sampleTime
		+ bandwidthFrequency * bandwidthFrequency * sampleTime * sampleTime) / denominator;
}

void ResponseModeller::ComputeCoefficients(const double& bandwidthFrequency,
	double& a, double& b)
{
	const double denominator(bandwidthFrequency * sampleTime + 2.0);
	a = bandwidthFrequency * sampleTime / denominator;
	b = (bandwidthFrequency * sampleTime - 2.0) / denominator;
}

double ResponseModeller::Unwind(const double& value) const
{
	assert(rolloverPoint > 0.0);
	assert(value >= 0.0);

	const unsigned int rolloverCount(value / rolloverPoint);
	return value - rolloverCount * rolloverPoint;
}
