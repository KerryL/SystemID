// File:  responseModeller.cpp
// Date:  7/10/2017
// Auth:  K. Loux
// Desc:  Model fitter class.

// Local headers
#include "responseModeller.h"

// Standard C++ headers
#include <cassert>

void ResponseModeller::ComputeModelledResponse(const double& bandwidthFrequency,
	const double& dampingRatio)
{
	modelledResponse.resize(data.size());
	double a, b1, b2;
	ComputeCoefficients(bandwidthFrequency, dampingRatio, a, b1, b2);

	unsigned int i;
	for (i = 0; i < modelledResponse.size(); ++i)
	{
		modelledResponse[i].resize(data[i].size());
		modelledResponse[i][0] = data[i][0].response;
		modelledResponse[i][1] = data[i][1].response;

		unsigned int j;
		for (j = 2; j < data[i].size(); ++j)
			modelledResponse[i][j] = a * (data[i][j].input + 2 * data[i][j - 1].input + data[i][j - 2].input)
				- b1 * modelledResponse[i][j - 1] - b2 * modelledResponse[i][j - 2];
	}
}

void ResponseModeller::ComputeModelledResponse(const double& bandwidthFrequency)
{
	modelledResponse.resize(data.size());
	double a, b;
	ComputeCoefficients(bandwidthFrequency, a, b);

	unsigned int i;
	for (i = 0; i < modelledResponse.size(); ++i)
	{
		modelledResponse[i].resize(data[i].size());
		modelledResponse[i][0] = data[i][0].response;

		unsigned int j;
		for (j = 1; j < data[i].size(); ++j)
			modelledResponse[i][j] = a * (data[i][j].input + data[i][j - 1].input) - b * modelledResponse[i][j - 1];
	}
}

double ResponseModeller::ComputeModelError(const Eigen::VectorXd& parameters)
{
	if (modelType == ModelType::FirstOrder)
		ComputeModelledResponse(parameters(0));
	else if (modelType == ModelType::SecondOrder)
		ComputeModelledResponse(parameters(0), parameters(1));
	else// if (modelType == ModelType::NthOrder)
		ComputeModelledResponse(GetNumeratorCoefficients(parameters), GetDenominatorCoefficients(parameters));

	maximumError = 0.0;
	double totalError(0.0);
	unsigned int i;
	for (i = 0; i < data.size(); ++i)
	{
		unsigned int j;
		for (j = 0; j < data[i].size(); ++j)
		{
			double error(fabs(data[i][j].response - modelledResponse[i][j]));
			totalError += error;
			if (error > maximumError)
				maximumError = error;
		}
	}

	return totalError;
}

void ResponseModeller::ComputeModelledResponse(const std::vector<double>& sNum,
	const std::vector<double>& sDen)
{
	assert(sNum.size() == sDen.size() + 1);

	modelledResponse.resize(data.size());
	std::vector<double> numCoef;
	std::vector<double> denCoef;
	ComputeCoefficients(sNum, sDen, numCoef, denCoef);

	for (unsigned int i = 0; i < modelledResponse.size(); ++i)
	{
		modelledResponse[i].resize(data[i].size());
		for (unsigned int j = 0; j < sDen.size(); ++j)
			modelledResponse[i][j] = data[i][0].response;

		for (unsigned int j = sDen.size(); j < data[i].size(); ++j)
		{
			modelledResponse[i][j] = numCoef[0] * data[i][j].input;
			for (unsigned int k = 0; k < sDen.size(); ++k)
				modelledResponse[i][j] += numCoef[k + 1] * data[i][j - k - 1].input - denCoef[k] * modelledResponse[i][j - k - 1];
		}
	}
}

// Converts from continuous time coefficients to difference equation coefficients
void ResponseModeller::ComputeCoefficients(const std::vector<double>& sNum,
	const std::vector<double>& sDen, std::vector<double>& numCoef, std::vector<double>& denCoef)
{
	numCoef.resize(sNum.size());
	denCoef.resize(sDen.size());
	// TODO:  Implement
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

std::vector<double> ResponseModeller::GetNumeratorCoefficients(const Eigen::VectorXd& parameters)
{
	assert(parameters.size() % 2 == 1);
	std::vector<double> num((parameters.size() + 1) / 2);
	for (unsigned int i = 0; i < num.size(); ++i)
		num[i] = parameters(i);
	return num;
}

std::vector<double> ResponseModeller::GetDenominatorCoefficients(const Eigen::VectorXd& parameters)
{
	assert(parameters.size() % 2 == 1);
	std::vector<double> den((parameters.size() - 1) / 2);
	for (unsigned int i = 0; i < den.size(); ++i)
		den[i] = parameters(i + den.size() + 1);
	return den;
}
