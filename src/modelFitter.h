// File:  modelFitter.h
// Date:  7/10/2017
// Auth:  K. Loux
// Desc:  Model fitter class.

#ifndef MODEL_FITTER_H_
#define MODEL_FITTER_H_

// Local headers
#include "nelderMead.h"
#include "responseModeller.h"
#include "slice.h"

// Standard C++ headers
#include <vector>

class ModelFitter
{
public:
	ModelFitter(const unsigned int& iterationLimit, const double& rolloverPoint)
		: iterationLimit(iterationLimit), rolloverPoint(rolloverPoint) {}

	bool DetermineParameters(const std::vector<std::vector<Slice>>& data,
		double& bandwidthFrequency, double& sampleTime);

	bool DetermineParameters(const std::vector<std::vector<Slice>>& data,
		double& bandwidthFrequency, double& dampingRatio, double& sampleTime);

	bool DetermineParameters(const std::vector<std::vector<Slice>>& data,
		std::vector<double>& numerator, std::vector<double>& denominator,
		const unsigned int& order, double& sampleTime);

	unsigned int GetIterationCount() const { return iterationCount; }
	double GetMaximumError() const { return maximumError; }

private:
	static double ComputeSampleTime(const std::vector<std::vector<Slice>>& data);

	const unsigned int iterationLimit;
	const double rolloverPoint;

	unsigned int iterationCount;
	double maximumError = 0.0;

	std::vector<double> BuildTriangularCoefficients(const unsigned int& order);
};

#endif// MODEL_FITTER_H_
