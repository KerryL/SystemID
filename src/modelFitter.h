// File:  modelFitter.h
// Date:  7/10/2017
// Auth:  K. Loux
// Desc:  Model fitter class.

#ifndef MODEL_FITTER_H_
#define MODEL_FITTER_H_

// Standard C++ headers
#include <vector>

class ModelFitter
{
public:
	ModelFitter(const unsigned int& iterationLimit, const double& rolloverPoint)
		: iterationLimit(iterationLimit), rolloverPoint(rolloverPoint) {}

	struct Slice
	{
		Slice(const double& time, const double& input, const double& response)
			: time(time), input(input), response(response) {}

		double time;// [sec]
		double input;
		double response;
	};

	bool DetermineParameters(const std::vector<Slice>& data,
		double& bandwidthFrequency, double& sampleTime);

	bool DetermineParameters(const std::vector<Slice>& data,
		double& bandwidthFrequency, double& dampingRatio, double& sampleTime);

	unsigned int GetIterationCount() const { return iterationCount; }
	double GetMaximumError() const { return maximumError; }

private:
	static double ComputeSampleTime(const std::vector<Slice>& data);

	const unsigned int iterationLimit;
	const double rolloverPoint;

	unsigned int iterationCount;
	double maximumError = 0.0;
};

#endif// MODEL_FITTER_H_
