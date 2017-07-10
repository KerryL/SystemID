// File:  responseModeller.h
// Date:  7/10/2017
// Auth:  K. Loux
// Desc:  Model fitter class.

#ifndef RESPONSE_MODELLER_H_
#define RESPONSE_MODELLER_H_

// Local headers
#include "modelFitter.h"

// Eigen headers
#include <Eigen/Eigen>

class ResponseModeller
{
public:
	enum class ModelType
	{
		FirstOrder,
		SecondOrder
	};

	ResponseModeller(const std::vector<ModelFitter::Slice>& data,
		const double& sampleTime, const double& rolloverPoint, const ModelType& modelType)
		: data(data), sampleTime(sampleTime), rolloverPoint(rolloverPoint), modelType(modelType) {}

	double ComputeModelError(const Eigen::VectorXd& parameters);
	double GetMaximumError() const { return maximumError; }

private:
	const std::vector<ModelFitter::Slice>& data;
	const double sampleTime;
	const double rolloverPoint;
	ModelType modelType;
	double maximumError = 0.0;

	std::vector<double> modelledResponse;

	void ComputeModelledResponse(const double& bandwidthFrequency,
		const double& dampingRatio);
	void ComputeModelledResponse(const double& bandwidthFrequency);

	void ComputeCoefficients(const double& bandwidthFrequency,
		const double& dampingRatio, double& a, double& b1, double& b2);
	void ComputeCoefficients(const double& bandwidthFrequency,
		double& a, double& b);

	double Unwind(const double& value) const;
};

#endif// RESPONSE_MODELLER_H_
