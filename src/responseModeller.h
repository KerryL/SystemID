// File:  responseModeller.h
// Date:  7/10/2017
// Auth:  K. Loux
// Desc:  Model fitter class.

#ifndef RESPONSE_MODELLER_H_
#define RESPONSE_MODELLER_H_

// Local headers
#include "slice.h"

// Eigen headers
#include <Eigen/Eigen>

// Standard C++ headers
#include <vector>

class ResponseModeller
{
public:
	enum class ModelType
	{
		FirstOrder,
		SecondOrder,
		NthOrder,
		UserModel
	};

	ResponseModeller(const std::vector<std::vector<Slice>>& data,
		const double& sampleTime, const double& rolloverPoint, const ModelType& modelType)
		: data(data), sampleTime(sampleTime), rolloverPoint(rolloverPoint), modelType(modelType) {}
	ResponseModeller(const std::vector<std::vector<Slice>>& data,
		const double& sampleTime, const double& rolloverPoint, const std::string& userNum,
		const std::string& userDen, const std::vector<std::string>& userParams)
		: data(data), sampleTime(sampleTime), rolloverPoint(rolloverPoint), modelType(ModelType::UserModel),
		userNum(userNum), userDen(userDen), userParams(userParams) {}

	double ComputeModelError(const Eigen::VectorXd& parameters);
	double GetMaximumError() const { return maximumError; }

private:
	const std::vector<std::vector<Slice>>& data;
	const double sampleTime;
	const double rolloverPoint;
	const ModelType modelType;
	const std::string userNum;
	const std::string userDen;
	const std::vector<std::string> userParams;
	double maximumError = 0.0;

	std::vector<std::vector<double>> modelledResponse;

	void ComputeModelledResponse(const double& bandwidthFrequency,
		const double& dampingRatio);
	void ComputeModelledResponse(const double& bandwidthFrequency);
	void ComputeModelledResponse(const std::vector<double>& sNum, const std::vector<double>& sDen);

	void ComputeCoefficients(const double& bandwidthFrequency,
		const double& dampingRatio, double& a, double& b1, double& b2);
	void ComputeCoefficients(const double& bandwidthFrequency,
		double& a, double& b);

	std::vector<double> GetNumeratorCoefficients(const Eigen::VectorXd& parameters);
	std::vector<double> GetDenominatorCoefficients(const Eigen::VectorXd& parameters);
};

#endif// RESPONSE_MODELLER_H_
