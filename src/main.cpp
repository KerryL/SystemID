// File:  main.cpp
// Date:  7/7/2017
// Auth:  K. Loux
// Desc:  Entry point for system identification tests.

// Local headers
#include "nelderMead.h"

// Standard C++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Slice
{
	Slice(const double& time, const double& input, const double& response)
		: time(time), input(input), response(response) {}

	double time;
	double input;
	double response;
};

bool ReadInputData(const std::string& fileName, std::vector<Slice>& data)
{
	std::ifstream inFile(fileName.c_str());
	if (!inFile.is_open() || !inFile.good())
	{
		std::cerr << "Failed to open '" << fileName << "' for input\n";
		return false;
	}

	// Skip first line for headers
	std::string line;
	std::getline(inFile, line);
	unsigned int lineCount(2);
	while (std::getline(inFile, line))
	{
		std::istringstream ss(line), tokenStream;
		std::string token;
		double time, input, response;

		std::getline(ss, token, ',');
		tokenStream.str(token);
		if ((tokenStream >> time).fail())
		{
			std::cerr << "Failed to parse time at row " << lineCount << std::endl;
			return false;
		}

		std::getline(ss, token, ',');
		tokenStream.clear();
		tokenStream.str(token);
		if ((tokenStream >> input).fail())
		{
			std::cerr << "Failed to parse input at row " << lineCount << std::endl;
			return false;
		}

		std::getline(ss, token, ',');
		tokenStream.clear();
		tokenStream.str(token);
		if ((tokenStream >> response).fail())
		{
			std::cerr << "Failed to parse response at row " << lineCount << std::endl;
			return false;
		}

		data.push_back(Slice(time, input, response));
		++lineCount;
	}

	return true;
}

double ComputeSampleTime(const std::vector<Slice>& data)
{
	// Assume equally spaced time steps
	return data.back().time / (data.size() - 1);
}

class ResponseModeller
{
public:
	ResponseModeller(const std::vector<Slice>& data, const double& sampleTime)
		: data(data), sampleTime(sampleTime) {}

	double ComputeModelError(const Eigen::VectorXd& parameters)
	{
		ComputeModelledResponse(parameters(0), parameters(1));
		double error(0.0);
		unsigned int i;
		for (i = 0; i < data.size(); ++i)
			error += fabs(data[i].response - modelledResponse[i]);
		return error;
	}

private:
	const std::vector<Slice>& data;
	const double sampleTime;

	std::vector<double> modelledResponse;
	void ComputeModelledResponse(const double& bandwidthFrequency,
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

	void ComputeCoefficients(const double& bandwidthFrequency,
		const double& dampingRatio, double& a, double& b1, double& b2)
	{
		const double denominator(4 + 4 * dampingRatio * bandwidthFrequency * sampleTime
			+ bandwidthFrequency * bandwidthFrequency * sampleTime * sampleTime);
		a = (bandwidthFrequency * bandwidthFrequency * sampleTime * sampleTime) / denominator;
		b1 = (2 * bandwidthFrequency * bandwidthFrequency * sampleTime * sampleTime - 8) / denominator;
		b2 = (4 - 4 * dampingRatio * bandwidthFrequency * sampleTime
			+ bandwidthFrequency * bandwidthFrequency * sampleTime * sampleTime) / denominator;
	}
};

void DetermineParameters(const std::vector<Slice>& data,
	double& bandwidthFrequency, double& dampingRatio, double& sampleTime)
{
	sampleTime = ComputeSampleTime(data);

	ResponseModeller modeller(data, sampleTime);
	auto objectiveFunction = std::bind(&ResponseModeller::ComputeModelError, modeller, std::placeholders::_1);
	NelderMead<2> optimization(objectiveFunction, 1000);
	Eigen::VectorXd initialGuess(2);
	initialGuess(0) = 1.0;
	initialGuess(1) = 1.0;
	optimization.SetInitialGuess(initialGuess);

	Eigen::VectorXd parameters(optimization.Optimize());
	bandwidthFrequency = parameters(0);
	dampingRatio = parameters(1);

	std::cout << "Optimization complete after " << optimization.GetIterationCount() << std::endl;
}

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		std::cout << "Usage:  " << argv[0] << " <input file>\n"
			<< "    Input file must be formatted into three columns separated by ',':\n"
			<< "    Time, Input, Response\n" << std::endl;
		return 1;
	}

	std::vector<Slice> data;
	if (!ReadInputData(argv[1], data))
		return 1;

	double bandwidthFrequency, dampingRatio, sampleTime;
	DetermineParameters(data, bandwidthFrequency, dampingRatio, sampleTime);
	std::cout << "Sample time = " << sampleTime << std::endl;
	std::cout << "Bandwidth frequency = " << bandwidthFrequency / 2.0 / M_PI << std::endl;
	std::cout << "Damping ratio = " << dampingRatio << std::endl;

	return 0;
}

