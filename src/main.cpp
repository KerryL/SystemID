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

class ResponseModeller
{
public:
	ResponseModeller(const std::vector<Slice>& data) : data(data) {}

	double ComputeModelError(const Eigen::VectorXd& parameters)
	{
		return 0.0;
	}

private:
	const std::vector<Slice>& data;
};

void DetermineParameters(const std::vector<Slice>& data,
	double& bandwidthFrequency, double& dampingRatio)
{
	ResponseModeller modeller(data);
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

	double bandwidthFrequency, dampingRatio;
	DetermineParameters(data, bandwidthFrequency, dampingRatio);
	std::cout << "Bandwidth frequency = " << bandwidthFrequency / 2.0 / M_PI << std::endl;
	std::cout << "Damping ratio = " << dampingRatio << std::endl;

	return 0;
}

