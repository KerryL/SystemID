// File:  main.cpp
// Date:  7/7/2017
// Auth:  K. Loux
// Desc:  Entry point for system identification tests.

// Local headers
#include "modelFitter.h"

// Standard C++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

bool ReadInputData(const std::string& fileName, std::vector<ModelFitter::Slice>& data)
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

		data.push_back(ModelFitter::Slice(time, input, response));
		++lineCount;
	}

	return true;
}

void Unwind(std::vector<ModelFitter::Slice>& data, const double& rollover)
{
	unsigned int i;
	for (i = 1; i < data.size(); ++i)
	{
		if (fabs(data[i].input - rollover - data[i - 1].input) < fabs(data[i].input - data[i - 1].input))
		{
			unsigned int j;
			for (j = i; j < data.size(); ++j)
				data[j].input -= rollover;
		}
		else if (fabs(data[i].input + rollover - data[i - 1].input) < fabs(data[i].input - data[i - 1].input))
		{
			unsigned int j;
			for (j = i; j < data.size(); ++j)
				data[j].input += rollover;
		}

		if (fabs(data[i].response - rollover - data[i - 1].response) < fabs(data[i].response - data[i - 1].response))
		{
			unsigned int j;
			for (j = i; j < data.size(); ++j)
				data[j].response -= rollover;
		}
		else if (fabs(data[i].response + rollover - data[i - 1].response) < fabs(data[i].response - data[i - 1].response))
		{
			unsigned int j;
			for (j = i; j < data.size(); ++j)
				data[j].response += rollover;
		}
	}
}

bool ProcessArguments(const std::vector<std::string>& args,
	std::string& inputFileName, double& rolloverPoint)
{
	if ((args.size() != 2 && args.size() != 4) ||
		(args.size() == 4 && args[1].compare("--rollover") != 0))
	{
		std::cout << "Usage:  " << args.front() << " [--rollover <rolloverPoint>] <input file>\n"
			<< "    Input file must be formatted into three columns separated by ',':\n"
			<< "    Time (must be in seconds), Input, Response\n"
			<< "    If rollover argument is omitted, no rollover correction is performed." << std::endl;
		return false;
	}

	inputFileName = args.back();// Always the last argument

	rolloverPoint = 0.0;
	if (args.size() == 4)
	{
		std::istringstream ss(args[2]);
		if ((ss >> rolloverPoint).fail())
		{
			std::cerr << "Invalid rollover specification:  '" << args[2] << "'\n";
			return false;
		}

		if (rolloverPoint < 0.0)
		{
			std::cerr << "Rollover point must be positive\n";
			return false;
		}
	}

	return true;
}

int main(int argc, char *argv[])
{
	std::string inputFileName;
	double rolloverPoint;

	std::vector<std::string> args(argv, argv + argc);
	if (!ProcessArguments(args, inputFileName, rolloverPoint))
		return 1;

	std::vector<ModelFitter::Slice> data;
	if (!ReadInputData(inputFileName, data))
		return 1;

	std::cout << "Found " << data.size() << " records in " << inputFileName << std::endl;

	const unsigned int iterationLimit(1000);
	ModelFitter fitter(iterationLimit, rolloverPoint);
	if (rolloverPoint > 0.0)
	{
		std::cout << "Unwinding data at " << rolloverPoint << std::endl;
		Unwind(data, rolloverPoint);
	}

	double bandwidthFrequency, sampleTime;

	std::cout << "\nFitting first-order model..." << std::endl;
	if (fitter.DetermineParameters(data, bandwidthFrequency, sampleTime))
	{
		std::cout << "  Optimization complete after " << fitter.GetIterationCount() << " iterations" << std::endl;
		std::cout << "  Sample time = " << sampleTime << " sec" << std::endl;
		std::cout << "  Bandwidth frequency = " << bandwidthFrequency / 2.0 / M_PI << " Hz" << std::endl;
		std::cout << "  Maximum error = " << fitter.GetMaximumError() << std::endl;
	}
	else
		std::cout << "  Solution failed to converge after " << iterationLimit << " iterations" << std::endl;

	double dampingRatio;
	std::cout << "\nFitting second-order model..." << std::endl;
	if (fitter.DetermineParameters(data, bandwidthFrequency, dampingRatio, sampleTime))
	{
		std::cout << "  Optimization complete after " << fitter.GetIterationCount() << " iterations" << std::endl;
		std::cout << "  Sample time = " << sampleTime << " sec" << std::endl;
		std::cout << "  Bandwidth frequency = " << bandwidthFrequency / 2.0 / M_PI << " Hz" << std::endl;
		std::cout << "  Damping ratio = " << dampingRatio << std::endl;
		std::cout << "  Maximum error = " << fitter.GetMaximumError() << std::endl;
	}
	else
		std::cout << "  Solution failed to converge after " << iterationLimit << " iterations" << std::endl;

	return 0;
}

