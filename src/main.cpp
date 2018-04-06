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

bool ParseValue(const std::string& token, double& value)
{
	std::istringstream ss(token);
	return !(ss >> value).fail();
}

bool ReadInputData(const std::string& fileName, std::vector<ModelFitter::Slice>& data,
	const int& commandColumn, const int& responseColumn)
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
		std::istringstream ss(line);
		std::string token;
		double time, input, response;

		int column(0), assignedDataCount(0);
		while (std::getline(ss, token, ','))
		{
			if (column == 0)
			{
				if (!ParseValue(token, time))
				{
					std::cerr << "Failed to parse time at row " << lineCount << std::endl;
					return false;
				}
				++assignedDataCount;
			}
			else if (column == commandColumn)
			{
				if (!ParseValue(token, input))
				{
					std::cerr << "Failed to parse input at row " << lineCount << std::endl;
					return false;
				}
				++assignedDataCount;
			}
			else if (column == responseColumn)
			{
				if (!ParseValue(token, response))
				{
					std::cerr << "Failed to parse response at row " << lineCount << std::endl;
					return false;
				}
				++assignedDataCount;
			}

			if (assignedDataCount == 3)
				break;

			++column;
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

void PrintUsage(const std::string& appName)
{
	std::cout << "Usage:  " << appName << " [--rollover <rolloverPoint>] [--inCol <index> --rspCol <index>] <input file>\n"
		<< "    Input file must be formatted into three columns separated by ',':\n"
		<< "    Time (must be in seconds), Input, Response\n"
		<< "    If data file is not in this format or has more than three columns,\n"
		<< "    --inCol and --rspCol argumetns must also be provided.\n"
		<< "    If rollover argument is omitted, no rollover correction is performed." << std::endl;
}

bool ProcessRolloverArgument(const std::string& arg, double& rollover)
{
	std::istringstream ss(arg);
	if ((ss >> rollover).fail())
	{
		std::cerr << "Invalid rollover specification:  '" << arg << "'\n";
		return false;
	}

	if (rollover < 0.0)
	{
		std::cerr << "Rollover point must be positive\n";
		return false;
	}

	return true;
}

bool ProcessColumnArgument(const std::string& arg, int& column)
{
	std::istringstream ss(arg);
	if ((ss >> column).fail())
	{
		std::cerr << "Invalid column specification:  '" << arg << "'\n";
		return false;
	}

	if (column == 0)
	{
		std::cerr << "Column specification must be strictly positive (zero is reserved for time data)\n";
		return false;
	}
	else if (column < 0)
	{
		std::cerr << "Column specification must be strictly positive\n";
		return false;
	}

	return true;
}

bool ProcessArguments(const std::vector<std::string>& args,
	std::string& inputFileName, double& rolloverPoint, int& commandColumn, int& responseColumn)
{
	rolloverPoint = 0.0;
	commandColumn = 1;
	responseColumn = 2;

	inputFileName = args.back();// Always the last argument
	if (args.size() == 2)
		return true;

	unsigned int argIndex;
	for (argIndex = 1; argIndex < args.size() - 1; argIndex += 2)
	{
		if (argIndex == args.size() - 2)
		{
			std::cerr << "Unexpected input arguments\n";
			PrintUsage(args.front());
			return false;
		}

		if (args[argIndex].compare("--rollover") == 0)
		{
			if (!ProcessRolloverArgument(args[argIndex + 1], rolloverPoint))
				return false;
		}
		else if (args[argIndex].compare("--inCol") == 0)
		{
			if (!ProcessColumnArgument(args[argIndex + 1], commandColumn))
				return false;
		}
		else if (args[argIndex].compare("--rspCol") == 0)
		{
			if (!ProcessColumnArgument(args[argIndex + 1], responseColumn))
				return false;
		}
		else
		{
			std::cerr << "Unexpected argument '" << args[argIndex] << "'\n";
			PrintUsage(args.front());
			return false;
		}
	}

	return true;
}

int main(int argc, char *argv[])
{
	std::string inputFileName;
	double rolloverPoint;
	int inputColumn;
	int responseColumn;

	std::vector<std::string> args(argv, argv + argc);
	if (!ProcessArguments(args, inputFileName, rolloverPoint, inputColumn, responseColumn))
		return 1;

	std::vector<ModelFitter::Slice> data;
	if (!ReadInputData(inputFileName, data, inputColumn, responseColumn))
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

