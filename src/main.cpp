// File:  main.cpp
// Date:  7/7/2017
// Auth:  K. Loux
// Desc:  Entry point for system identification tests.

// Local headers
#include "modelFitter.h"
#include "expressionTree.h"

// Standard C++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <set>

bool ParseValue(const std::string& token, double& value)
{
	std::istringstream ss(token);
	return !(ss >> value).fail();
}

bool ReadInputDataFile(const std::string& fileName, std::vector<Slice>& data,
	const int& commandColumn, const int& responseColumn, const double& timeFactor)
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
		double time(0.0), input, response;

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

		data.push_back(Slice(time * timeFactor, input, response));
		++lineCount;
	}

	return true;
}

bool ReadInputData(const std::vector<std::string>& fileName, std::vector<std::vector<Slice>>& data,
	const int& commandColumn, const int& responseColumn, const double& timeFactor)
{
	unsigned int i(0);
	for (const auto& f : fileName)
	{
		if (!ReadInputDataFile(f, data[i], commandColumn, responseColumn, timeFactor))
			return false;
		++i;
	}

	return true;
}

void Unwind(std::vector<Slice>& data, const double& rollover)
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
	std::cout
		<< "Usage:  " << appName << " [--rollover <rolloverPoint>] [--order <order>]\n"
		<< "            [--time-factor <factor>] [--inCol <index> --rspCol <index>]\n"
		<< "            [--bwGuess <bandwidth>] <first input file> ...\n\n"

		<< "    Input files must be formatted into three columns separated by ',':\n"
		<< "        Time, Input, Response\n"
		<< "        Time must be in seconds or the time-factor argument must be supplied\n"
		<< "        such that time * factor results units of seconds.\n"
		<< "        If data files are not in this format or have more than three columns,\n"
		<< "        --inCol and --rspCol argumetns must also be provided.  Column 0 must\n"
		<< "        contain the time data.  All input files must be in the same format.\n\n"

		<< "    If rollover argument is omitted, no rollover correction is performed.\n\n"

		<< "    If time-factor argument is omitted, no scaling is applied to the time data.\n\n"

		<< "    If order argument is supplied, an attempt is made to fit an arbitary\n"
		<< "    transfer function of the specified order (n zeros and n poles).\n"

		<< "    The solver requires an initial guess for the parameters.  Convergence is\n"
		<< "    sometimes sensitive to the choice of the initial guess.  If you run into\n"
		<< "    convergence issues, you can suggest a reasonable starting point with the\n"
		<< "    bwGuess argument.  Units for bwGuess are Hertz.  The damping ratio for\n"
		<< "    the initial guess is always zero.  If omitted, the initial bandwidth is\n"
		<< "    about 10 Hz.  Does not apply to user-specified models.\n\n"

		<< "    The solver can be used to adjust parameters for user-specified models.  To\n"
		<< "    specify a model, use the --model_num and --model_den arguments to give the\n"
		<< "    numerator and denominators of the continuous-time transfer function.  Use\n"
		<< "    'p' to prefix parameter names.  For example, to fit a second-order high-\n"
		<< "    pass filter, specify arguments like this:\n"
		<< "      --modelNum s^2 --modelDen s^2+2*pOmega=30*pZeta=0.7*s+pOmega^2\n"
		<< "    Do not use spaces within the model arguments.  All multiplications must be\n"
		<< "    specified explicitly (i.e. '5s' will fail where '5*s' will be accepted).\n"
		<< "    Initial guesses can be specified for all parameters using 'pName=value'\n"
		<< "    syntax, for example 'pOmega=30.'  If no initial guess is specified, the\n"
		<< "    initial value for that parameter will be 1.\n\n"
		<< std::endl;
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

bool ProcessOrderArgument(const std::string& arg, unsigned int& order)
{
	std::istringstream ss(arg);
	if ((ss >> order).fail())
	{
		std::cerr << "Invalid order specification:  '" << arg << "'\n";
		return false;
	}

	if (order == 0)
	{
		std::cerr << "Order must be strictly positive\n";
		return false;
	}

	return true;
}

bool ProcessTimeFactorArgument(const std::string& arg, double& timeFactor)
{
	std::istringstream ss(arg);
	if ((ss >> timeFactor).fail())
	{
		std::cerr << "Invalid time factor specification:  '" << arg << "'\n";
		return false;
	}

	if (timeFactor <= 0.0)
	{
		std::cerr << "Time factor must be strictly positive\n";
		return false;
	}

	return true;
}

bool ProcessBandwidthArgument(const std::string& arg, double& bandwidth)
{
	std::istringstream ss(arg);
	if ((ss >> bandwidth).fail())
	{
		std::cerr << "Invalid bandwidth specification:  '" << arg << "'\n";
		return false;
	}

	if (bandwidth <= 0.0)
	{
		std::cerr << "Bandwidth must be strictly positive\n";
		return false;
	}

	bandwidth *= 2.0 * M_PI;

	return true;
}

bool IsOperatorOrParenthese(const char& c)
{
	return c == '-' || c == '+' || c == '*' ||
		c == '/' || c == '^' || c == '(' || c == ')';
}

bool IsNumber(const std::string& s)
{
	std::istringstream ss(s);
	double dummy;
	return !(ss >> dummy).fail() && ss.eof();
}

std::string::size_type GetNextEndPosition(const std::string& s)
{
	for (std::string::size_type p = 1; p < s.length(); ++p)
	{
		if (IsOperatorOrParenthese(s[p]))
			return p;
	}

	return std::string::npos;
}

struct Parameter
{
	std::string name;
	double initialGuess = 1.0;

	bool operator<(const Parameter& other) const
	{
		return name < other.name;
	}
};

bool ParseParameter(const std::string& s, Parameter& p)
{
	const auto equalPos(s.find('='));
	p.name = s.substr(0, equalPos);
	if (equalPos == std::string::npos)
		return true;

	std::istringstream ss(s.substr(equalPos + 1));
	if ((ss >> p.initialGuess).fail())
	{
		std::cerr << "Failed to parse initial guess from " << s << '\n';
		return false;
	}

	return true;
}

std::string StripInitialGuesses(std::string model)
{
	std::string::size_type equalPos(0);
	while (equalPos = model.find('=', equalPos), equalPos != std::string::npos)
	{
		const auto nextEnd(GetNextEndPosition(model.substr(equalPos)));
		std::string temp(model.substr(0, equalPos));
		if (nextEnd == std::string::npos)
			model = temp;
		else
			model = temp + model.substr(equalPos + nextEnd);
	}

	return model;
}

bool ExpressionIsValid(std::string model, const std::set<Parameter>& params)
{
	for (const auto& p : params)
		model = Utilities::ReplaceAllOccurrences(model, p.name, "1");

	ExpressionTree e;
	std::string expression;
	return e.Solve(model, expression).empty();
}

bool ProcessModelArgument(const std::string& arg, std::string& model, std::set<Parameter>& params)
{
	if (arg.empty())// TODO:  Trim first?
		return false;

	for (std::string::size_type p = 0; p < arg.length();)
	{
		const auto sLen(GetNextEndPosition(arg.substr(p)));
		const auto subStr(arg.substr(p, sLen));

		if (subStr.front() == 'p')
		{
			Parameter param;
			if (!ParseParameter(subStr, param))
				return false;

			auto it(params.begin());
			for (; it != params.end(); ++it)
			{
				if (it->name.compare(param.name) == 0)
					break;
			}

			if (it == params.end())
				params.insert(param);
		}
		else if (subStr.front() == '(')
		{
			++p;
			continue;
		}
		
		p += subStr.length() + 1;
	}

	model = StripInitialGuesses(arg);
	return ExpressionIsValid(model, params);
}

struct Configuration
{
	std::vector<std::string> inputFileNames;

	double rolloverPoint = 0.0;

	int commandColumn = 1;
	int responseColumn = 2;

	unsigned int order = 0;
	std::string modelNum;
	std::string modelDen;
	std::set<Parameter> modelParams;

	double timeFactor = 1.0;// [sec/input units]
	double bandwidthGuess = 30.0;// [rad/sec]
};

bool ProcessArguments(const std::vector<std::string>& args, Configuration& configuration)
{
	unsigned int argIndex;
	for (argIndex = 1; argIndex < args.size(); ++argIndex)
	{
		if (args[argIndex].substr(0, 2).compare("--") == 0 &&
			(argIndex == args.size() - 2 ||// Shouldn't have non-file-name argument identifier as second-to-last argument
			!configuration.inputFileNames.empty()))// After we start adding file names, don't allow non-file-name arguments
		{
			std::cerr << "Unexpected input arguments\n";
			return false;
		}

		if (args[argIndex].compare("--rollover") == 0)
		{
			if (!ProcessRolloverArgument(args[argIndex + 1], configuration.rolloverPoint))
				return false;
			++argIndex;
		}
		else if (args[argIndex].compare("--inCol") == 0)
		{
			if (!ProcessColumnArgument(args[argIndex + 1], configuration.commandColumn))
				return false;
			++argIndex;
		}
		else if (args[argIndex].compare("--rspCol") == 0)
		{
			if (!ProcessColumnArgument(args[argIndex + 1], configuration.responseColumn))
				return false;
			++argIndex;
		}
		else if (args[argIndex].compare("--order") == 0)
		{
			if (!ProcessOrderArgument(args[argIndex + 1], configuration.order))
				return false;
			++argIndex;
		}
		else if (args[argIndex].compare("--time-factor") == 0)
		{
			if (!ProcessTimeFactorArgument(args[argIndex + 1], configuration.timeFactor))
				return false;
			++argIndex;
		}
		else if (args[argIndex].compare("--bwGuess") == 0)
		{
			if (!ProcessBandwidthArgument(args[argIndex + 1], configuration.bandwidthGuess))
				return false;
			++argIndex;
		}
		else if (args[argIndex].compare("--modelNum") == 0)
		{
			if (!ProcessModelArgument(args[argIndex + 1], configuration.modelNum, configuration.modelParams))
				return false;
			++argIndex;
		}
		else if (args[argIndex].compare("--modelDen") == 0)
		{
			if (!ProcessModelArgument(args[argIndex + 1], configuration.modelDen, configuration.modelParams))
				return false;
			++argIndex;
		}
		else
			configuration.inputFileNames.push_back(args[argIndex]);
	}

	return !configuration.inputFileNames.empty() &&
		configuration.modelNum.empty() == configuration.modelDen.empty();
}

std::string AssembleSExpression(const std::vector<double> coefficients)
{
	std::ostringstream ss;
	for (unsigned int c = 0; c < coefficients.size(); ++c)
	{
		if (!ss.str().empty())
		{
			if (coefficients[c] > 0)
				ss << " + ";
			else
				ss << " - ";
		}

		if (fabs(coefficients[c]) != 1.0)
			ss << fabs(coefficients[c]);
		const unsigned int power(coefficients.size() - c - 1);
		if (power == 1)
			ss << "s";
		else if (power > 1)
			ss << "s^" << power;
	}

	return ss.str();
}

void PrintTransferFunction(const std::vector<double> numerator, const std::vector<double>& denominator)
{
	assert(numerator.size() == denominator.size() + 1);

	std::vector<double> entireDenominator(denominator);
	entireDenominator.insert(entireDenominator.begin(), 1.0);

	const auto numString(AssembleSExpression(numerator));
	const auto denString(AssembleSExpression(entireDenominator));

	std::cout << "  Transfer Function = " << numString << std::endl;
	std::cout << "                      " << std::string(std::max(numString.length(), denString.length()), '-') << std::endl;
	std::cout << "                      " << denString << std::endl;
}

int main(int argc, char *argv[])
{
	Configuration configuration;

	std::vector<std::string> args(argv, argv + argc);
	if (!ProcessArguments(args, configuration))
	{
		PrintUsage(argv[0]);
		return 1;
	}

	std::vector<std::vector<Slice>> data(configuration.inputFileNames.size());
	if (!ReadInputData(configuration.inputFileNames, data, configuration.commandColumn, configuration.responseColumn, configuration.timeFactor))
		return 1;

	const unsigned int recordCount([&data]()
	{
		unsigned int total(0);
		for (const auto& d : data)
			total += d.size();
		return total;
	}());
	std::cout << "Found " << recordCount << " records in " << configuration.inputFileNames.size() << " files" << std::endl;

	const unsigned int iterationLimit(100000);
	ModelFitter fitter(iterationLimit, configuration.rolloverPoint, configuration.bandwidthGuess);
	if (configuration.rolloverPoint > 0.0)
	{
		std::cout << "Unwinding data at " << configuration.rolloverPoint << std::endl;
		for (auto& d : data)
			Unwind(d, configuration.rolloverPoint);
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

	if (configuration.order > 0)
	{
		std::vector<double> nOrderNum;
		std::vector<double> nOrderDen;

		std::cout << "\nFitting model for transfer function of order " << configuration.order << "..." << std::endl;
		if (fitter.DetermineParameters(data, nOrderNum, nOrderDen, configuration.order, sampleTime))
		{
			std::cout << "  Optimization complete after " << fitter.GetIterationCount() << " iterations" << std::endl;
			std::cout << "  Sample time = " << sampleTime << " sec" << std::endl;
			PrintTransferFunction(nOrderNum, nOrderDen);
			std::cout << "  Maximum error = " << fitter.GetMaximumError() << std::endl;
		}
		else
			std::cout << "  Solution failed to converge after " << iterationLimit << " iterations" << std::endl;
	}

	if (!configuration.modelNum.empty() && !configuration.modelDen.empty() && !configuration.modelParams.empty())
	{
		std::cout << "\nFitting parameters for model given as\n\n    "
			<< configuration.modelNum << "\n    "
			<< std::string(std::max(configuration.modelNum.length(), configuration.modelDen.length()), '-') << "\n    "
			<< configuration.modelDen << '\n' << std::endl;
		std::cout << "Model includes " << configuration.modelParams.size() << " parameters" << std::endl;
		std::map<std::string, double> modelParams;
		for (const auto& p : configuration.modelParams)
			modelParams[p.name] = p.initialGuess;

		if (fitter.DetermineParameters(data, configuration.modelNum, configuration.modelDen, modelParams, sampleTime))
		{
			std::cout << "  Optimization complete after " << fitter.GetIterationCount() << " iterations" << std::endl;
			std::cout << "  Sample time = " << sampleTime << " sec" << std::endl;
			auto it(modelParams.begin());
			for (; it != modelParams.end(); ++it)
				std::cout << "  " << it->first << " = " << it->second << std::endl;
			std::cout << "  Maximum error = " << fitter.GetMaximumError() << std::endl;
		}
		else
			std::cout << "  Solution failed to converge after " << iterationLimit << " iterations" << std::endl;
	}

	return 0;
}
