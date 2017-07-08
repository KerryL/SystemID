// File:  optimizer.h
// Date:  8/4/2016
// Auth:  K. Loux
// Desc:  Base class for optimization algorithms.

#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

// Standard C++ headers
#include <functional>

// Eigen headers
#include <Eigen/Eigen>

class Optimizer
{
public:
	typedef std::function<double(const Eigen::VectorXd&)> ObjectiveFunction;
	Optimizer(ObjectiveFunction objectiveFunction, const unsigned int& iterationLimit);
	virtual ~Optimizer() {}

	virtual Eigen::VectorXd Optimize() const = 0;

protected:
	const ObjectiveFunction objectiveFunction;
	const unsigned int& iterationLimit;
};

#endif// OPTIMIZER_H_