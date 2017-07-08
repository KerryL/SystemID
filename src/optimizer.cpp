// File:  optimizer.cpp
// Date:  8/4/2016
// Auth:  K. Loux
// Desc:  Base class for optimization algorithms.

// Local headers
#include "optimizer.h"

//==========================================================================
// Class:			Optimizer
// Function:		Optimizer
//
// Description:		Constructor for Optimizer class.
//
// Input Arguments:
//		objectiveFunction	= ObjectiveFunction
//		iterationLimit		= const unsigned int&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//==========================================================================
Optimizer::Optimizer(ObjectiveFunction objectiveFunction,
	const unsigned int& iterationLimit) : objectiveFunction(objectiveFunction),
	iterationLimit(iterationLimit)
{
}
