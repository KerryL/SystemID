/*=============================================================================
                                   LibPlot2D
                       Copyright Kerry R. Loux 2011-2016

                  This code is licensed under the GPLv2 License
                    (http://opensource.org/licenses/GPL-2.0).
=============================================================================*/

// File:  expressionTree.cpp
// Date:  5/6/2011
// Auth:  K. Loux
// Desc:  Handles user-specified mathematical operations on datasets.

// Local headers
#include "expressionTree.h"
#include "idMath.h"

//=============================================================================
// Class:			ExpressionTree
// Function:		Constant Declarations
//
// Description:		Constant declarations for ExpressionTree class.
//
// Input Arguments:
//		None
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
const unsigned int ExpressionTree::mPrintfPrecision = 15;

//=============================================================================
// Class:			ExpressionTree
// Function:		Solve
//
// Description:		Main solving method for the tree.
//
// Input Arguments:
//		expression			= std::string containing the expression to parse
//
// Output Arguments:
//		solvedExpression	= std::string& containing the evaluated data
//
// Return Value:
//		std::string, empty for success, error string if unsuccessful
//
//=============================================================================
std::string ExpressionTree::Solve(std::string expression,
	std::string &solvedExpression)
{
	if (!ParenthesesBalanced(expression))
		return "Imbalanced parentheses!";

	std::string errorString;
	errorString = ParseExpression(expression).c_str();

	if (!errorString.empty())
		return errorString;

	errorString = EvaluateExpression(solvedExpression);

	return errorString;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		ParenthesesBalanced
//
// Description:		Checks to see if the expression has balanced parentheses.
//
// Input Arguments:
//		expression	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		bool, true if parentheses are balanced, false otherwise
//
//=============================================================================
bool ExpressionTree::ParenthesesBalanced(const std::string &expression) const
{
	unsigned int leftCount(0), rightCount(0);
	auto location = expression.find("(");

	while (location != std::string::npos)
	{
		++leftCount;
		location = expression.find("(", location + 1);
	}

	location = expression.find(")");

	while (location != std::string::npos)
	{
		++rightCount;
		location = expression.find(")", location + 1);
	}

	if (leftCount != rightCount)
		return false;

	return true;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		ParseExpression
//
// Description:		Parses the expression and produces a queue of Reverse
//					Polish Notation values and operations.  Implements the
//					shunting-yard algorithm as described by Wikipedia.
//
// Input Arguments:
//		expression	= const std::string& to be parsed
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string containing error descriptions or an empty string on success
//
//=============================================================================
std::string ExpressionTree::ParseExpression(const std::string &expression)
{
	std::stack<std::string> operatorStack;
	unsigned int i, advance;
	bool lastWasOperator(true);
	std::string errorString;

	for (i = 0; i < expression.length(); ++i)
	{
		if (expression.substr(i, 1).empty())
			continue;

		errorString = ParseNext(expression.substr(i), lastWasOperator,
			advance, operatorStack);
		if (!errorString.empty())
			return errorString;
		i += advance - 1;
	}

	if (!EmptyStackToQueue(operatorStack))
		errorString = "Imbalanced parentheses!";

	return errorString;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		ParseNext
//
// Description:		Parses the expression and processes the next item.
//
// Input Arguments:
//		expression	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string containing any errors
//
//=============================================================================
std::string ExpressionTree::ParseNext(const std::string &expression,
	bool &lastWasOperator, unsigned int &advance,
	std::stack<std::string> &operatorStack)
{
	bool thisWasOperator(false);
	if (NextIsNumber(expression, &advance, lastWasOperator))
		mOutputQueue.push(expression.substr(0, advance));
	else if (NextIsS(expression, &advance))
		mOutputQueue.push(expression.substr(0, advance));
	else if (NextIsOperator(expression, &advance))
	{
		ProcessOperator(operatorStack, expression.substr(0, advance));
		thisWasOperator = true;
	}
	else if (expression[0] == '(')
	{
		if (!lastWasOperator)
			operatorStack.push("*");
		operatorStack.push(std::string(1, expression[0]));
		advance = 1;
		thisWasOperator = true;
	}
	else if (expression[0] == ')')
	{
		ProcessCloseParenthese(operatorStack);
		advance = 1;
	}
	else
		return "Unrecognized character:  '" + expression.substr(0, 1) + "'.";
	lastWasOperator = thisWasOperator;
	return std::string();
}

//=============================================================================
// Class:			ExpressionTree
// Function:		ProcessOperator
//
// Description:		Processes the next operator in the expression, adding it
//					to the appropriate stack.  This method enforces the order
//					of operations.
//
// Input Arguments:
//		operatorStack	= std::stack<std::string>&
//		s				= const std::string& representing the next operator
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
void ExpressionTree::ProcessOperator(std::stack<std::string> &operatorStack,
	const std::string &s)
{
	// Handle operator precedence
	while (!operatorStack.empty())
	{
		if ((!NextIsOperator(operatorStack.top()) ||
			!OperatorShift(operatorStack.top(), s)))
			break;
		PopStackToQueue(operatorStack);
	}
	operatorStack.push(s);
}

//=============================================================================
// Class:			ExpressionTree
// Function:		ProcessCloseParenthese
//
// Description:		Adjusts the stacks in response to encountering a close
//					parenthese.
//
// Input Arguments:
//		operatorStack	= std::stack<std::string>&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
void ExpressionTree::ProcessCloseParenthese(
	std::stack<std::string> &operatorStack)
{
	while (!operatorStack.empty())
	{
		if (operatorStack.top().compare("(") == 0)
			break;
		PopStackToQueue(operatorStack);
	}

	if (operatorStack.empty())
	{
		assert(false);
		// Should never happen due to prior parenthese balance checks
		//return _T("Imbalanced parentheses!");
	}

	operatorStack.pop();
}

//=============================================================================
// Class:			ExpressionTree
// Function:		EvaluateExpression
//
// Description:		Evaluates the expression in the queue using Reverse Polish
//					Notation.
//
// Input Arguments:
//		None
//
// Output Arguments:
//		results	= std::string&
//
// Return Value:
//		std::string containing a description of any errors, or empty string on
//		success
//
//=============================================================================
std::string ExpressionTree::EvaluateExpression(std::string &results)
{
	std::string next, errorString;

	std::stack<double> doubleStack;
	std::stack<std::string> stringStack;
	std::stack<bool> useDoubleStack;

	while (!mOutputQueue.empty())
	{
		next = mOutputQueue.front();
		mOutputQueue.pop();

		if (!EvaluateNext(next, doubleStack, stringStack,
			useDoubleStack, errorString))
			return errorString;
	}

	if (useDoubleStack.size() > 1)
		return "Not enough operators!";

	if (useDoubleStack.top())
	{
		std::ostringstream ss;
		ss.precision(IDMath::GetPrecision(doubleStack.top(), mPrintfPrecision));
		ss << std::fixed << doubleStack.top();
		results = ss.str();
	}
	else
		results = stringStack.top();

	return "";
}

//=============================================================================
// Class:			ExpressionTree
// Function:		PopStackToQueue
//
// Description:		Removes the top entry of the stack and puts it in the queue.
//
// Input Arguments:
//		stack	= std::stack<std::string>& to be popped
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
void ExpressionTree::PopStackToQueue(std::stack<std::string> &stack)
{
	mOutputQueue.push(stack.top());
	stack.pop();
}

//=============================================================================
// Class:			ExpressionTree
// Function:		EmptyStackToQueue
//
// Description:		Empties the contents of the stack into the queue.
//
// Input Arguments:
//		stack	= std::stack<std::string>& to be emptied
//
// Output Arguments:
//		None
//
// Return Value:
//		bool, true for success, false otherwise (imbalance parentheses)
//
//=============================================================================
bool ExpressionTree::EmptyStackToQueue(std::stack<std::string> &stack)
{
	while (!stack.empty())
	{
		if (stack.top().compare("(") == 0)
			return false;
		PopStackToQueue(stack);
	}

	return true;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		NextIsNumber
//
// Description:		Determines if the next portion of the expression is a
//					number. Some cleverness is required to tell the difference
//					between a minus sign and a negative sign (minus sign would
//					return false).
//
// Input Arguments:
//		s				= const std::string& containing the expression
//		lastWasOperator	= const bool& indicating whether or not the last thing
//						  on the stack is an operator
//
// Output Arguments:
//		stop	= unsigned int* (optional) indicating length of number
//
// Return Value:
//		bool, true if a number is next in the expression
//
//=============================================================================
bool ExpressionTree::NextIsNumber(const std::string &s, unsigned int *stop,
	const bool &lastWasOperator)
{
	if (s.length() == 0)
		return false;

	bool foundDecimal = s[0] == '.';
	if (foundDecimal ||
		(int(s[0]) >= int('0') && int(s[0]) <= int('9')) ||
		(s[0] == '-' && lastWasOperator && NextIsNumber(s.substr(1), nullptr, false)))
	{
		unsigned int i;
		for (i = 1; i < s.length(); ++i)
		{
			if (s[i] == '.')
			{
				if (foundDecimal)
					return false;
				foundDecimal = true;
			}
			else if (int(s[i]) < int('0') || int(s[i]) > int('9'))
				break;
		}

		if (stop)
			*stop = i;
		return true;
	}

	return false;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		NextIsS
//
// Description:		Determines if the next portion of the expression is
//					complex frequency (s) or discrete time (z).
//
// Input Arguments:
//		s		= const std::string& containing the expression
//
// Output Arguments:
//		stop	= unsigned int* (optional) indicating length
//
// Return Value:
//		bool, true if a dataset is next in the expression
//
//=============================================================================
bool ExpressionTree::NextIsS(const std::string &s, unsigned int *stop)
{
	if (s[0] == 's' || s[0] == 'z')
	{
		if (s.length() > 1 &&
			((s[1] >= 'a' && s[1] <= 'z') ||
			((s[1] >= 'A' && s[1] <= 'Z'))))
			return false;
		
		if (stop)
			*stop = 1;
		return true;
	}

	return false;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		NextIsOperator
//
// Description:		Determines if the next portion of the expression is an
//					operator.
//
// Input Arguments:
//		s		= const std::string& containing the expression
//
// Output Arguments:
//		stop	= unsigned int* (optional) indicating length of operator
//
// Return Value:
//		bool, true if an operator is next in the expression
//
//=============================================================================
bool ExpressionTree::NextIsOperator(const std::string &s, unsigned int *stop)
{
	if (s.length() == 0)
		return false;

	if (s[0] == '+' ||// From least precedence
		s[0] == '-' ||
		s[0] == '*' ||
		s[0] == '/' ||
		s[0] == '%' ||
		s[0] == '^')// To most precedence
	{
		if (stop)
			*stop = 1;
		return true;
	}

	return false;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		OperatorShift
//
// Description:		Determines if the new operator requires a shift in
//					operator placement.
//
// Input Arguments:
//		stackString	= const std::string& containing the expression
//		newString	= const std::string& containing the expression
//
// Output Arguments:
//		None
//
// Return Value:
//		bool, true if shifting needs to occur
//
//=============================================================================
bool ExpressionTree::OperatorShift(const std::string &stackString,
	const std::string &newString) const
{
	unsigned int stackPrecedence = GetPrecedence(stackString);
	unsigned int newPrecedence = GetPrecedence(newString);

	if (stackPrecedence == 0 || newPrecedence == 0)
		return false;

	if (IsLeftAssociative(newString[0]))
	{
		if (newPrecedence <= stackPrecedence)
			return true;
	}
	else if (newPrecedence < stackPrecedence)
		return true;

	return false;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		GetPrecedence
//
// Description:		Determines the precedence of the specified operator
//					(higher values are performed first)
//
// Input Arguments:
//		s	= const std::string& containing the operator
//
// Output Arguments:
//		None
//
// Return Value:
//		unsigned int representing the precedence
//
//=============================================================================
unsigned int ExpressionTree::GetPrecedence(const std::string &s) const
{
	if (s.length() != 1)
		return 0;

	if (s[0] == '+' ||
		s[0] == '-')
		return 2;
	else if (s[0] == '*' ||
		s[0] == '/' ||
		s[0] == '%')
		return 3;
	else if (s[0] == '^')
		return 4;

	return 0;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		IsLeftAssociative
//
// Description:		Determines if the specified operator is left or right
//					associative.
//
// Input Arguments:
//		c	= const char&
//
// Output Arguments:
//		None
//
// Return Value:
//		bool, true if left associative
//
//=============================================================================
bool ExpressionTree::IsLeftAssociative(const char &c) const
{
	switch (c)
	{
	case '^':
		return false;

	default:
		return true;
	}
}

//=============================================================================
// Class:			ExpressionTree
// Function:		PushToStack
//
// Description:		Pushes the specified value onto the stack.
//
// Input Arguments:
//		value			= const double&
//		doubleStack		= std::stack<double>&
//		useDoubleStack	= std::stack<bool>&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
void ExpressionTree::PushToStack(const double &value,
	std::stack<double> &doubleStack, std::stack<bool> &useDoubleStack) const
{
	doubleStack.push(value);
	useDoubleStack.push(true);
}

//=============================================================================
// Class:			ExpressionTree
// Function:		PushToStack
//
// Description:		Pushes the specified dataset onto the stack.
//
// Input Arguments:
//		s				= const std::string&
//		stringStack		= std::stack<std::string>&
//		useDoubleStack	= std::stack<bool>&
//
// Output Arguments:
//		None
//
// Return Value:
//		None
//
//=============================================================================
void ExpressionTree::PushToStack(const std::string &s,
	std::stack<std::string> &stringStack, std::stack<bool> &useDoubleStack) const
{
	stringStack.push(s);
	useDoubleStack.push(false);
}

//=============================================================================
// Class:			ExpressionTree
// Function:		PopFromStack
//
// Description:		Pops the next value from the top of the appropriate stack.
//
// Input Arguments:
//		doubleStack		= std::stack<double>&
//		stringStack		= std::stack<std::string>&
//		useDoubleStack	= std::stack<bool>&
//
// Output Arguments:
//		string			= std::string&
//		value			= double&
//
// Return Value:
//		bool, true if a double was popped, false otherwise
//
//=============================================================================
bool ExpressionTree::PopFromStack(std::stack<double> &doubleStack,
	std::stack<std::string> &stringStack, std::stack<bool> &useDoubleStack,
	std::string& string, double &value) const
{
	bool useDouble = useDoubleStack.top();
	useDoubleStack.pop();

	if (useDouble)
	{
		assert(!doubleStack.empty());
		value = doubleStack.top();
		doubleStack.pop();
	}
	else
	{
		assert(!stringStack.empty());
		string = stringStack.top();
		stringStack.pop();
	}

	return useDouble;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		ApplyOperation
//
// Description:		Applies the specified operation to the specified operands.
//
// Input Arguments:
//		operation	= const std::string& describing the function to apply
//		first		= const double&
//		second		= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		double containing the result of the operation
//
//=============================================================================
double ExpressionTree::ApplyOperation(const std::string &operation,
	const double &first, const double &second) const
{
	if (operation.compare("+") == 0)
		return second + first;
	else if (operation.compare("-") == 0)
		return second - first;
	else if (operation.compare("*") == 0)
		return second * first;
	else if (operation.compare("/") == 0)
		return second / first;
	else if (operation.compare("%") == 0)
		return fmod(second, first);
	else if (operation.compare("^") == 0)
		return pow(second, first);

	assert(false);
	return 0.0;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		ApplyOperation
//
// Description:		Applies the specified operation to the specified operands.
//
// Input Arguments:
//		operation	= const std::string& describing the function to apply
//		first		= const std::string&
//		second		= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string containing the result of the operation
//
//=============================================================================
std::string ExpressionTree::ApplyOperation(const std::string &operation,
	const std::string &first, const std::string &second) const
{
	/*if (operation.compare("+") == 0)
		return second + first;
	else if (operation.compare("-") == 0)
		return second - first;
	else */if (operation.compare("*") == 0)
		return StringMultiply(first, second);
	return second + operation + first;

/*	assert(false);
	return "";*/
}

//=============================================================================
// Class:			ExpressionTree
// Function:		ApplyOperation
//
// Description:		Applies the specified operation to the specified operands.
//
// Input Arguments:
//		operation	= const std::string& describing the function to apply
//		first		= const std::string&
//		second		= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string containing the result of the operation
//
//=============================================================================
std::string ExpressionTree::ApplyOperation(const std::string &operation,
	const std::string &first, const double &second) const
{
	if (operation.compare("+") == 0)
		return StringAdd(first, second);
	else if (operation.compare("-") == 0)
		return StringSubtract(first, second);
	else if (operation.compare("*") == 0)
		return StringMultiply(first, second);

	assert(false);
	return "";
}

//=============================================================================
// Class:			ExpressionTree
// Function:		ApplyOperation
//
// Description:		Applies the specified operation to the specified operands.
//
// Input Arguments:
//		operation	= const std::string& describing the function to apply
//		first		= const double&
//		second		= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string containing the result of the operation
//
//=============================================================================
std::string ExpressionTree::ApplyOperation(const std::string &operation,
	const double &first, const std::string &second) const
{
	if (operation.compare("+") == 0)
		return StringAdd(first, second);
	else if (operation.compare("-") == 0)
		return StringSubtract(first, second);
	else if (operation.compare("*") == 0)
		return StringMultiply(first, second);
	else if (operation.compare("/") == 0)
		return StringDivide(first, second);
	else if (operation.compare("^") == 0)
		return StringPower(first, second);

	assert(false);
	return "";
}

//=============================================================================
// Class:			ExpressionTree
// Function:		EvaluateOperator
//
// Description:		Evaluates the operator specified.
//
// Input Arguments:
//		operator		= const std::string& describing the function to apply
//		doubleStack		= std::stack<double>&
//		stringStack		= std::stack<std::string>&
//		useDoubleStack	= std::stack<bool>&
//
// Output Arguments:
//		errorString		= std::string&
//
// Return Value:
//		bool, true for success, false otherwise
//
//=============================================================================
bool ExpressionTree::EvaluateOperator(const std::string &operation, std::stack<double> &doubleStack,
	std::stack<std::string> &stringStack, std::stack<bool> &useDoubleStack, std::string &errorString) const
{
	double value1, value2;
	std::string string1, string2;

	if (useDoubleStack.size() < 2)
		return EvaluateUnaryOperator(operation, doubleStack, stringStack, useDoubleStack, errorString);
	else if (PopFromStack(doubleStack, stringStack, useDoubleStack, string1, value1))
	{
		if (PopFromStack(doubleStack, stringStack, useDoubleStack, string2, value2))
			PushToStack(ApplyOperation(operation, value1, value2), doubleStack, useDoubleStack);
		else
			PushToStack(ApplyOperation(operation, value1, string2), stringStack, useDoubleStack);
	}
	else if (PopFromStack(doubleStack, stringStack, useDoubleStack, string2, value2))
		PushToStack(ApplyOperation(operation, string1, value2), stringStack, useDoubleStack);
	else
		PushToStack(ApplyOperation(operation, string1, string2), stringStack, useDoubleStack);

	return true;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		EvaluateUnaryOperator
//
// Description:		Evaluates the operator specified.  The only unary operator
//					we recognize is minus (negation).
//
// Input Arguments:
//		operator		= const std::string& describing the function to apply
//		doubleStack		= std::stack<double>&
//		stringStack		= std::stack<std::string>&
//		useDoubleStack	= std::stack<bool>&
//
// Output Arguments:
//		errorString		= std::string&
//
// Return Value:
//		bool, true for success, false otherwise
//
//=============================================================================
bool ExpressionTree::EvaluateUnaryOperator(const std::string &operation, std::stack<double> &doubleStack,
	std::stack<std::string> &stringStack, std::stack<bool> &useDoubleStack, std::string &errorString) const
{
	if (operation.compare("-") != 0)
	{
		errorString = "Attempting to apply operator without two operands!";
		return false;
	}

	double value;
	std::string string;
	if (PopFromStack(doubleStack, stringStack, useDoubleStack, string, value))
		PushToStack(ApplyOperation("*", -1.0, value), doubleStack, useDoubleStack);
	else
		PushToStack(ApplyOperation("*", -1.0, string), stringStack, useDoubleStack);

	return true;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		EvaluateNumber
//
// Description:		Evaluates the number specified.
//
// Input Arguments:
//		number			= const std::string& describing the function to apply
//		doubleStack		= std::stack<double>&
//		setStack		= std::stack<Dataset2D>&
//		useDoubleStack	= std::stack<bool>&
//
// Output Arguments:
//		errorString		= std::string&
//
// Return Value:
//		bool, true for success, false otherwise
//
//=============================================================================
bool ExpressionTree::EvaluateNumber(const std::string &number, std::stack<double> &doubleStack,
	std::stack<bool> &useDoubleStack, std::string &errorString) const
{
	double value;

	std::istringstream ss(number);
	if ((ss >> value).fail())
	{
		std::ostringstream oss;
		oss << "Could not convert " << number << " to a number.";
		errorString = oss.str();
		return false;
	}

	PushToStack(value, doubleStack, useDoubleStack);

	return true;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		EvaluateNext
//
// Description:		Determines how to evaluate the specified term and takes
//					appropriate action.
//
// Input Arguments:
//		next			= const std::string&
//		doubleStack		= std::stack<double>&
//		stringStack		= std::stack<std::string>&
//		useDoubleStack	= std::stack<bool>&
//
// Output Arguments:
//		errorString		= std::string&
//
// Return Value:
//		bool, true for valid operation, false otherwise
//
//=============================================================================
bool ExpressionTree::EvaluateNext(const std::string &next, std::stack<double> &doubleStack,
		std::stack<std::string> &stringStack, std::stack<bool> &useDoubleStack, std::string &errorString) const
{
	if (NextIsNumber(next))
		return EvaluateNumber(next, doubleStack, useDoubleStack, errorString);
	else if(NextIsOperator(next))
		return EvaluateOperator(next, doubleStack, stringStack, useDoubleStack, errorString);
	else if (NextIsS(next))
	{
		PushToStack(next, stringStack, useDoubleStack);
		return true;
	}
	else
		errorString = "Unable to evaluate '" + next + "'.";

	return false;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		StringAdd
//
// Description:		Performs arithmatic on the arguments, returns a string.
//
// Input Arguments:
//		first	= const std::string&
//		second	= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string ExpressionTree::StringAdd(const std::string &first, const double &second) const
{
	std::ostringstream ss;
	ss.precision(IDMath::GetPrecision(second, mPrintfPrecision));
	ss << std::fixed << second << '+' << first;
	return ss.str();
}

//=============================================================================
// Class:			ExpressionTree
// Function:		StringAdd
//
// Description:		Performs arithmatic on the arguments, returns a string.
//
// Input Arguments:
//		first	= const double&
//		second	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string ExpressionTree::StringAdd(const double &first, const std::string &second) const
{
	std::ostringstream ss;
	ss.precision(IDMath::GetPrecision(first, mPrintfPrecision));
	ss << second << '+' << std::fixed << first;
	return ss.str();
}

//=============================================================================
// Class:			ExpressionTree
// Function:		StringSubtract
//
// Description:		Performs arithmatic on the arguments, returns a string.
//
// Input Arguments:
//		first	= const std::string&
//		second	= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string ExpressionTree::StringSubtract(const std::string &first, const double &second) const
{
	std::ostringstream ss;
	ss.precision(IDMath::GetPrecision(second, mPrintfPrecision));
	ss << std::fixed << second << '-' << first;
	return ss.str();
}

//=============================================================================
// Class:			ExpressionTree
// Function:		StringSubtract
//
// Description:		Performs arithmatic on the arguments, returns a string.
//
// Input Arguments:
//		first	= const double&
//		second	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string ExpressionTree::StringSubtract(const double &first, const std::string &second) const
{
	std::ostringstream ss;
	ss.precision(IDMath::GetPrecision(first, mPrintfPrecision));
	ss << second << '-' << std::fixed << first;
	return ss.str();
}

//=============================================================================
// Class:			ExpressionTree
// Function:		StringMultiply
//
// Description:		Performs arithmatic on the arguments, returns a string.
//
// Input Arguments:
//		first	= const std::string&
//		second	= const double&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string ExpressionTree::StringMultiply(const std::string &first,
	const double &second) const
{
	std::vector<std::pair<int, double>> terms(FindPowersAndCoefficients(
		BreakApartTerms(first)));
	std::string expression;
	for (const auto& term : terms)
		AddToExpressionString(expression, term.second * second, term.first);

	return expression;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		StringMultiply
//
// Description:		Performs arithmatic on the arguments, returns a string.
//
// Input Arguments:
//		first	= const std::string&
//		second	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string ExpressionTree::StringMultiply(const std::string &first,
	const std::string &second) const
{
	std::vector<std::pair<int, double>> firstTerms(
		FindPowersAndCoefficients(BreakApartTerms(first)));
	std::vector<std::pair<int, double>> secondTerms(
		FindPowersAndCoefficients(BreakApartTerms(second)));
	std::vector<std::pair<int, double>> terms;
	for (const auto& firstTerm : firstTerms)
	{
		for (const auto& secondTerm : secondTerms)
			terms.push_back(std::pair<int, double>(
				firstTerm.first + secondTerm.first,
				firstTerm.second * secondTerm.second));
	}

	std::string expression;
	for (const auto& term : terms)
		AddToExpressionString(expression, term.second, term.first);

	return expression;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		StringMultiply
//
// Description:		Performs arithmatic on the arguments, returns a string.
//
// Input Arguments:
//		first	= const double&
//		second	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string ExpressionTree::StringMultiply(const double &first, const std::string &second) const
{
	return StringMultiply(second, first);
}

//=============================================================================
// Class:			ExpressionTree
// Function:		StringDivide
//
// Description:		Performs arithmatic on the arguments, returns a string.
//
// Input Arguments:
//		first	= const double&
//		second	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string ExpressionTree::StringDivide(const double &first, const std::string &second) const
{
	return StringMultiply(second, 1.0 / first);
}

//=============================================================================
// Class:			ExpressionTree
// Function:		StringPower
//
// Description:		Performs arithmatic on the arguments, returns a string.
//					For positive powers, expand and do the multiplication.
//					For negative powers (i.e. z-domain stuff), add them to the
//					string.  Assumes exponent is an integer.
//
// Input Arguments:
//		first	= const double&
//		second	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::string
//
//=============================================================================
std::string ExpressionTree::StringPower(const double &first, const std::string &second) const
{
	if (first < 0.0)
	{
		std::ostringstream ss;
		ss << second << '^' << static_cast<int>(first);
		return ss.str();
	}

	std::string result(second);
	unsigned int i;
	for (i = 1; i < static_cast<unsigned int>(first); ++i)
		result = StringMultiply(result, second);

	return result;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		BreakApartTerms
//
// Description:		Breaks apart all the terms in the string expression.  Be
//					wary of negative signs preceded by another operator!
//
// Input Arguments:
//		s	= const std::string&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::vector<std::string>
//
//=============================================================================
std::vector<std::string> ExpressionTree::BreakApartTerms(const std::string &s)
{
	std::vector<std::string> terms;
	std::string::size_type start(0), end(0);
	while (end != std::string::npos)
	{
		end = FindEndOfNextTerm(s, start);

		if (start > 0 && s.substr(start - 1, 1).compare("-") == 0)
		{
			if (end != std::string::npos)
				terms.push_back(s.substr(start - 1, end + 1));
			else
				terms.push_back(s.substr(start - 1));
		}
		else
			terms.push_back(s.substr(start, end));

		start += end + 1;
	}

	return terms;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		FindEndOfNextTerm
//
// Description:		Finds the end of the next term in the string.
//
// Input Arguments:
//		s		= const std::string&
//		start	= const unsigned int&
//
// Output Arguments:
//		None
//
// Return Value:
//		unsigned int
//
//=============================================================================
unsigned int ExpressionTree::FindEndOfNextTerm(const std::string &s, const unsigned int &start)
{
	std::string::size_type end, plusEnd, minusEnd;

	plusEnd = s.substr(start).find('+', 1);
	minusEnd = s.substr(start).find('-', 1);

	if (minusEnd < plusEnd && start + minusEnd > 0 && NextIsOperator(s.substr(start + minusEnd - 1, 1)))
	{
		auto nextMinus(s.substr(start + minusEnd + 1).find('-'));
		if (nextMinus != std::string::npos)
			minusEnd += nextMinus + 1;
		else
			minusEnd = nextMinus;
	}
	end = std::min(plusEnd, minusEnd);

	if (end != std::string::npos && NextIsOperator(s.substr(start + end - 1)))
	{
		plusEnd = s.substr(start + end).find('+');
		minusEnd = s.substr(start + end).find('-');
		end += std::min(plusEnd, minusEnd);
	}

	return end;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		FindPowersAndCoefficients
//
// Description:		Breaks a (previously separated) set of terms into a coefficient
//					and a power of the algebraic variable.
//
// Input Arguments:
//		terms	= const std::vector<std::string>&
//
// Output Arguments:
//		None
//
// Return Value:
//		std::vector<std::pair<int, double>>
//
//=============================================================================
std::vector<std::pair<int, double>> ExpressionTree::FindPowersAndCoefficients(const std::vector<std::string> &terms)
{
	std::vector<std::pair<int, double>> processedTerms;
	double temp;
	for (const auto& term : terms)
	{
		int count(0);
		std::string::size_type start(0), end(0);
		double coefficient(1.0);
		while (end != std::string::npos)
		{
			end = term.substr(start).find('*');
			std::istringstream ss(term.substr(start, end));
			if (!(ss >> temp).fail())
				coefficient = temp;
			else
			{
				if (term[0] == '-' && coefficient == 1.0)
				{
					coefficient = -1.0;
					++start;
					if (end != std::string::npos)
						++end;
				}

				count += GetTermPower(term.substr(start), start, end);
			}
			start += end + 1;
		}

		processedTerms.push_back(std::pair<int, double>(count, coefficient));
	}

	return processedTerms;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		GetTermPower
//
// Description:		Returns the value of the power for the specified term
//					(power of s or z).
//
// Input Arguments:
//		s		= const std::string&
//
// Output Arguments:
//		start	= std::string::size_type&
//		end		= std::string::size_type&
//
// Return Value:
//		int
//
//=============================================================================
int ExpressionTree::GetTermPower(const std::string &s,
	std::string::size_type &start, std::string::size_type &end)
{
	std::string::size_type power;
	if (s[0] == 's' || s[0] == 'z')
	{
		power = s.find('^');
		if (power == std::string::npos)
			return 1;
		else
		{
			start += power + 1;
			end = s.find('*');
			std::istringstream ss(s.substr(power + 1, end));
			int powerInt;
			if (!(ss >> powerInt).fail())
				return powerInt;
		}
	}

	return 0;
}

//=============================================================================
// Class:			ExpressionTree
// Function:		AddToExpressionString
//
// Description:		Adds the next term to the string.  Handles signed terms,
//					cleans up for terms with coefficient == 1.0, etc.
//
// Input Arguments:
//		coefficient	= const double&
//		power		= const int&
//
// Output Arguments:
//		expression	= std::string&
//
// Return Value:
//		None
//
//=============================================================================
void ExpressionTree::AddToExpressionString(std::string &expression,
	const double &coefficient, const int &power) const
{
	if (!expression.empty())
	{
		if (coefficient > 0.0)
			expression.append("+");
		else
			expression.append("-");
	}

	if (coefficient == 1.0 && power != 0)
	{
		if (power == 1)
			expression.append("s");
		else
		{
			std::ostringstream ss;
			ss << "s^" << power;
			expression.append(ss.str());
		}
	}
	else
	{
		std::ostringstream ss;
		ss.precision(IDMath::GetPrecision(coefficient, mPrintfPrecision));
		if (expression.empty())
			ss << std::fixed << coefficient;
		else
			ss << std::fixed << fabs(coefficient);

		if (power == 1)
			ss << "*s";
		else if (power != 0)
			ss << "*s^" << power;
		expression.append(ss.str());
	}
}
