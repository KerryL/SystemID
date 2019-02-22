// File:  expressionTree.h
// Date:  2/8/2019
// Auth:  K. Loux
// Desc:  Handles user-specified mathematical operations on datasets.

#ifndef EXPRESSION_TREE_H_
#define EXPRESSION_TREE_H_

// Standard C++ headers
#include <queue>
#include <stack>
#include <string>
#include <memory>

/// Class for processing user-specified mathematical operations.  Uses a
/// shunting yard algorithm to build and evaluate expression trees from
/// user-specified strings.
class ExpressionTree
{
public:
	/// Solves the specified expression by simplifying and combining like
	/// terms.
	///
	/// \param expression             Expression to evaluate.
	/// \param solvedExpression [out] The result of the simplified expression.
	///
	/// \returns A description of any parsing/evaluation errors, or an empty
	///          string for success.
	std::string Solve(std::string expression, std::string &solvedExpression);

	/// Breaks the specified expression string into separate terms.
	///
	/// \param s Expression string.
	///
	/// \returns A list of separate terms.
	static std::vector<std::string> BreakApartTerms(const std::string &s);

	/// Processes each term to extract the value of the coefficient and the
	/// power to which the variable is raised.
	///
	/// \param terms List of terms to process.
	///
	/// \returns Values that describe each term.
	static std::vector<std::pair<int, double>> FindPowersAndCoefficients(
		const std::vector<std::string> &terms);

private:
	static const unsigned int mPrintfPrecision;

	double mXAxisFactor;

	std::queue<std::string> mOutputQueue;

	std::string ParseExpression(const std::string &expression);
	std::string ParseNext(const std::string &expression, bool &lastWasOperator,
		unsigned int &advance, std::stack<std::string> &operatorStack);
	std::string EvaluateExpression(std::string &results);

	void ProcessOperator(std::stack<std::string> &operatorStack, const std::string &s);
	void ProcessCloseParenthese(std::stack<std::string> &operatorStack);

	static bool NextIsNumber(const std::string &s, unsigned int *stop = nullptr, const bool &lastWasOperator = true);
	static bool NextIsOperator(const std::string &s, unsigned int *stop = nullptr);
	static bool NextIsS(const std::string &s, unsigned int *stop = nullptr);

	static unsigned int FindEndOfNextTerm(const std::string &s, const unsigned int &start);
	static int GetTermPower(const std::string &s, std::string::size_type &start, std::string::size_type &end);

	bool IsLeftAssociative(const char &c) const;
	bool OperatorShift(const std::string &stackString, const std::string &newString) const;

	void PopStackToQueue(std::stack<std::string> &stack);
	bool EmptyStackToQueue(std::stack<std::string> &stack);
	unsigned int GetPrecedence(const std::string &s) const;

	void PushToStack(const double &value, std::stack<double> &doubleStack,
		std::stack<bool> &useDoubleStack) const;

	double ApplyOperation(const std::string &operation, const double &first, const double &second) const;

	bool EvaluateNumber(const std::string &number, std::stack<double> &doubleStack,
		std::stack<bool> &useDoubleStack, std::string &errorString) const;

	void PushToStack(const std::string &s, std::stack<std::string> &stringStack,
		std::stack<bool> &useDoubleStack) const;
	bool PopFromStack(std::stack<double> &doubleStack, std::stack<std::string> &stringStack,
		std::stack<bool> &useDoubleStack, std::string& string, double &value) const;

	bool EvaluateNext(const std::string &next, std::stack<double> &doubleStack,
		std::stack<std::string> &stringStack, std::stack<bool> &useDoubleStack, std::string &errorString) const;
	bool EvaluateOperator(const std::string &operation, std::stack<double> &doubleStack,
		std::stack<std::string> &stringStack, std::stack<bool> &useDoubleStack, std::string &errorString) const;
	bool EvaluateUnaryOperator(const std::string &operation, std::stack<double> &doubleStack,
		std::stack<std::string> &stringStack, std::stack<bool> &useDoubleStack, std::string &errorString) const;

	std::string ApplyOperation(const std::string &operation, const std::string &first, const std::string &second) const;
	std::string ApplyOperation(const std::string &operation, const std::string &first, const double &second) const;
	std::string ApplyOperation(const std::string &operation, const double &first, const std::string &second) const;

	bool ParenthesesBalanced(const std::string &expression) const;

	std::string StringAdd(const std::string &first, const double &second) const;
	std::string StringAdd(const double &first, const std::string &second) const;

	std::string StringSubtract(const std::string &first, const double &second) const;
	std::string StringSubtract(const double &first, const std::string &second) const;

	std::string StringMultiply(const std::string &first, const double &second) const;
	std::string StringMultiply(const std::string &first, const std::string &second) const;
	std::string StringMultiply(const double &first, const std::string &second) const;

	std::string StringDivide(const double &first, const std::string &second) const;

	std::string StringPower(const double &first, const std::string &second) const;

	void AddToExpressionString(std::string &expression, const double &coefficient, const int &power) const;
};

#endif// EXPRESSION_TREE_H_
