#include <iostream>

#include "snlUnitTest.hpp"

using namespace std;

snlUnitTest::~snlUnitTest()
{

}

snlUnitTest::snlUnitTest(const char* unitTestName) : _name{unitTestName}, _result{true}, _level{0}, _numChildUnitTests{0}
{
}

bool snlUnitTest::run()
{
	cout << _indentTabs << "Running:" << _name << "\n";

	// Run main tests first.
	runTests();
	// Post test run hook.
	postRunTests();

	// Run child tests after main tests so as to give main tests a chance to test for more broader conditions.
	bool childResult;
	for(int index = 0; index < _numChildUnitTests; index++)
	{
		childResult = _childUnitTests[index] -> run();
		_result = _result && childResult;
		// Stop on first failure.
		if(!childResult) break;
	}

	// Output the overall result.
	if(_result)
	{
		cout <<  _indentTabs << _name << " : " << PASS_STRING << "\n";
	}
	else
	{
		cout <<  _indentTabs << _name << " : " << FAIL_STRING << "\n";
	}

	return _result;
}

bool snlUnitTest::getResult()
{
	return _result;
}

void snlUnitTest::addChildUnitTest(snlUnitTest* childUnitTest)
{
	_childUnitTests[_numChildUnitTests++] = childUnitTest;
	childUnitTest -> setLevel(_level + 1);
}

void snlUnitTest::notifyTestResult(const char* testName, bool result, const char* resultMessage)
{
	_result = _result && result;

	// Tab the output in level in from this unit tests tab indent.
	cout << "\t" << _indentTabs;

	// Output the result.
	if(result)
	{
		cout << testName << " : " << PASS_STRING  << " : " << resultMessage << "\n";
	}
	else
	{
		cout << testName << " : " << FAIL_STRING  << " : " << resultMessage << "\n";
	}
}

void snlUnitTest::setLevel(int level)
{
	_level = level;

	// Pre-create tab output according to the nesting level.
	_indentTabs = "";
	for(int level = 0; level < _level; level++)
	{
		_indentTabs += "\t";
	}
}
