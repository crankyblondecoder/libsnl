#include "snlVectorUnitTest.hpp"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	cout << "\n" << "*** Starting libSNL Unit Tests ***" << "\n\n";

	snlVectorUnitTest vectorUnitTest;
	vectorUnitTest.run();

	cout << "\n" << "*** libSNL Unit Tests complete ***" << "\n";
}