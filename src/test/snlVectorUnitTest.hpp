#include "snlUnitTest.hpp"

#include "../snlVector.h"

class snlVectorUnitTest : public snlUnitTest
{
	public:

		snlVectorUnitTest() : snlUnitTest("snlVector Test") {}

	protected:

		virtual void runTests()
		{
			testConstruction();
		}

	private:

		void testConstruction()
		{
			bool allPassed = true;

			{
				snlVector vec;
				if(vec.elements[0] != 0.0 || vec.elements[1] != 0.0 || vec.elements[2] != 0.0 || vec.elements[3] != 0.0)
				{
					notifyTestResult("Default constructor", false, "Incorrect initialisation.");
					allPassed = false;
				}
			}

			{
				snlVector vec(1, 2, 3, 4);
				if(vec.elements[0] != 1 || vec.elements[1] != 2 || vec.elements[2] != 3 || vec.elements[3] != 4)
				{
					notifyTestResult("Element init constructor", false, "Incorrect initialisation.");
					allPassed = false;
				}
			}

			// Catch all for all passed.
			if(allPassed) notifyTestResult("Construction tests", true, "All construction tests passed.");
		}
};