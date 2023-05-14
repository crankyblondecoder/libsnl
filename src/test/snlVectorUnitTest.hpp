#include <cmath>

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
				if(vec.components[0] != 0.0 || vec.components[1] != 0.0 || vec.components[2] != 0.0 || vec.components[3] != 0.0)
				{
					notifyTestResult("Default constructor", false, "Incorrect initialisation.");
					allPassed = false;
				}
			}

			{
				snlVector vec(1, 2, 3, 4);
				if(vec.components[0] != 1 || vec.components[1] != 2 || vec.components[2] != 3 || vec.components[3] != 4)
				{
					notifyTestResult("Element init constructor", false, "Incorrect initialisation.");
					allPassed = false;
				}
			}

			{
				snlVector from(1, 2, 3, 4);
				snlVector vec;
				vec.assign(from);
				if(vec.components[0] != 1 || vec.components[1] != 2 || vec.components[2] != 3 || vec.components[3] != 4)
				{
					notifyTestResult("Assign", false, "Incorrect assignment.");
					allPassed = false;
				}
			}

			{
				snlVector vec2(2, 3, 4, 5);
				snlVector vec1(1, 2, 3, 4);
				snlVector vec;
				vec.diff(vec1, vec2);
				if(vec.components[0] != 1 || vec.components[1] != 1 || vec.components[2] != 1 || vec.components[3] != 1)
				{
					notifyTestResult("Diff", false, "Incorrect difference.");
					allPassed = false;
				}
			}

			{
				snlVector vec2(2, 3, 4, 5);
				snlVector vec1(1, 2, 3, 4);
				snlVector vec;
				vec.crossProduct(vec1, vec2);
				if(vec.components[0] != -1 || vec.components[1] != 2 || vec.components[2] != -1 || vec.components[3] != 0)
				{
					notifyTestResult("CrossProduct", false, "Incorrect result.");
					allPassed = false;
				}
			}

			{
				snlVector vec2(2, 3, 4, 5);
				snlVector vec1(1, 2, 3, 4);
				if(vec1.dot(vec2) != 40)
				{
					notifyTestResult("Dot", false, "Incorrect result.");
					allPassed = false;
				}
			}

			{
				snlVector vec(1, 2, 3, 4);
				if(vec.lengthSqrd() != 30)
				{
					notifyTestResult("LengthSqrd", false, "Incorrect result.");
					allPassed = false;
				}
			}

			{
				snlVector vec(2, 4, 8, 16);
				if(vec.length() != 18.439088914585774)
				{
					notifyTestResult("Length", false, "Incorrect result.");
					allPassed = false;
				}
			}

			{
				double deg_30 = M_PI / 6.0;
				double deg_60 = M_PI / 3.0;

				snlVector vec1(cos(deg_30), sin(deg_30), 0, 0);
				snlVector vec2(cos(deg_60), sin(deg_60), 0, 0);

				if(vec1.calcAbsCos(vec2) != deg_30)
				{
					notifyTestResult("CalcAbsCos", false, "Incorrect result.");
					allPassed = false;
				}
			}

			// Catch all for all passed.
			if(allPassed) notifyTestResult("Construction tests", true, "All construction tests passed.");
		}
};