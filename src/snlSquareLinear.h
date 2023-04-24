// *** Class for the solution of square linear algebraic equations ***

#ifndef SNL_SQUARELINEAR_H
#define SNL_SQUARELINEAR_H

#include "snlPoint.h"

#ifdef SGI_MIPS

    #include <iostream.h>
    #include <math.h>

#else

    #include <iostream>
    #include <cmath>

    using namespace std;

#endif

class snlSquareLinear
{
    public:

        snlSquareLinear();
        ~snlSquareLinear();
        snlSquareLinear ( int numUnknowns, int numRightHandSides, double* coefficients, double* rightHandSides );

        void solve();

        void print();
        void printCoeffs();
        void printRhSides();

    private:

        int        num_unknowns; // Number of unknowns corresponding to coefficient matrix.
        int        num_rhSides;  // Number of right hand sides.

        double*    coeffs;  // Coefficient matrix array.
        double*    rhSides;  // Array of right hands sides.
};

#endif
