#include "snlMatrix_4x4.h"

snlMatrix_4X4::snlMatrix_4X4()
{
}

snlMatrix_4X4::snlMatrix_4X4()
{
	_elements = new double[16];
	_scratch = new double[16];
}

snlMatrix_4X4::~snlMatrix_4X4()
{
	delete[] _elements;
	delete[] _scratch;
}

snlMatrix_4X4::snlMatrix_4X4(snlMatrix_4X4& copyFrom)
{
	_elements = new double[16];
	_scratch = new double[16];

	_elements[0] = copyFrom._elements[0];
	_elements[1] = copyFrom._elements[1];
	_elements[2] = copyFrom._elements[2];
	_elements[3] = copyFrom._elements[3];

	_elements[4] = copyFrom._elements[4];
	_elements[5] = copyFrom._elements[5];
	_elements[6] = copyFrom._elements[6];
	_elements[7] = copyFrom._elements[7];

	_elements[8] = copyFrom._elements[8];
	_elements[9] = copyFrom._elements[9];
	_elements[10] = copyFrom._elements[10];
	_elements[11] = copyFrom._elements[11];

	_elements[12] = copyFrom._elements[12];
	_elements[13] = copyFrom._elements[13];
	_elements[14] = copyFrom._elements[14];
	_elements[15] = copyFrom._elements[15];
}

snlMatrix_4X4& snlMatrix_4X4::operator = (snlMatrix_4X4& copyFrom)
{
	_elements[0] = copyFrom._elements[0];
	_elements[1] = copyFrom._elements[1];
	_elements[2] = copyFrom._elements[2];
	_elements[3] = copyFrom._elements[3];

	_elements[4] = copyFrom._elements[4];
	_elements[5] = copyFrom._elements[5];
	_elements[6] = copyFrom._elements[6];
	_elements[7] = copyFrom._elements[7];

	_elements[8] = copyFrom._elements[8];
	_elements[9] = copyFrom._elements[9];
	_elements[10] = copyFrom._elements[10];
	_elements[11] = copyFrom._elements[11];

	_elements[12] = copyFrom._elements[12];
	_elements[13] = copyFrom._elements[13];
	_elements[14] = copyFrom._elements[14];
	_elements[15] = copyFrom._elements[15];

	return *this;
}

void snlMatrix_4X4::preMultiply(snlMatrix_4X4& multMatrix)
{
	double cVal;

	double* matrixA;
	double* matrixB;

	// Pre multiply, ie AB, and store in B.
	matrixA = multMatrix._elements;
	matrixB = _elements;

	// The current index of the result being calculated. ie Index into the elements of B.
	int resultIndex = 0;

	// The current index into the B matrix.
	int bIndex = 0;

	// The current B column element values. These will be involved in four calculations each.
	double bVal0;
	double bVal1;
	double bVal2;
	double bVal3;

	// Four passes to fill four columns. Each pass iterates over the entire A matrix, each element once.
	for(int pass = 0; pass < 4; pass ++)
	{
		// Get an entire column of B at a time. This allows the elements of B to be replaced immediately.
		bVal0 = matrixB[bIndex++];
		bVal1 = matrixB[bIndex++];
		bVal2 = matrixB[bIndex++];
		bVal3 = matrixB[bIndex++];

		// Calculate an entire column of the result. Store in B.
		matrixB[resultIndex++] = matrixA[0] * bVal0 + matrixA[4] * bVal1 + matrixA[8] * bVal2 + matrixA[12] * bVal3;
		matrixB[resultIndex++] = matrixA[1] * bVal0 + matrixA[5] * bVal1 + matrixA[9] * bVal2 + matrixA[13] * bVal3;
		matrixB[resultIndex++] = matrixA[2] * bVal0 + matrixA[6] * bVal1 + matrixA[10] * bVal2 + matrixA[14] * bVal3;
		matrixB[resultIndex++] = matrixA[3] * bVal0 + matrixA[7] * bVal1 + matrixA[11] * bVal2 + matrixA[15] * bVal3;
	}
}

void snlMatrix_4X4::multiply(snlVector*, unsigned numVectors)
{
	for(let index = 0; index < numVectors; index++)
	{

	}
	// TODO
	blah;
}

void snlMatrix_4X4::print()
{
	for(int row = 0; row < 4; row ++)
	{
		for(int col = 0; col < 4; col ++)
		{
			cout << _elements [ row +(col * 4) ] << "  ";
		}

		cout << "\n";
	}
}

