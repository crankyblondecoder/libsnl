#include "snlMatrix_4x4.h"

snlMatrix_4X4::snlMatrix_4X4()
{
	ident();
}

snlMatrix_4X4::snlMatrix_4X4(double initialElements[16])
{
	for(int index = 0; index < 16; index++) {

		_elements[index] = initialElements[index];
	}
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
		// Note: The column major ordering makes this quite simple because the result can be immediately stored directly
		//       in this matrix.
		matrixB[resultIndex++] = matrixA[0] * bVal0 + matrixA[4] * bVal1 + matrixA[8] * bVal2 + matrixA[12] * bVal3;
		matrixB[resultIndex++] = matrixA[1] * bVal0 + matrixA[5] * bVal1 + matrixA[9] * bVal2 + matrixA[13] * bVal3;
		matrixB[resultIndex++] = matrixA[2] * bVal0 + matrixA[6] * bVal1 + matrixA[10] * bVal2 + matrixA[14] * bVal3;
		matrixB[resultIndex++] = matrixA[3] * bVal0 + matrixA[7] * bVal1 + matrixA[11] * bVal2 + matrixA[15] * bVal3;
	}
}

void snlMatrix_4X4::multiply(snlVector& vector)
{
	// Simply converts to array version.
	multiply(&vector, 1);
}

void snlMatrix_4X4::multiply(snlVector* vectors, unsigned numVectors)
{
	// Cache of vector elements
	double vElem0;
	double vElem1;
	double vElem2;
	double vElem3;

	// Matrix element index.
	unsigned mElemIndex;

	// Current vector elements
	double* vElems;

	for(unsigned vecIndex = 0; vecIndex < numVectors; vecIndex++)
	{
		vElems = vectors[vecIndex].elements;

		// Vector elements are progressively written over so cache them ahead of time.
		vElem0 = vElems[0];
		vElem1 = vElems[1];
		vElem2 = vElems[2];
		vElem3 = vElems[3];

		// Start at the beginning of the first matrix column.
		mElemIndex = 0;

		vElems[0] = _elements[mElemIndex++] * vElem0;
		vElems[1] = _elements[mElemIndex++] * vElem0;
		vElems[2] = _elements[mElemIndex++] * vElem0;
		vElems[3] = _elements[mElemIndex++] * vElem0;

		vElems[0] += _elements[mElemIndex++] * vElem1;
		vElems[1] += _elements[mElemIndex++] * vElem1;
		vElems[2] += _elements[mElemIndex++] * vElem1;
		vElems[3] += _elements[mElemIndex++] * vElem1;

		vElems[0] += _elements[mElemIndex++] * vElem2;
		vElems[1] += _elements[mElemIndex++] * vElem2;
		vElems[2] += _elements[mElemIndex++] * vElem2;
		vElems[3] += _elements[mElemIndex++] * vElem2;

		vElems[0] += _elements[mElemIndex++] * vElem3;
		vElems[1] += _elements[mElemIndex++] * vElem3;
		vElems[2] += _elements[mElemIndex++] * vElem3;
		vElems[3] += _elements[mElemIndex++] * vElem3;
	}
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

