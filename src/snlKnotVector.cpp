#include "snlKnotVector.h"

snlKnotVector::~snlKnotVector ()
{
	if(__knots) delete[] __knots;
}

snlKnotVector::snlKnotVector(const snlKnotVector& vector)
{
	__copyFrom(vector);
}

snlKnotVector::snlKnotVector(knot* knotArrayToUse, unsigned size, int degree, int snlKnotVectorType, bool copy)
{
	unsigned index;

	if(copy)
	{
		// Create new knot vector array from supplied data of "size".
		__knots = new knot[size];

		// Copy supplied data into array.
		if(knotArrayToUse)
			for(index = 0; index < size; index ++)
				*( __knots + index) = *( knotArrayToUse + index);
	}
	else
		__knots = knotArrayToUse;

	__vectorSize = size;

	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		__deg = degree;
	}
	else
	{
		__deg = SNL_KNOT_VECTOR_MAX_DEG;
	}

	__kvType = snlKnotVectorType;
}

snlKnotVector::snlKnotVector(knot startVal, knot endVal, unsigned numKnots, int degree)
{
	unsigned index;

	__kvType = open;

	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		__deg = degree;
	}
	else
	{
		__deg = SNL_KNOT_VECTOR_MAX_DEG;
	}

	// Calculate Spacing between knots.
	knot step = (endVal - startVal) / (knot)(numKnots - (2 * degree) - 1);

	// Create new array.
	__knots = new knot[numKnots];

	// Fill knot vector with data.
	for(index = 0; index < numKnots; index ++)
	{
		if(index < (unsigned) degree)
		{
			__knots[index] = startVal;
		}
		else if(index >(numKnots - 1 - degree))
		{
			__knots[index] = endVal;
		}
		else
		{
			__knots[index] = startVal +(step * ((knot)(index - degree)));
		}
	}

	__vectorSize = numKnots;
}

snlKnotVector::snlKnotVector(int size, int degree, knot* params)
{
	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		__deg = degree;
	}
	else
	{
		__deg = SNL_KNOT_VECTOR_MAX_DEG;
	}

	__kvType = open;

	__vectorSize = size + degree + 1;

	__knots = new knot[__vectorSize];

	int index;

	// Start clamp.
	for(index = 0; index <= degree; index ++)
		__knots[index] = 0.0;

	// End clamp.
	for(index = size; index < (int) __vectorSize; index ++)
		__knots[index] = 1.0;

	// Internal knots.
	for(index = 1; index < size - degree; index ++)
	{
		knot sum = 0.0;

		for(int paramIndex = index; paramIndex < index + degree; paramIndex ++)
			sum += params[paramIndex];

		__knots[index + degree] = sum / (double) degree;
	}
}

snlKnotVector::snlKnotVector(int degree)
{
	__kvType = open;

	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		__deg = degree;
	}
	else
	{
		__deg = SNL_KNOT_VECTOR_MAX_DEG;
	}

	__vectorSize =(degree + 1) * 2;

	__knots = new knot[__vectorSize];

	int upperIndex = degree + 1;

	for(int index = 0; index <= degree; index ++)
	{
		__knots[index] = 0.0;
		__knots[upperIndex ++] = 1.0;
	}
}

snlKnotVector& snlKnotVector::operator=(const snlKnotVector& knotVectToCopy)
{
	if(this != &knotVectToCopy)
	{
		if(__knots) delete[] __knots;

		__copyFrom(knotVectToCopy);
	}

	return *this;
}

void snlKnotVector::__copyFrom(const snlKnotVector& vector)
{
	__vectorSize = vector.__vectorSize;

	__knots = new knot[__vectorSize];

	__deg = vector.__deg;

	for(unsigned index = 0; index < __vectorSize; index ++)
		__knots[index] = vector.__knots[index];

	__kvType = vector.getType();
}

unsigned snlKnotVector::findSpan(knot param) const
{
	unsigned count;
	unsigned span = 0;

	if(param > __knots[__vectorSize - 1]) param = __knots[__vectorSize - 1];

	if(param == __knots[__vectorSize - 1])
	{
		// Allow clamped end value to be a valid parameter.
		// Not strictly correct but works better in practice.

		// Step backwards through knot array until first non-zero length span found.
		for( count =(__vectorSize - 1); count > 0 ; count --)
		{
			if(param <= __knots[count] && param > __knots[count - 1])
				span = count - 1;
		}
	}
	else
	{
		for( count = 0; count <(__vectorSize - 1); count ++)
		{
			if(param >= __knots[count] && param < __knots[count + 1])
				span = count;
		}
	}

	return span;
}

knot snlKnotVector::getKnotVal(unsigned index) const
{
	return __knots[index];
}

const knot* snlKnotVector::getKnotPtr(unsigned index)
{
	return __knots + index;
}

unsigned snlKnotVector::getSize() const
{
	return __vectorSize;
}

int snlKnotVector::getDegree()
{
	return __deg;
}

void snlKnotVector::setDegree(int degree)
{
	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		__deg = degree;
	}
	else
	{
		__deg = SNL_KNOT_VECTOR_MAX_DEG;
	}
}

bool snlKnotVector::equals(const snlKnotVector& knotVect) const
{
	if(__deg != knotVect.__deg) return false;

	if(__vectorSize != knotVect.__vectorSize) return false;

	for(unsigned index = 0; index < __vectorSize; index ++)
		if(__knots[index] != knotVect.__knots[index]) return false;

	return true;
}

void snlKnotVector::insertKnot(knot param, int numTimes)
{
	unsigned        index;

	if(! __knots) return;

	unsigned span = findSpan(param);

	knot* newKnots = new knot[__vectorSize + numTimes];

	// Copy up to insertion point.
	for(index = 0; index <= span; index ++)
		newKnots[index] = __knots[index];

	// Add in new knot.
	for(int count = 0; count < numTimes; count ++)
		newKnots[span + count + 1] = param;

	// Copy rest of old knots to new vector.
	for(index = span + numTimes + 1; index < __vectorSize + numTimes; index ++)
		newKnots[index] = __knots[index - numTimes];

	delete[] __knots;

	__knots = newKnots;

	__vectorSize += numTimes;
}

void snlKnotVector::removeKnot(unsigned spanIndex)
{
	unsigned index;

	knot * newKnots = new knot[__vectorSize - 1];

	// Copy up to removal point.
	for(index = 0; index < spanIndex; index ++)
		newKnots[index] = __knots[index];

	// Copy remainder of knots. Skip knot to be removed.
	for(index = spanIndex; index < __vectorSize - 1; index ++)
		newKnots[index] = __knots[index + 1];

	delete[] __knots;

	__knots = newKnots;

	__vectorSize --;
}

void snlKnotVector::grow(unsigned bySize)
{
	knot* newKnots = new knot[__vectorSize + bySize];

	for(unsigned index = 0; index < __vectorSize; index ++)
		newKnots[index] = __knots[index];

	delete[] __knots;

	__knots = newKnots;

	__vectorSize += bySize;
}

void snlKnotVector::increaseMultiplicity(unsigned spanIndex, int numKnotsToAdd)
{
	unsigned index;

	if(!__knots) return;

	knot param = __knots[spanIndex];

	knot* newKnots = new knot[__vectorSize + numKnotsToAdd];

	// Copy up to insertion point.
	for(index = 0; index <= spanIndex; index ++)
		newKnots[index] = __knots[index];

	// Add in new knots.
	for(index = spanIndex + 1; index <= spanIndex + numKnotsToAdd; index ++)
		newKnots[index] = param;

	// Copy rest of old knots to new vector.
	for(index = spanIndex + numKnotsToAdd + 1; index < __vectorSize + numKnotsToAdd; index ++)
		newKnots[index] = __knots[index - numKnotsToAdd];

	delete[] __knots;

	__knots = newKnots;

	__vectorSize += numKnotsToAdd;
}

int snlKnotVector::getType() const
{
	return __kvType;
}

const knot* snlKnotVector::getKnotArray()
{
	return __knots;
}

knot snlKnotVector::max() const
{
	return __knots[__vectorSize - 1];
}

knot snlKnotVector::min() const
{
	return __knots[0];
}

unsigned snlKnotVector::numSpans() const
{
	unsigned numSpans = 0;

	for(unsigned index = 0; index <(__vectorSize - 1); index ++)
	{
		if(__knots[index + 1] > __knots[index]) numSpans++;
	}

	return numSpans;
}

unsigned snlKnotVector::firstSpan() const
{
	for(unsigned index = 0; index <(__vectorSize - 1); index ++)
	{
		if(__knots[index + 1] > __knots[index]) return index;
	}

	return 0;
}

unsigned snlKnotVector::nextSpan(unsigned spanIndex) const
{
	for(unsigned index = spanIndex + 1; index <(__vectorSize - 1); index ++)
	{
		if(__knots[index + 1] > __knots[index]) return index;
	}

	return 0;
}

unsigned snlKnotVector::previousSpan(unsigned spanIndex) const
{
	for(unsigned index = spanIndex - 1; index >= 0; index --)
	{
		if(__knots[index + 1] > __knots[index]) return index;
	}

	return 0;
}

int snlKnotVector::findMultiplicity(unsigned index) const
{
	// Find last index that has required knot value.

	unsigned cIndex = index;

	while(__knots[cIndex] == __knots[cIndex + 1]) cIndex ++;

	// Count multiples backwards.

	int multi = 1;

	while(__knots[cIndex] == __knots[cIndex - 1])
	{
		multi ++;

		if(cIndex == 1)
			break;
		else
			cIndex --;
	}

	return multi;
}

int snlKnotVector::findMultiplicity(knot param) const
{
	// Find span parameter belongs to.
	unsigned span = findSpan(param);

	// Find value of knot associated with span.
	knot spanKnotVal = getKnotVal(span);

	// Return multiplicity.
	if(param == spanKnotVal)
		return findMultiplicity(span);

	return 0;
}

void snlKnotVector::truncate(knot param, bool keepLast)
{
	unsigned start, end;

	if(keepLast)
	{
		start = findSpan(param);
		start = previousSpan(start) + 1;  // Go to start of knots if more than one at param.
		end = __vectorSize - 1;
	}
	else
	{
		start = 0;
		end = findSpan(param);
	}

	// Generate new knot vector and populate.

	unsigned newSize = end - start + 2;  // Add one point for correct clamping.

	knot* newKnots = new knot[newSize];

	unsigned copyFrom = start;

	if(keepLast)
	{
		for(unsigned index = 1; index < newSize; index ++, copyFrom ++)
			newKnots[index] = __knots[copyFrom];

		newKnots[0] = param;
	}
	else
	{
		for(unsigned index = 0; index < newSize - 1; index ++, copyFrom ++)
			newKnots[index] = __knots[copyFrom];

		newKnots[newSize - 1] = param;
	}

	delete[] __knots;

	__knots = newKnots;

	__vectorSize = newSize;
}

void snlKnotVector::reparameterise(knot startKnot, knot endKnot)
{
	double oldLength = __knots[__vectorSize - 1] - __knots[0];
	double newLength = endKnot - startKnot;

	double oldStartKnot = __knots[0];

	for(unsigned index = 0; index < __vectorSize; index ++)
		__knots[index] =(( __knots[index] - oldStartKnot) / oldLength) * newLength + startKnot;
}

void snlKnotVector::reverse()
{
	unsigned midPoint = __vectorSize / 2;

	unsigned swapIndex = __vectorSize - 1;

	for(unsigned index = 0; index < midPoint; index ++)
	{
		knot trans = __knots[index];  // Transfer value.
		__knots[index] = __knots[swapIndex];
		__knots[swapIndex --] = trans;
	}
}

void snlKnotVector::join(snlKnotVector* knotVector)
{
	if(knotVector -> __deg != __deg) return;

	unsigned newSize = __vectorSize +(knotVector -> __vectorSize) - __deg - 2;

	unsigned oldSize = __vectorSize;

	grow(newSize - __vectorSize);

	// Copy new points into knot array.

	unsigned copyFromIndex = __deg + 1;  // The appended vector loses deg + 1 knots from the start of it's array.

	// This knot vector looses one knot at end of the array

	for(unsigned index = oldSize - 1; index < newSize; index ++)
		__knots[index] = knotVector -> __knots[copyFromIndex ++];
}

void snlKnotVector::evalBasis(knot param, basis* retArray)
{
	if(param == __knots[0] && __kvType == open)
	{
		// First basis function value is 1 the rest are zero.

		retArray[0] = 1.0;

		for(int index = 1; index < __deg + 1; index ++) retArray[index] = 0.0;

		return;
	}

	if(param == __knots[__vectorSize - 1] && __kvType == open)
	{
		// Last basis function value is 1 the rest are zero.

		retArray[__deg] = 1.0;

		for(int index = 0; index < __deg; index ++) retArray[index] = 0.0;

		return;
	}

	unsigned spanIndex = findSpan(param);

	// Using the maximum size allowed saves on memory allocation overhead.
	// Keeping them local in the stack is thread safe as well.
	basis right[SNL_KNOT_VECTOR_MAX_DEG_PLUS_1];
	basis left[SNL_KNOT_VECTOR_MAX_DEG_PLUS_1];

	basis saved, temp;
	unsigned index, level;

	retArray[0] = 1.0;

	for(level = 1; level <= (unsigned) __deg; level ++)
	{
		left[level] = param - __knots[spanIndex + 1 - level];
		right[level] = __knots[spanIndex + level] - param;

		saved = 0.0;

		for(index = 0; index < level; index ++)
		{
			temp = retArray[index] /(right[index + 1] + left[level - index]);
			retArray[index] = saved + right[index + 1] * temp;
			saved = left[level - index] * temp;
		}

		retArray[level] = saved;
	}
}

void snlKnotVector::evalBasisDeriv(knot param, int deriv, basis* retArray)
{
	// Assume retArray is size [(deriv + 1) * (__deg + 1)];

	// Using the maximum size allowed saves on memory allocation overhead.
	// Keeping them local in the stack is thread safe as well.
	basis right[SNL_KNOT_VECTOR_MAX_DEG_PLUS_1];
	basis left[SNL_KNOT_VECTOR_MAX_DEG_PLUS_1];
	basis derivSaved[SNL_KNOT_VECTOR_MAX_DEG];

	basis saved, temp;

	unsigned index, level;
	int count;

	unsigned spanIndex = findSpan(param);

	unsigned overFlow;  // Just in case deriv is bigger than degree.

	if(deriv > __deg)
	{
		overFlow = deriv - __deg;
		deriv = __deg;
	}
	else
	{
		overFlow = 0;
	}

	retArray[0] = 1.0;

	for(level = 1; level <= (unsigned) __deg; level ++)
	{
		left[level] = param - __knots[spanIndex + 1 - level];
		right[level] = __knots[spanIndex + level] - param;

		saved = 0.0;

		for(count = 0; count < deriv; count ++)
			derivSaved[count] = 0.0;

		for(index = 0; index < level; index ++)
		{
			temp = retArray[index] / (right[index + 1] + left[level - index]);

			retArray[index]  = saved + right[index + 1] * temp;

			saved = left[level - index] * temp;

			// Process first order derivatives as needed.
			if(level > (unsigned)(__deg - deriv))
			{
				retArray[index +(( __deg - level + 1) * (__deg + 1))] = level * (derivSaved[0] - temp);
				derivSaved[0] = temp;
			}

			// Process other order derivatives.
			for(count = __deg - level + 2; count <= deriv; count ++)
			{
				temp = retArray[index + (count *(__deg + 1))] / (right[index + 1] + left[level - index]);

				retArray[index +(count * (__deg + 1))] = level * (derivSaved[count - 1] - temp);

				derivSaved[count - 1] = temp;
			}
		}

		retArray[level] = saved;

		// Add last first order derivative at this level.
		if(level > (unsigned)(__deg - deriv))
		{
			retArray[level + (( __deg - level + 1) * (__deg + 1))] = level * derivSaved[0];
		}

		// Add last other order derivatives at this level.
		for(count = __deg - level + 2; count <= deriv; count ++)
		{
			retArray[index +(count * (__deg + 1))] = level * derivSaved[count - 1];
		}
	}

	if(overFlow)
	{
		for(index =(deriv + 1) * (__deg + 1); index <(deriv + overFlow + 1) * (__deg + 1); index ++)
		{
			retArray[index] = 0.0;
		}
	}
}

double* snlKnotVector::calcRemovalAlphas(unsigned span)
{
	// Find multiplicity of knot at index.
	unsigned multi = findMultiplicity(span);

	// Calculate the number of equations.
	unsigned numEqns = __deg - multi + 1;

	knot rParam = getKnotVal(span);

	double* alpha = new double[numEqns];

	unsigned count = 0;

	for(unsigned index = span - __deg; index <=(span - multi); index ++)
	{
		alpha[count ++]  =(rParam - (__knots[index])) /(__knots[index + __deg + 1] - __knots[index]);
	}

	return alpha;
}

void snlKnotVector::print()
{
	cout << "degree: " << __deg << " Knots: ";

	for(unsigned index = 0; index < __vectorSize; index ++)
		cout << __knots[index] << "  ";

	cout << "\n";
}

void snlKnotVector::print_cpp()
{
	cout << "{ ";

	for(unsigned index = 0; index < __vectorSize; index ++)
	{
		cout << __knots[index];

		if(index < __vectorSize - 1)
			cout << ", ";
	}

	cout << " };";
}

