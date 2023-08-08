#include "snlKnotVector.h"

snlKnotVector::~snlKnotVector ()
{
	if(_knots) delete[] _knots;
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
		_knots = new knot[size];

		// Copy supplied data into array.
		if(knotArrayToUse)
			for(index = 0; index < size; index ++)
				*( _knots + index) = *( knotArrayToUse + index);
	}
	else
		_knots = knotArrayToUse;

	_vectorSize = size;

	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		_deg = degree;
	}
	else
	{
		_deg = SNL_KNOT_VECTOR_MAX_DEG;
	}

	_kvType = snlKnotVectorType;
}

snlKnotVector::snlKnotVector(knot startVal, knot endVal, unsigned numKnots, int degree)
{
	unsigned index;

	_kvType = open;

	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		_deg = degree;
	}
	else
	{
		_deg = SNL_KNOT_VECTOR_MAX_DEG;
	}

	// Calculate Spacing between knots.
	knot step = (endVal - startVal) / (knot)(numKnots - (2 * degree) - 1);

	// Create new array.
	_knots = new knot[numKnots];

	// Fill knot vector with data.
	for(index = 0; index < numKnots; index ++)
	{
		if(index < (unsigned) degree)
		{
			_knots[index] = startVal;
		}
		else if(index >(numKnots - 1 - degree))
		{
			_knots[index] = endVal;
		}
		else
		{
			_knots[index] = startVal +(step * ((knot)(index - degree)));
		}
	}

	_vectorSize = numKnots;
}

snlKnotVector::snlKnotVector(int size, int degree, knot* params)
{
	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		_deg = degree;
	}
	else
	{
		_deg = SNL_KNOT_VECTOR_MAX_DEG;
	}

	_kvType = open;

	_vectorSize = size + degree + 1;

	_knots = new knot[_vectorSize];

	int index;

	// Start clamp.
	for(index = 0; index <= degree; index ++)
		_knots[index] = 0.0;

	// End clamp.
	for(index = size; index < (int) _vectorSize; index ++)
		_knots[index] = 1.0;

	// Internal knots.
	for(index = 1; index < size - degree; index ++)
	{
		knot sum = 0.0;

		for(int paramIndex = index; paramIndex < index + degree; paramIndex ++)
			sum += params[paramIndex];

		_knots[index + degree] = sum / (double) degree;
	}
}

snlKnotVector::snlKnotVector(int degree)
{
	_kvType = open;

	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		_deg = degree;
	}
	else
	{
		_deg = SNL_KNOT_VECTOR_MAX_DEG;
	}

	_vectorSize =(degree + 1) * 2;

	_knots = new knot[_vectorSize];

	int upperIndex = degree + 1;

	for(int index = 0; index <= degree; index ++)
	{
		_knots[index] = 0.0;
		_knots[upperIndex ++] = 1.0;
	}
}

snlKnotVector& snlKnotVector::operator=(const snlKnotVector& knotVectToCopy)
{
	if(this != &knotVectToCopy)
	{
		if(_knots) delete[] _knots;

		__copyFrom(knotVectToCopy);
	}

	return *this;
}

void snlKnotVector::__copyFrom(const snlKnotVector& vector)
{
	_vectorSize = vector._vectorSize;

	_knots = new knot[_vectorSize];

	_deg = vector._deg;

	for(unsigned index = 0; index < _vectorSize; index ++)
		_knots[index] = vector._knots[index];

	_kvType = vector.getType();
}

unsigned snlKnotVector::findSpan(knot param) const
{
	unsigned count;
	unsigned span = 0;

	if(param > _knots[_vectorSize - 1]) param = _knots[_vectorSize - 1];

	if(param == _knots[_vectorSize - 1])
	{
		// Allow clamped end value to be a valid parameter.
		// Not strictly correct but works better in practice.

		// Step backwards through knot array until first non-zero length span found.
		for( count =(_vectorSize - 1); count > 0 ; count --)
		{
			if(param <= _knots[count] && param > _knots[count - 1])
				span = count - 1;
		}
	}
	else
	{
		for( count = 0; count <(_vectorSize - 1); count ++)
		{
			if(param >= _knots[count] && param < _knots[count + 1])
				span = count;
		}
	}

	return span;
}

knot snlKnotVector::getKnotVal(unsigned index) const
{
	return _knots[index];
}

const knot* snlKnotVector::getKnotPtr(unsigned index)
{
	return _knots + index;
}

unsigned snlKnotVector::getSize() const
{
	return _vectorSize;
}

int snlKnotVector::getDegree()
{
	return _deg;
}

void snlKnotVector::setDegree(int degree)
{
	if(degree < SNL_KNOT_VECTOR_MAX_DEG)
	{
		_deg = degree;
	}
	else
	{
		_deg = SNL_KNOT_VECTOR_MAX_DEG;
	}
}

bool snlKnotVector::equals(const snlKnotVector& knotVect) const
{
	if(_deg != knotVect._deg) return false;

	if(_vectorSize != knotVect._vectorSize) return false;

	for(unsigned index = 0; index < _vectorSize; index ++)
		if(_knots[index] != knotVect._knots[index]) return false;

	return true;
}

void snlKnotVector::insertKnot(knot param, int numTimes)
{
	unsigned        index;

	if(! _knots) return;

	unsigned span = findSpan(param);

	knot* newKnots = new knot[_vectorSize + numTimes];

	// Copy up to insertion point.
	for(index = 0; index <= span; index ++)
		newKnots[index] = _knots[index];

	// Add in new knot.
	for(int count = 0; count < numTimes; count ++)
		newKnots[span + count + 1] = param;

	// Copy rest of old knots to new vector.
	for(index = span + numTimes + 1; index < _vectorSize + numTimes; index ++)
		newKnots[index] = _knots[index - numTimes];

	delete[] _knots;

	_knots = newKnots;

	_vectorSize += numTimes;
}

void snlKnotVector::removeKnot(unsigned spanIndex)
{
	unsigned index;

	knot * newKnots = new knot[_vectorSize - 1];

	// Copy up to removal point.
	for(index = 0; index < spanIndex; index ++)
		newKnots[index] = _knots[index];

	// Copy remainder of knots. Skip knot to be removed.
	for(index = spanIndex; index < _vectorSize - 1; index ++)
		newKnots[index] = _knots[index + 1];

	delete[] _knots;

	_knots = newKnots;

	_vectorSize --;
}

void snlKnotVector::grow(unsigned bySize)
{
	knot* newKnots = new knot[_vectorSize + bySize];

	for(unsigned index = 0; index < _vectorSize; index ++)
		newKnots[index] = _knots[index];

	delete[] _knots;

	_knots = newKnots;

	_vectorSize += bySize;
}

void snlKnotVector::increaseMultiplicity(unsigned spanIndex, int numKnotsToAdd)
{
	unsigned index;

	if(!_knots) return;

	knot param = _knots[spanIndex];

	knot* newKnots = new knot[_vectorSize + numKnotsToAdd];

	// Copy up to insertion point.
	for(index = 0; index <= spanIndex; index ++)
		newKnots[index] = _knots[index];

	// Add in new knots.
	for(index = spanIndex + 1; index <= spanIndex + numKnotsToAdd; index ++)
		newKnots[index] = param;

	// Copy rest of old knots to new vector.
	for(index = spanIndex + numKnotsToAdd + 1; index < _vectorSize + numKnotsToAdd; index ++)
		newKnots[index] = _knots[index - numKnotsToAdd];

	delete[] _knots;

	_knots = newKnots;

	_vectorSize += numKnotsToAdd;
}

int snlKnotVector::getType() const
{
	return _kvType;
}

const knot* snlKnotVector::getKnotArray()
{
	return _knots;
}

knot snlKnotVector::max() const
{
	return _knots[_vectorSize - 1];
}

knot snlKnotVector::min() const
{
	return _knots[0];
}

unsigned snlKnotVector::numSpans() const
{
	unsigned numSpans = 0;

	for(unsigned index = 0; index <(_vectorSize - 1); index ++)
	{
		if(_knots[index + 1] > _knots[index]) numSpans++;
	}

	return numSpans;
}

unsigned snlKnotVector::firstSpan() const
{
	for(unsigned index = 0; index <(_vectorSize - 1); index ++)
	{
		if(_knots[index + 1] > _knots[index]) return index;
	}

	return 0;
}

unsigned snlKnotVector::nextSpan(unsigned spanIndex) const
{
	for(unsigned index = spanIndex + 1; index <(_vectorSize - 1); index ++)
	{
		if(_knots[index + 1] > _knots[index]) return index;
	}

	return 0;
}

unsigned snlKnotVector::previousSpan(unsigned spanIndex) const
{
	for(unsigned index = spanIndex - 1; index >= 0; index --)
	{
		if(_knots[index + 1] > _knots[index]) return index;
	}

	return 0;
}

int snlKnotVector::findMultiplicity(unsigned index) const
{
	// Find last index that has required knot value.

	unsigned cIndex = index;

	while(_knots[cIndex] == _knots[cIndex + 1]) cIndex ++;

	// Count multiples backwards.

	int multi = 1;

	while(_knots[cIndex] == _knots[cIndex - 1])
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
		end = _vectorSize - 1;
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
			newKnots[index] = _knots[copyFrom];

		newKnots[0] = param;
	}
	else
	{
		for(unsigned index = 0; index < newSize - 1; index ++, copyFrom ++)
			newKnots[index] = _knots[copyFrom];

		newKnots[newSize - 1] = param;
	}

	delete[] _knots;

	_knots = newKnots;

	_vectorSize = newSize;
}

void snlKnotVector::reparameterise(knot startKnot, knot endKnot)
{
	double oldLength = _knots[_vectorSize - 1] - _knots[0];
	double newLength = endKnot - startKnot;

	double oldStartKnot = _knots[0];

	for(unsigned index = 0; index < _vectorSize; index ++)
		_knots[index] =(( _knots[index] - oldStartKnot) / oldLength) * newLength + startKnot;
}

void snlKnotVector::reverse()
{
	unsigned midPoint = _vectorSize / 2;

	unsigned swapIndex = _vectorSize - 1;

	for(unsigned index = 0; index < midPoint; index ++)
	{
		knot trans = _knots[index];  // Transfer value.
		_knots[index] = _knots[swapIndex];
		_knots[swapIndex --] = trans;
	}
}

void snlKnotVector::join(snlKnotVector* knotVector)
{
	if(knotVector -> _deg != _deg) return;

	unsigned newSize = _vectorSize +(knotVector -> _vectorSize) - _deg - 2;

	unsigned oldSize = _vectorSize;

	grow(newSize - _vectorSize);

	// Copy new points into knot array.

	unsigned copyFromIndex = _deg + 1;  // The appended vector loses deg + 1 knots from the start of it's array.

	// This knot vector looses one knot at end of the array

	for(unsigned index = oldSize - 1; index < newSize; index ++)
		_knots[index] = knotVector -> _knots[copyFromIndex ++];
}

void snlKnotVector::evalBasis(knot param, basis* retArray)
{
	if(param == _knots[0] && _kvType == open)
	{
		// First basis function value is 1 the rest are zero.

		retArray[0] = 1.0;

		for(int index = 1; index < _deg + 1; index ++) retArray[index] = 0.0;

		return;
	}

	if(param == _knots[_vectorSize - 1] && _kvType == open)
	{
		// Last basis function value is 1 the rest are zero.

		retArray[_deg] = 1.0;

		for(int index = 0; index < _deg; index ++) retArray[index] = 0.0;

		return;
	}

	unsigned spanIndex = findSpan(param);

	// Using the maximum size allowed saves on memory allocation overhead.
	// Keeping them local in the stack is thread safe as well.
	basis right[SNL_KNOT_VECTOR_MAX_NUM_BASIS_VALS];
	basis left[SNL_KNOT_VECTOR_MAX_NUM_BASIS_VALS];

	basis saved, temp;
	unsigned index, level;

	retArray[0] = 1.0;

	for(level = 1; level <= (unsigned) _deg; level ++)
	{
		left[level] = param - _knots[spanIndex + 1 - level];
		right[level] = _knots[spanIndex + level] - param;

		saved = 0.0;

		for(index = 0; index < level; index ++)
		{
			temp = retArray[index] / (right[index + 1] + left[level - index]);
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
	basis right[SNL_KNOT_VECTOR_MAX_NUM_BASIS_VALS];
	basis left[SNL_KNOT_VECTOR_MAX_NUM_BASIS_VALS];
	basis derivSaved[SNL_KNOT_VECTOR_MAX_DEG];

	basis saved, temp;

	unsigned index, level;
	int count;

	unsigned spanIndex = findSpan(param);

	unsigned overFlow;  // Just in case deriv is bigger than degree.

	if(deriv > _deg)
	{
		overFlow = deriv - _deg;
		deriv = _deg;
	}
	else
	{
		overFlow = 0;
	}

	retArray[0] = 1.0;

	for(level = 1; level <= (unsigned) _deg; level ++)
	{
		left[level] = param - _knots[spanIndex + 1 - level];
		right[level] = _knots[spanIndex + level] - param;

		saved = 0.0;

		for(count = 0; count < deriv; count ++)
			derivSaved[count] = 0.0;

		for(index = 0; index < level; index ++)
		{
			temp = retArray[index] / (right[index + 1] + left[level - index]);

			retArray[index]  = saved + right[index + 1] * temp;

			saved = left[level - index] * temp;

			// Process first order derivatives as needed.
			if(level > (unsigned)(_deg - deriv))
			{
				retArray[index +(( _deg - level + 1) * (_deg + 1))] = level * (derivSaved[0] - temp);
				derivSaved[0] = temp;
			}

			// Process other order derivatives.
			for(count = _deg - level + 2; count <= deriv; count ++)
			{
				temp = retArray[index + (count *(_deg + 1))] / (right[index + 1] + left[level - index]);

				retArray[index +(count * (_deg + 1))] = level * (derivSaved[count - 1] - temp);

				derivSaved[count - 1] = temp;
			}
		}

		retArray[level] = saved;

		// Add last first order derivative at this level.
		if(level > (unsigned)(_deg - deriv))
		{
			retArray[level + (( _deg - level + 1) * (_deg + 1))] = level * derivSaved[0];
		}

		// Add last other order derivatives at this level.
		for(count = _deg - level + 2; count <= deriv; count ++)
		{
			retArray[index +(count * (_deg + 1))] = level * derivSaved[count - 1];
		}
	}

	if(overFlow)
	{
		for(index =(deriv + 1) * (_deg + 1); index <(deriv + overFlow + 1) * (_deg + 1); index ++)
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
	unsigned numEqns = _deg - multi + 1;

	knot rParam = getKnotVal(span);

	double* alpha = new double[numEqns];

	unsigned count = 0;

	for(unsigned index = span - _deg; index <=(span - multi); index ++)
	{
		alpha[count ++]  =(rParam - (_knots[index])) /(_knots[index + _deg + 1] - _knots[index]);
	}

	return alpha;
}

void snlKnotVector::print()
{
	cout << "degree: " << _deg << " Knots: ";

	for(unsigned index = 0; index < _vectorSize; index ++)
		cout << _knots[index] << "  ";

	cout << "\n";
}

void snlKnotVector::print_cpp()
{
	cout << "{ ";

	for(unsigned index = 0; index < _vectorSize; index ++)
	{
		cout << _knots[index];

		if(index < _vectorSize - 1)
			cout << ", ";
	}

	cout << " };";
}

