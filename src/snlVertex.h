#ifndef SNLVERTEX_H
#define SNLVERTEX_H

#include "snlPoint.h"
#include "snlKnotVector.h"

class snlVertex : public snlPoint
{
	public:

		snlVertex();

		snlVertex(double x, double y, double z, double w);

		snlVector& normal();

		/** Set vertex's normal. */
		void normal(snlVector& setTo);

		/** Set U parameter vertex was evaluated at. */
		void evalParamU(knot value);
		/** Get U parameter vertex was evaluated at. */
		knot evalParamU();

		/** Set V parameter vertex was evaluated at. */
		void evalParamV(knot value);
		/** Get V parameter vertex was evaluated at. */
		knot evalParamV();

		int flag; // Meaning depends on last operation performed on vertex.

	protected:

		/** Normal vector associated with vertex. */
		snlVector _norm;

		/** U Parameter value that vertex was evaluated at.*/
		knot _paramU;
		/** V Parameter value that vertex was evaluated at. Not used if a curve */
		knot _paramV;

	private:

};

#endif

