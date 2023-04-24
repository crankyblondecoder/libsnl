

#ifndef SNLVERTEXNET_H
#define SNLVERTEXNET_H

#include "snlVertex.h"
#include "snlCtrlPoint.h"

#include <iostream>
#include <cmath>

using namespace std;

/**
 * 2D Vertex Net.
 *
 * Natural orientation of the network is:
 *
 *   Maximum V
 *      |
 *      |
 *      |
 *      |
 *      0 ------- Maximum U
 *
 * You are looking at the front face. Unless you are inside the monitor ;-)
 */
class snlVertexNet
{
	public:

		snlVertexNet();
		virtual ~snlVertexNet();

		snlVertexNet ( const snlVertexNet& copyFrom );

		void vertexNet ( const snlCtrlPoint* ctrlPts, int sizeU, int sizeV );
		void vertexNet ( const snlCtrlPoint* ctrlPts, int numPts );

		void vertexNet ( const snlPoint* pts, int sizeU, int sizeV );
		void vertexNet ( const snlPoint* pts, int numPts );

		void appendRow ( const snlCtrlPoint* ctrlPts );

		snlVertex* vertex ( int index );
		snlVertex* vertex ( int U, int V );

		snlVertex* vertexes();

		int size() const;
		int sizeU() const;
		int sizeV() const;

		void calcNormals();

	private:

		snlVertex* vertex_net;  // 2D Vertex network.

		int size_u; // Size of network in first dimension.
		int size_v; // Size of network in second dimension.
};

#endif
