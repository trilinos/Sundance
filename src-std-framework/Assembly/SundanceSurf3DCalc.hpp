/*
 * SundanceSurf3DCalc.hpp
 *
 *  Created on: Oct 25, 2011
 *      Author: benk
 */

#ifndef SUNDANCESURF3DCALC_HPP_
#define SUNDANCESURF3DCALC_HPP_

#include "SundanceParametrizedCurve.hpp"
#include "SundanceOut.hpp"
#include "SundancePoint.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceCellType.hpp"
#include "SundanceMesh.hpp"

namespace Sundance {

/** This class contains the computational methods which are necessary for the 3D surface and cut-cell
 * calculations.*/
class SundanceSurf3DCalc {
public:

	/** empty Ctor */
	SundanceSurf3DCalc() {;}

	/** empty Dtor */
	virtual ~SundanceSurf3DCalc() {;}

	/** method to compute the intersection surface between one curve and a brick cell.
	 * The resulted surface is represented by a triangle surface.
	 * @param maxCellType [IN]
	 * @param maxCellLID [IN]
	 * @param mesh [IN]
	 * @param paramCurve [IN] the surface object
	 * @param intersectPoints [OUT] the array with the intersection points in the reference coordinate
	 * @param brickPoints [OUT] 12 points of the brick cell
	 * @param edgeIndex [OUT] the intersection points are by default on the edges of the brick
	 * @param triangleIndex [OUT] contains the triangles, the size of this vector is 3*nrTriags. */
    static void getCurveQuadPoints(CellType  maxCellType ,
    		                       int maxCellLID ,
    		                       const Mesh& mesh ,
								   const ParametrizedCurve& paramCurve,
								   const Array<Point>& brickPoints,
								   Array<Point>& intersectPoints ,
								   Array<int>& edgeIndex,
								   Array<int>& triangleIndex);


	static const int edegIndex[12][2];

	static const int faceEdges[6][4];

	static const int edgeFaces[12][2];

	static const int edgeNeighboredges[12][12];

};


}
#endif /* SUNDANCESURF3DCALC_HPP_ */
