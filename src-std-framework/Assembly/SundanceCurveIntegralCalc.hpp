/*
 * SundanceCurveIntagralCalc.hpp
 *
 *  Created on: Sep 2, 2010
 *      Author: benk
 */

#ifndef SUNDANCECURVEINTEGRALCALC_HPP_
#define SUNDANCECURVEINTEGRALCALC_HPP_

#include "SundanceParametrizedCurve.hpp"
#include "SundanceOut.hpp"
#include "SundancePoint.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceCellType.hpp"
#include "SundanceMesh.hpp"

namespace Sundance {

/** Class to compute the intersection/quadrature points of a cell with a curve in 2/3D */
class CurveIntegralCalc {
public:

	/** Empty Ctor*/
	CurveIntegralCalc();

	virtual ~CurveIntegralCalc() {;}

    /**
         [in]  maxCellType <br>
         [in]  const Array<Point>& cellPoints (the physical points of the cell) <br>
         [in]  paramCurve <br>
         [in]  Quadraturefamily , quadrature 1D for curve, or 2D for surface <br>
         [out] Array<Point>& curvePoints <br>
         [out] Array<Point>& curveDerivs  <br>
         [out] Array<Point>& curveNormals <br> */

    static void getCurveQuadPoints(CellType  maxCellType ,
    		                       int maxCellLID ,
    		                       const Mesh& mesh ,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );

    /** this method shows if special methods for polygon can be used */
    static bool usePolygoneCurveQuadrature(
    		const CellType& cellType ,
    		const Mesh& mesh ,
    		const ParametrizedCurve& paramCurve);

    /** test if we are in 3D for the Brick surface quadrature */
    static bool use3DSurfQuadrature(
    		const CellType& cellType ,
    		const Mesh& mesh ){
    	if ( (mesh.cellType(mesh.spatialDim()) == BrickCell) && (cellType == BrickCell) ) { return true; }
    	else {return false;}
    }

private :

    /** the simplest method which considers only the intersection points on the
     * facets, and inside the element we use */
    static void getCurveQuadPoints_line(
                                   int  maxCellLID ,
                                   CellType  maxCellType ,
                                   const Mesh& mesh ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );

    /** projects the point from 2D surf to 3D points of the element (needed for surf integrals)*/
    static void get3DRealCoordsOnSurf(const Point &refP ,
    		                          const Array<Point>& cellPoints,
    		                          const Array<int> &triangles ,
    		                          const int nrTriag ,
    		                          const Array<Point> edgeIntersectPoint,
    		                          int &triagIndex ,
    		                          Point &realPoint);


    /** polygons in 1D curves, uses the information from the polygon */
    static void getCurveQuadPoints_polygon(
                                   int  maxCellLID ,
                                   CellType  maxCellType ,
                                   const Mesh& mesh ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );

    /** special method for 3D surface quadrature similar to the one in getCurveQuadPoints_line ,
     * but it uses special Quadrature class, up to 4 triangles per brick cell*/
    static void getSurfQuadPoints(
                                   int  maxCellLID ,
                                   CellType  maxCellType ,
                                   const Mesh& mesh ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );
};

}

#endif /* SUNDANCECURVEINTEGRALCALC_HPP_ */
