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

    /** 0 - line , 1 - polygon */
    static void setMethod(int methodNr) { method_option = methodNr; }

private :

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

    static void get3DRealCoordsOnSurf(const Point &refP ,
    		                          const Array<Point>& cellPoints,
    		                          const Array<int> &triangles ,
    		                          const int nrTriag ,
    		                          const Array<Point> edgeIntersectPoint,
    		                          int &triagIndex ,
    		                          Point &realPoint);

    static void getCurveQuadPoints_polyline(
                                   int  maxCellLID ,
                                   CellType  maxCellType ,
                                   const Mesh& mesh ,
                                   const Array<Point>& cellPoints,
                                   const ParametrizedCurve& paramCurve,
                                   const QuadratureFamily& quad ,
                                   Array<Point>& curvePoints ,
                                   Array<Point>& curveDerivs ,
                                   Array<Point>& curveNormals  );

    /** the int "flag" which chooses among different curve discretization methods*/
    static int method_option;

};

}

#endif /* SUNDANCECURVEINTEGRALCALC_HPP_ */
