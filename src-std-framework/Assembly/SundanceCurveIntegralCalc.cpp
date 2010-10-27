/*
 * SundanceCurveIntagralCalc.cpp
 *
 *  Created on: Sep 2, 2010
 *      Author: benk
 */

#include "SundanceCurveIntegralCalc.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCellJacobianBatch.hpp"

using namespace Sundance;

// the int "flag" which chooses among different curve discretization methods
// default value is zero, which is the simplest discretization (by a simple line or plane)
int CurveIntegralCalc::method_option = 0;


CurveIntegralCalc::CurveIntegralCalc() {
	//nothing to do
}

void CurveIntegralCalc::getCurveQuadPoints(CellType  maxCellType ,
							   int maxCellLID ,
							   const Mesh& mesh ,
							   const ParametrizedCurve& paramCurve,
							   const QuadratureFamily& quad ,
							   Array<Point>& curvePoints ,
							   Array<Point>& curveDerivs ,
							   Array<Point>& curveNormals ){

	int verb = 0;
	Tabs tabs;

    // get all the points from the maxDimCell
	int nr_point = mesh.numFacets( mesh.spatialDim() , 0 , 0);
	Array<Point> maxCellsPoints(nr_point);
	int tmp_i , point_LID;

	SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints nr points per cell: " << nr_point)
	for (int jj = 0 ; jj < nr_point ; jj++){
		point_LID = mesh.facetLID( mesh.spatialDim() , maxCellLID , 0 , jj , tmp_i );
		maxCellsPoints[jj] = mesh.nodePosition(point_LID);
		SUNDANCE_MSG3(verb, tabs << " max cell point p[" << jj << "]:"<< maxCellsPoints[jj]);
	}

	// based on the static value decide which method to choose
	switch (method_option){
	case 0:{
		   // call the simple line method
		   CurveIntegralCalc::getCurveQuadPoints_line( maxCellLID , maxCellType , mesh , maxCellsPoints , paramCurve,
							    quad , curvePoints , curveDerivs , curveNormals);
		   break;}
	case 1:{
			// later here could be more sophisticated method called instead of the simple line/plane method
		    CurveIntegralCalc::getCurveQuadPoints_polyline( maxCellLID , maxCellType , mesh , maxCellsPoints , paramCurve,
								    quad , curvePoints , curveDerivs , curveNormals);
			break; }
	}

}

void CurveIntegralCalc::getCurveQuadPoints_line(
                                   int  maxCellLID ,
                                   CellType  maxCellType ,
                                   const Mesh& mesh ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals ){

	int verb = 0;
	Tabs tabs;

    // select the dimension of the curve
	switch (paramCurve.getCurveDim()){
	  case 1: {

		  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, curve has ONE dimnesion")
		  // 1D curve integral in 2D or 3D

		  // first determine the two intersection points with the edges
		  Point startPoint(0.0,0.0);
		  Point endPoint(0.0,0.0);
	      int nrPoints , total_points = 0 ;

		  switch (maxCellType){
		    case QuadCell:{
		    	// loop over each edge and detect intersection point
		    	// there can be only one
		    	TEST_FOR_EXCEPTION( cellPoints.size() != 4 ,
		    			RuntimeError ," CurveIntegralCalc::getCurveQuadPoints , QuadCell must have 4 points" );
				SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, on QuadCell")
		    	Array<Point> result(0);
		    	int edegIndex[4][2] = { {0,1} , {0,2} , {1,3} , {2,3} };

		    	// loop over the edges
		    	for (int jj = 0 ; jj < 4 ; jj++ ){
					paramCurve.returnIntersectPoints(cellPoints[edegIndex[jj][0]], cellPoints[edegIndex[jj][1]], nrPoints , result);
					// test if we have intersection point
			    	if (nrPoints > 0){
			    		if (total_points == 0) startPoint = result[0];
			    		else endPoint = result[0];
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point:" << result[0]);
			    	}
					total_points += nrPoints;
					SUNDANCE_MSG3(verb, tabs << "ind:" << jj << ", nr Int points :" << nrPoints
							<< " , start:" << startPoint << ", end:"<< endPoint);
					TEST_FOR_EXCEPTION( nrPoints > 1 ,
							RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , QuadCell one edge " << jj << " , can have only one intersection point" );
		    	}
		    	// test if we have too much intersection points
		    	TEST_FOR_EXCEPTION( total_points > 2 ,RuntimeError , " CurveIntegralCalc::getCurveQuadPoints total_points > 2 : " << total_points );
				SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_line, end finding intersection points")
		    } break;
		    case TriangleCell:{
		    	TEST_FOR_EXCEPTION( cellPoints.size() != 3 , RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , TriangleCell must have 3 points" );
		    	Array<Point> result;
		    	int edegIndex[3][2] = { {0,1} , {0,2} , {1,2} };

		    	// loop over the edges
		    	for (int jj = 0 ; jj < 3 ; jj++ ){
					paramCurve.returnIntersectPoints(cellPoints[edegIndex[jj][0]], cellPoints[edegIndex[jj][1]], nrPoints , result);
					// test if we have intersection point
			    	if (nrPoints > 0){
			    		if (total_points == 0) startPoint = result[0];
			    		else endPoint = result[0];
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point:" << result[0]);
			    	}
					total_points += nrPoints;
					TEST_FOR_EXCEPTION( nrPoints > 1 ,
							RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , TriangleCell one edge " << jj << " , can have only one intersection point" );
		    	}
		    	// test if we have too much intersection points
		    	TEST_FOR_EXCEPTION( total_points > 2 ,RuntimeError , " CurveIntegralCalc::getCurveQuadPoints total_points > 2 : " << total_points );
		    } break;
		    default : {
		    	TEST_FOR_EXCEPTION( true , RuntimeError , "CurveIntegralCalc::getCurveQuadPoints , Unknown Cell in 2D" );
		    }
		  }
		  // test for to many intersection points
		  TEST_FOR_EXCEPTION( total_points > 2 ,RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , no much intersection points: " << total_points);

		  // -> having the intersection points now we have to determine the:
		  // quad points and the gradients at the points(which norm will be in the integration)
		  // and the normalized normal vector (normalized since we need the sin and cos for Nitsche )
		  // -> as a very simple approach we consider the curve as a line between intersection points "startPoint -> endPoint"

		  // in the X direction the line should always have increasing values
		  if (startPoint[0] > endPoint[0]){
			  Point tmp = startPoint;
			  startPoint = endPoint;
			  endPoint = tmp;
		  }
		  SUNDANCE_MSG3(verb, tabs << "start end and points , start:" << startPoint << ", end:"<< endPoint)

		  // get the quadrature points for the line
		  Array<Point> quadPoints;
		  Array<double> quadWeights;
		  quad.getPoints( LineCell , quadPoints , quadWeights );
		  int nr_quad_points = quadPoints.size();

		  // The intersection points , we distribute the points along the line
		  curvePoints.resize(nr_quad_points);
		  SUNDANCE_MSG3(verb, tabs << " setting reference quadrature points" );
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  SUNDANCE_MSG3(verb, tabs << " nr:" << ii << " Line quad point:" << quadPoints[ii] << ", line:" << (endPoint - startPoint));
			  curvePoints[ii] = startPoint + quadPoints[ii][0]*(endPoint - startPoint);

			  // we transform the intersection points to the reference cell
			  // this should work for unstructured case
			  // - get the Jacobian of the cells
			  // - get the inverse of this Jacobian
			  // - apply this Jacobian to the physical points (pos0 = (0,0))of the cell

			  Array<int> cellLID(1); cellLID[0] = maxCellLID;
		      CellJacobianBatch JBatch;
		      JBatch.resize(1, 2, 2);
			  mesh.getJacobians(2, cellLID , JBatch);
			  double pointc[2] = { curvePoints[ii][0] - cellPoints[0][0] , curvePoints[ii][1] - cellPoints[0][1]};
			  JBatch.applyInvJ(0, pointc , 1 , false);
			  curvePoints[ii][0] = pointc[0];
			  curvePoints[ii][1] = pointc[1];

			  SUNDANCE_MSG3(verb, tabs << " quad point nr:" << ii << " = " << curvePoints[ii]);
		  }


		  // The derivatives at points, for this simple method are simple and constant over the whole quadrature
		  curveDerivs.resize(nr_quad_points);
		  Point dist_point(endPoint[0] - startPoint[0] , endPoint[1] - startPoint[1]);
		  SUNDANCE_MSG3(verb, tabs << " setting derivative values points" );
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  curveDerivs[ii] = endPoint - startPoint;
			  SUNDANCE_MSG3(verb, tabs << " curve Derivatives point nr:" << ii << " = " << curveDerivs[ii]);
		  }


		  // calculate the norms to the curve
		  curveNormals.resize(nr_quad_points);

		  // this is a shorter variant
		  SUNDANCE_MSG3(verb, tabs << " setting normal values at points" );
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  Point norm_vec = Point(-curveDerivs[ii][1], curveDerivs[ii][0]);
			  norm_vec = norm_vec / ::sqrt(norm_vec*norm_vec);
			  double relative_length = 1e-2 * ::sqrt((endPoint - startPoint)*(endPoint - startPoint));
			  if ( paramCurve.curveEquation( startPoint + relative_length*norm_vec ) > 0 ) {
				  curveNormals[ii] = -norm_vec;
			  }
			  else {
				  curveNormals[ii] = norm_vec;
			  }
			  SUNDANCE_MSG3(verb, tabs << " curve Normals point nr:" << ii << " = " << curveNormals[ii]);
		  }


	  } break;

//========================== END 1D curve in 2D context ==========================

	  case 2: {
		  // 2D curve integral in 3D

		  // here we consider a simple surface, as it is in bilinear interpolation

		  TEST_FOR_EXCEPTION( true, RuntimeError,"CurveIntagralCalc::getCurveQuadPoints , surface integral not implemented yet ");

	  } break;

//========================== END 2D curve in 3D context ==========================

	  default: {
        // throw exception
		TEST_FOR_EXCEPTION( true, RuntimeError,"CurveIntagralCalc::getCurveQuadPoints , curve dimension must be 1 or two ");
	  }
	}
}


void CurveIntegralCalc::getCurveQuadPoints_polyline(
        int  maxCellLID ,
        CellType  maxCellType ,
        const Mesh& mesh ,
        const Array<Point>& cellPoints,
		   const ParametrizedCurve& paramCurve,
		   const QuadratureFamily& quad ,
		   Array<Point>& curvePoints ,
		   Array<Point>& curveDerivs ,
		   Array<Point>& curveNormals ) {

	int verb = 0;
	Tabs tabs;

    // Switch according to dimension
	switch (paramCurve.getCurveDim()){

	  case 1: {

		  // 1D curve integral in 2D
		  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline, curve is 1D")

		  // first determine the two intersection points with the edges
		  Point startPoint(0.0,0.0);
		  Point endPoint(0.0,0.0);
	      int total_points(0);
	      Array<int> nIntEdge;

		  switch (maxCellType){

		    case QuadCell:{
				SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline, on QuadCell")
		    	TEST_FOR_EXCEPTION( cellPoints.size() != 4 ,
		    			RuntimeError ,"CurveIntegralCalc::getCurveQuadPoints_polyline, QuadCell must consist of 4 points" );

		    	// Loop over each edge and detect intersection points
				Array<Point> result;
				nIntEdge.resize(4);
		    	int edgeIndex[4][2] = { {0,1} , {0,2} , {1,3} , {2,3} };

		    	//  2 __________ 3
		    	//   |     3    |
		    	//   |          |
		    	//  1|          |2
		    	//   |__________|
		    	//  0      0     1

		    	// For each edge ...
		    	for (int jj = 0 ; jj < 4 ; jj++ ) {
		    		paramCurve.returnIntersectPoints(cellPoints[edgeIndex[jj][0]], cellPoints[edgeIndex[jj][1]], nIntEdge[jj] , result);

					// Test for intersection points
			    	if (nIntEdge[jj] > 0) {
			    		if (total_points == 0) {
			    			startPoint = result[0];
				    		SUNDANCE_MSG3(verb, tabs << "Found intersection point: " << result[0] << " Take it as start point!")
			    		}
			    		else {
			    			endPoint = result[0];
				    		SUNDANCE_MSG3(verb, tabs << "Found intersection point: " << result[0] << " Take it as end point!")
			    		}
						total_points += nIntEdge[jj];
			    	}
					SUNDANCE_MSG3(verb, tabs << "Edge index: " << jj << ", Nr intersections: " << nIntEdge[jj])
					TEST_FOR_EXCEPTION( nIntEdge[jj] > 1 ,
							RuntimeError , " CurveIntegralCalc::getCurveQuadPoints_polyline, QuadCell: Edge " << jj << " , cannot handle too many intersection points" );
		    	}

		    	// Test for too many intersection points
		    	TEST_FOR_EXCEPTION( total_points > 2 ,RuntimeError , " CurveIntegralCalc::getCurveQuadPoints_polyline total_points > 2 : " << total_points );
				SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Finding intersection points finished: start point "
						<< startPoint << " end point " << endPoint)
		    } break;

		    default : {
		    	TEST_FOR_EXCEPTION( true , RuntimeError , "CurveIntegralCalc::getCurveQuadPoints_polyline , Unhandled 2D cell" );
		    }
		  }

		  // -> Intersection points given, now we have to determine:
		  // quad points and the gradients at the points (which norm will be in the integration)
		  // and the normalized normal vector (normalized since we need the sin and cos for Nitsche)
		  // As an approximation we consider a piecewise linear curve

		  // we follow the line along the abscissa
		  if (startPoint[0] > endPoint[0]) {
			  Point foo = startPoint;
			  startPoint = endPoint;
			  endPoint = foo;
		  }
		  Point dir_vec = endPoint - startPoint;

		  // get the quadrature points
		  Array<Point> quadPoints;
		  Array<double> quadWeights;
		  quad.getPoints( LineCell, quadPoints, quadWeights );
		  int nr_quad_points = quadPoints.size();

		  // Distribute quadrature points on curve
		  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Determine quadrature points on curve")

		  int nCuts(0);
		  Array<Point> curvePts;
		  curvePoints.resize(nr_quad_points);
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  Point pivotPoint = startPoint + quadPoints[ii][0] * dir_vec;
			  if ( quadPoints[ii][0] == 0. || quadPoints[ii][0] == 1. ) {
				  curvePoints[ii] = pivotPoint;

				  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Quad point " << ii << " is on line: " << curvePoints[ii])
				  continue;
			  }
			  Point norm_vec(-dir_vec[1], dir_vec[0]);
			  Point normEndPoint = pivotPoint + 2*norm_vec;
			  paramCurve.returnIntersectPoints(pivotPoint, normEndPoint, nCuts, curvePts);
			  if ( nCuts > 0 ) {
				  curvePoints[ii] = curvePts[0];

				  TEST_FOR_EXCEPTION( nCuts > 1, RuntimeError, "CurveIntegralCalc::getCurveQuadPoints_polyline: Too many intersections along search line");
				  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Quad point " << ii << " is to the left: " << curvePoints[ii])
			  }
			  else {
				  Point normEndPoint = pivotPoint - 2*norm_vec;
				  paramCurve.returnIntersectPoints(pivotPoint, normEndPoint, nCuts, curvePts);

				  TEST_FOR_EXCEPTION( nCuts == 0, RuntimeError, "CurveIntegralCalc::getCurveQuadPoints_polyline: No intersections along search line");
				  TEST_FOR_EXCEPTION( nCuts > 1, RuntimeError, "CurveIntegralCalc::getCurveQuadPoints_polyline: Too many intersections along search line");
				  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Quad point " << ii << " is to the right: " << curvePoints[ii])
				  curvePoints[ii] = curvePts[0];
			  }
		  }

		  // Set points to reference coordinates
		  startPoint[0] = -(cellPoints[0][0] - startPoint[0])/(cellPoints[3][0] - cellPoints[0][0]);
		  startPoint[1] = -(cellPoints[0][1] - startPoint[1])/(cellPoints[3][1] - cellPoints[0][1]);
		  endPoint[0] = -(cellPoints[0][0] - endPoint[0])/(cellPoints[3][0] - cellPoints[0][0]);
		  endPoint[1] = -(cellPoints[0][1] - endPoint[1])/(cellPoints[3][1] - cellPoints[0][1]);
		  for (int ii = 0 ; ii < nr_quad_points ; ii++ ) {
			  // we transform the intersection points to the reference cell
			  // this should work for unstructured case
			  // - get the Jacobian of the cells
			  // - get the inverse of this Jacobian
			  // - apply this Jacobian to the physical points (pos0 = (0,0))of the cell

			  Array<int> cellLID(1); cellLID[0] = maxCellLID;
		      CellJacobianBatch JBatch;
		      JBatch.resize(1, 2, 2);
			  mesh.getJacobians(2, cellLID , JBatch);
			  double pointc[2] = { curvePoints[ii][0] - cellPoints[0][0] , curvePoints[ii][1] - cellPoints[0][1]};
			  JBatch.applyInvJ(0, pointc , 1 , false);
			  curvePoints[ii][0] = pointc[0];
			  curvePoints[ii][1] = pointc[1];
			  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Reference quad point " << ii << " is: " << curvePoints[ii])
		  }

		  // The derivatives at points
		  curveDerivs.resize(nr_quad_points);

		  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Setting derivative values" );
		  if ( quadPoints[0][0] == 0. ) {
			  curveDerivs[0] = curvePoints[1] - curvePoints[0];
		  }
		  else {
			  curveDerivs[0] = (curvePoints[1] - startPoint)/2;
		  }
		  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Curve derivatives point nr: 0 = " << curveDerivs[0])

		  for (int ii = 1; ii < nr_quad_points-1 ; ii++) {
			  curveDerivs[ii] = (curvePoints[ii+1] - curvePoints[ii-1])/2.;

			  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Curve derivatives point nr: " << ii << " = " << curveDerivs[ii]);
		  }

		  if ( quadPoints[nr_quad_points-1][0] == 1. ) {
			  curveDerivs[nr_quad_points-1] = curvePoints[nr_quad_points-1] - curvePoints[nr_quad_points-2];
		  }
		  else {
			  curveDerivs[nr_quad_points-1] = (endPoint - curvePoints[nr_quad_points-2])/2;
		  }
		  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Curve derivatives point nr: " <<  nr_quad_points-1 << " = " << curveDerivs[nr_quad_points-1])

		  // Set the normalized normal vector at each point
		  curveNormals.resize(nr_quad_points);

		  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Determine normal values at points" )
		  for (int ii = 0 ; ii < nr_quad_points ; ii++) {
			  Point norm_vec = Point(-curveDerivs[ii][1], curveDerivs[ii][0]);
			  norm_vec = norm_vec / ::sqrt(norm_vec*norm_vec);
			  if ( paramCurve.curveEquation(curvePoints[ii]+1e-2*norm_vec) > 0 ) {
				  curveNormals[ii] = -norm_vec;
			  }
			  else {
				  curveNormals[ii] = norm_vec;
			  }
			  SUNDANCE_MSG3(verb, tabs << "CurveIntegralCalc::getCurveQuadPoints_polyline: Curve normal point nr:" << ii << " = " << curveNormals[ii]);
		  }
	  } break;

//========================== END 1D curve in 2D context ==========================

	  case 2: {
		  // 2D curve integral in 3D

		  // here we consider a simple surface, as it is in bilinear interpolation

		  TEST_FOR_EXCEPTION( true, RuntimeError,"CurveIntagralCalc::getCurveQuadPoints_polyline , surface integral not implemented yet ");

	  } break;

//========================== END 2D curve in 3D context ==========================

	  default: {
        // throw exception
		TEST_FOR_EXCEPTION( true, RuntimeError,"CurveIntagralCalc::getCurveQuadPoints_polyline , curve dimension must be 1 or two ");
	  }
	}
}
