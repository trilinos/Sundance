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
			    	if (nrPoints == 1){
			    		if (total_points == 0) startPoint = result[0];
			    		else endPoint = result[0];
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point:" << result[0]);
						total_points += nrPoints;
			    	}
					SUNDANCE_MSG3(verb, tabs << "ind:" << jj << ", nr Int points :" << nrPoints
							<< " , start:" << startPoint << ", end:"<< endPoint);
					// if we have more than one intersection point then just ignore that edge
					//TEST_FOR_EXCEPTION( nrPoints > 1 ,
					//		RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , QuadCell one edge " << jj << " , can have only one intersection point" );
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
			    	if (nrPoints == 1){
			    		if (total_points == 0) startPoint = result[0];
			    		else endPoint = result[0];
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point:" << result[0]);
						total_points += nrPoints;
			    	}
			    	// if we have more than one intersection point then just ignore that edge
					//TEST_FOR_EXCEPTION( nrPoints > 1 ,
					//		RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , TriangleCell one edge " << jj << " , can have only one intersection point" );
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
// ====================2D curve integral in 3D context =============================

		  // take the intersection with each edge,
		  // from these intersection points construct one 3D (convex) polygon (from triangle to hexadron)
		  //  - maybe divide the polygons into many triangles
		  // map the quadrature points somehow inside this surface

		  // here we consider a simple surface, as it is in bilinear interpolation
		  Tabs tabs;
		  switch (maxCellType){
		    case BrickCell:{
		    	int edegIndex[12][2] = { {0,1} , {0,2} , {0,4} , {1,3} , {1,5} , {2,3} , {2,6} , {3,7} ,
		    			                 {4,5} , {4,6} , {5,7} , {6,7} };
		    	int faceEdges[6][4]  = { {0,1,3,5} , {0,2,4,8} , {1,2,6,9} , {3,4,7,10} , {5,6,7,11} , {8,9,10,11} };
		    	// loop over the edges
		    	int t_inters_points = 0;
		    	int nrPoints = 0;
		    	Array<Point> edgeIntersectPoint(t_inters_points);
		    	Array<int> edgeIndex(t_inters_points);
		    	Array<int> edgePointI(12,-1);

		    	// loop over the edges and get the intersection points
		    	for (int jj = 0 ; jj < 12 ; jj++ ){
			    	Array<Point> result(0);
					paramCurve.returnIntersectPoints(cellPoints[edegIndex[jj][0]], cellPoints[edegIndex[jj][1]], nrPoints , result);
					// test if we have intersection point
			    	if (nrPoints > 0){
			    		edgeIntersectPoint.resize( t_inters_points + 1 );
			    		edgeIndex.resize( t_inters_points + 1 );
			    		edgeIntersectPoint[t_inters_points] = result[0];
			    		edgeIndex[t_inters_points] = jj;
			    		edgePointI[jj] = t_inters_points;
			    		SUNDANCE_MSG3(verb, tabs << "found Int. point [" << t_inters_points <<"]=" << result[0]);
			    		t_inters_points++;
			    	}
					SUNDANCE_MSG3(verb, tabs << "ind:" << jj << ", nr Int points :" << nrPoints );
					TEST_FOR_EXCEPTION( nrPoints > 1 ,
							RuntimeError , " CurveIntegralCalc::getCurveQuadPoints , QuadCell one edge " << jj << " , can have only one intersection point" );
		    	} // getting intersection points


		    	// -------- get the point which is nearest to the other points ---------
		    	Array<double> totDist(t_inters_points);
		    	// calculate the distance from each point to each point
		    	SUNDANCE_MSG3(verb, tabs << " Calculating distances ");
		    	for (int i = 0 ; i < t_inters_points ; i++){
		    		double sum = 0.0;
		    		for (int j = 0 ; j < t_inters_points ; j++){
		    			sum = sum + ::sqrt(edgeIntersectPoint[i]*edgeIntersectPoint[j]);
		    		}
		    		totDist[i] = sum;
		    	}
		    	// find the point which is mostly in the middle
		    	int minIndex = 0;
		    	for (int i = 1 ; i < t_inters_points ; i++){
		    		if (totDist[i] < totDist[minIndex]) minIndex = i;
		    	}
		    	SUNDANCE_MSG3(verb, tabs << " minIndex = " << minIndex );


		    	// --------- discover the polygon and the triangles-----------
		    	int intersTriangleInd = 0;
		    	int nrTriangles = t_inters_points-2;
		    	Array<int> triangles(3*nrTriangles,-1);
		    	Array<double> triangle_area(nrTriangles,0.0);
		    	double total_area = 0.0;
		    	SUNDANCE_MSG3(verb, tabs << " nrTriangles = " << nrTriangles );
		    	int intPointProc = -1;
		    	int intPointProc_lastNeigh = -1;

// -------------------------------------
		    	// get one neighboring intersection point of the "minIndex", this will be the starting direction of the polygon
		    	for (int jj = 0 ; jj < 6 ; jj++ )
		    	{
		    		int nrIpointPerFace = 0;
		    		int firstPI = -1 , secondPI = -1;
		    		// loop over each edge in one face
		    		for (int ii = 0 ; ii < 4 ; ii++){
		    			// if this edge has one intersection point
		    			if ( edgePointI[faceEdges[jj][ii]] >= 0){
		    				if (nrIpointPerFace > 0){
		    					secondPI = edgePointI[faceEdges[jj][ii]];
		    				}else{
		    					firstPI = edgePointI[faceEdges[jj][ii]];
		    				}
		    				TEST_FOR_EXCEPTION( nrIpointPerFace > 2 , RuntimeError,"getCurveQuadPoints , nrIpointPerFace > 2 , " << nrIpointPerFace );
		    				nrIpointPerFace++;
		    			}
		    		}
		    		TEST_FOR_EXCEPTION( nrIpointPerFace == 1 , RuntimeError,"getCurveQuadPoints , nrIpointPerFace == 1 , " << nrIpointPerFace );
		    		// Found the first neighboring point of "minIndex"
		    		if (( nrIpointPerFace > 1) && ((minIndex == firstPI) || (minIndex == secondPI) )){
			    			SUNDANCE_MSG3(verb, tabs << " Found  starting line:" << firstPI << " -> " << secondPI);
			    			intPointProc_lastNeigh = minIndex;
			    			if (minIndex == firstPI){
			    				intPointProc = secondPI;
			    			}else{
			    				intPointProc = firstPI;
			    			}
			    			// once we found then break the current for loop
			    			break;
		    		}
		    	}
		        TEST_FOR_EXCEPTION( intPointProc < 0 , RuntimeError,"getCurveQuadPoints , intPointProc < 0 , " << intPointProc );
		    	// store this as act. point and get the next one(which is not neighbor to "minIndex") which is not the one which we already had
		        // here we create all triangles, by traversing the polygon in one direction, so that the mapping will be continous
		        for (int pI = 0 ; pI < nrTriangles; pI++)
		        {
		        	// for each new
			    	for (int jj = 0 ; jj < 6 ; jj++ )
			    	{
			    		int nrIpointPerFace = 0;
			    		int firstPI = -1 , secondPI = -1;
			    		// loop over each edge in one face
			    		for (int ii = 0 ; ii < 4 ; ii++){
			    			// if this edge has one intersection point
			    			if ( edgePointI[faceEdges[jj][ii]] >= 0){
			    				if (nrIpointPerFace > 0){
			    					secondPI = edgePointI[faceEdges[jj][ii]];
			    				}else{
			    					firstPI = edgePointI[faceEdges[jj][ii]];
			    				}
			    				TEST_FOR_EXCEPTION( nrIpointPerFace > 2 , RuntimeError,"getCurveQuadPoints , nrIpointPerFace > 2 , " << nrIpointPerFace );
			    				nrIpointPerFace++;
			    			}
			    		}
			    		TEST_FOR_EXCEPTION( nrIpointPerFace == 1 , RuntimeError,"getCurveQuadPoints , nrIpointPerFace == 1 , " << nrIpointPerFace );
			    		// condition to find the neighbor of "intPointProc" but not the one which has been already found
			    		if (( nrIpointPerFace > 1) && ((intPointProc == firstPI) || (intPointProc == secondPI))
			    				&& ((intPointProc_lastNeigh != firstPI) && (intPointProc_lastNeigh != secondPI)))
			    		{
				    			SUNDANCE_MSG3(verb, tabs << " Found next line:" << firstPI << " -> " << secondPI);
				    			if (intPointProc == firstPI){
				    				intPointProc = secondPI;
					    			intPointProc_lastNeigh = firstPI;
				    			}else{
				    				intPointProc = firstPI;
					    			intPointProc_lastNeigh = secondPI;
				    			}
				    			// add triangle
			    				triangles[intersTriangleInd] = minIndex;
			    				triangles[intersTriangleInd+1] = intPointProc_lastNeigh;
			    				triangles[intersTriangleInd+2] = intPointProc;
			    				SUNDANCE_MSG3(verb, tabs << " Found triangle:" << minIndex << " , " << firstPI << " , " << secondPI);
			    				Point v1 = edgeIntersectPoint[firstPI] - edgeIntersectPoint[minIndex];
			    				Point v2 = edgeIntersectPoint[secondPI] - edgeIntersectPoint[minIndex];
			    				Point areaV = cross(v1,v2);
			    				SUNDANCE_MSG3(verb, tabs << " TriangINdex = " << intersTriangleInd << " , area = " << ::sqrt(areaV*areaV));
			    				triangle_area[ intersTriangleInd/3 ] = ::sqrt(areaV*areaV) / (double)2.0 ; // div by 2 to get the are of the triangle
			    				total_area = total_area + triangle_area[ intersTriangleInd/3 ];
			    				intersTriangleInd = intersTriangleInd + 3;
				    			// once we found the next triangle then break the current for loop
				    			break;
			    		}
			    	}
			    }
// ---------------------------
		    	SUNDANCE_MSG3(verb, tabs << " END Found triangle  total_area=" << total_area);


		    	// if we have only 1 triangle (3 intersection points) then add to the longest edge one middle point
		    	// so that we will have at least 2 triangles
		    	if (t_inters_points < 4){
		    		// 0,1,2
		    		SUNDANCE_MSG3(verb, tabs << " Add one additional point to have at least 2 triangles ");
		    		int longEdI1 = -1 , longEdI2 = -1;
		    		for (int ii = 0 ; ii < t_inters_points ; ii++)
		    			if ( (ii != minIndex) ){
		    				longEdI1 = ii;
		    				break;
		    			}
		    		for (int ii = 0 ; ii < t_inters_points ; ii++)
		    			if ( (ii != minIndex) && (ii != longEdI1) ){
		    				longEdI2 = ii;
		    				break;
		    			}
		    		TEST_FOR_EXCEPTION( (longEdI1 < 0)||(longEdI2 < 0), RuntimeError," (longEdI1 < 0)||(longEdI2 < 0):" << longEdI1 << "," <<longEdI2);

		    		//add the middle point
		    		edgeIntersectPoint.resize(t_inters_points+1);
		    		edgeIntersectPoint[t_inters_points] = edgeIntersectPoint[longEdI1]/2.0 + edgeIntersectPoint[longEdI2]/2.0;

		    		int new_p = t_inters_points;
		    		t_inters_points += 1;
		    		nrTriangles += 1;

		    		// we'll have two new triangles
    				triangles.resize(3*(nrTriangles));
    				triangles[0] = minIndex;
    				triangles[1] = longEdI1;
    				triangles[2] = new_p;
    				triangles[3] = minIndex;
    				triangles[4] = new_p;
    				triangles[5] = longEdI2;
                    // the are will be half of the original one
    				triangle_area[0] = triangle_area[0] / 2.0;
    				triangle_area.resize(2);
    				triangle_area[1] = triangle_area[0];
		    	}

		    	// now map the quadrature points to the real cell (triangles)
				Array<Point> quadPoints;
				Array<double> quadWeights;
				quad.getPoints( QuadCell , quadPoints , quadWeights );
				int nr_quad_points = quadPoints.size();

				// --------- the quadrature points per cell ------------
				curvePoints.resize(nr_quad_points);
				Array<int> triangleInd(nr_quad_points,-1);
				SUNDANCE_MSG3(verb, tabs << " setting reference quadrature points" );
				for (int ii = 0 ; ii < nr_quad_points ; ii++) {
					//SUNDANCE_MSG3(verb, tabs << " nr:" << ii << " Quad quad point:" << quadPoints[ii] << ", line:" << (endPoint - startPoint));
					int tInd = -1;
					Point tmp(0.0,0.0,0.0);
					get3DRealCoordsOnSurf( quadPoints[ii] , cellPoints , triangles ,
					        nrTriangles , edgeIntersectPoint, tInd , tmp );
					triangleInd[ii] = tInd;
					curvePoints[ii] = tmp;
					//set the real quad points
					SUNDANCE_MSG3(verb, tabs << " quad point nr:" << ii << " = " << curvePoints[ii]);
				}

				// ------- The derivatives at points, for this simple method are simple and constant over the whole quadrature -----
				curveDerivs.resize(nr_quad_points);
				SUNDANCE_MSG3(verb, tabs << " setting surface values points" );
				for (int ii = 0 ; ii < nr_quad_points ; ii++) {
					// set the vector according to formula (only the vector noem will be taken in account)
					int tInd = triangleInd[ii] ;
					curveDerivs[ii] = Point( nrTriangles*triangle_area[tInd] , 0.0, 0.0 );
					//curveDerivs[ii] = Point( total_area , 0.0, 0.0 ); //total_area
					SUNDANCE_MSG3(verb, tabs << " curve Derivatives point nr:" << ii << " = " << curveDerivs[ii]);
				}


				// calculate the norms to the curve
				curveNormals.resize(nr_quad_points);
				SUNDANCE_MSG3(verb, tabs << " setting normal values at points" );
				for (int ii = 0 ; ii < nr_quad_points ; ii++) {
					int tInd = triangleInd[ii];
					// get the normal to the triangle
				    Point norm_vec = cross( edgeIntersectPoint[triangles[3*tInd+1]] - edgeIntersectPoint[triangles[3*tInd]] ,
				    		edgeIntersectPoint[triangles[3*tInd+2]] - edgeIntersectPoint[triangles[3*tInd]]);
				    double relative_length = 1e-2 * ::sqrt(norm_vec*norm_vec);
				    norm_vec = norm_vec / ::sqrt(norm_vec*norm_vec);
				    if ( paramCurve.curveEquation( edgeIntersectPoint[triangles[3*tInd]] + relative_length*norm_vec ) > 0 ) {
					    curveNormals[ii] = -norm_vec;
				    }
				    else {
					    curveNormals[ii] = norm_vec;
				    }
				    SUNDANCE_MSG3(verb, tabs << " curve Normals point nr:" << ii << " = " << curveNormals[ii]);
				}

                break;
		    } // ----------- end brick cell --------------
		    case TetCell:{
		    	// todo: later when this will be needed
		    	TEST_FOR_EXCEPTION( true, RuntimeError,"CurveIntagralCalc::getCurveQuadPoints , not implemented for " << maxCellType );
		    	break;
		    }
		    default:{
		    	TEST_FOR_EXCEPTION( true, RuntimeError,"CurveIntagralCalc::getCurveQuadPoints , not implemented for " << maxCellType );
		    }
		  }
		  SUNDANCE_MSG3(verb, tabs << " END CurveIntagralCalc::getCurveQuadPoints");
	  } break;

//========================== END 2D curve in 3D context ==========================

	  default: {
        // throw exception
		TEST_FOR_EXCEPTION( true, RuntimeError,"CurveIntagralCalc::getCurveQuadPoints , curve dimension must be 1 or two ");
	  }
	}
}


void CurveIntegralCalc::get3DRealCoordsOnSurf(const Point &refP ,
		const Array<Point>& cellPoints,
        const Array<int> &triangles ,
        const int nrTriag ,
        const Array<Point> edgeIntersectPoint,
        int &triagIndex ,
        Point &realPoint){

	Tabs tabs;
	int verb = 0;
	realPoint[0] = 0.0; realPoint[1] = 0.0; realPoint[2] = 0.0;
	SUNDANCE_MSG3(verb, tabs << " start get3DRealCoordsOnSurf refP="<<refP );
	Point v1;
	Point v2;
	// these matrixes were calculated in octave
	double case1[2][4] = {{1 ,-1  ,0   ,1},{ 1 , 0 ,  -1  , 1}};
	double case2[3][4] = {{1 ,-2, 0, 2} , { 2 , -2, -1, 2},{ 1 , 0, -1, 1}};
	double case3[4][4] = {{1, -2,0 , 2} , { 2 , -2, -1,2},{ 2, -1,-2, 2},{ 2 ,0, -2,  1}};
	double *refInv = 0;
	switch (nrTriag){
		case 2:{
			if (refP[0] >= refP[1]){ // first triangle
				triagIndex = 0; v1 = Point(1.0,0.0); v2 = Point(1.0,1.0);
			}else{// second triangle
				triagIndex = 1; v1 = Point(1.0,1.0); v2 = Point(0.0,1.0);
			}
			refInv = case1[triagIndex];
			break;
		}
		case 3:{
			if (refP[0] >= 2*refP[1]){ // first triangle
				triagIndex = 0; v1 = Point(1.0,0.0); v2 = Point(1.0,0.5);
			}else if ((refP[0] <= 2*refP[1]) && (refP[0] >= refP[1])){ // second triangle
				triagIndex = 1; v1 = Point(1.0,0.5); v2 = Point(1.0,1.0);
			}else{ // third triangle
				triagIndex = 2; v1 = Point(1.0,1.0); v2 = Point(0.0,1.0);
			}
			refInv = case2[triagIndex];
			break;
		}
		case 4:{
			if (refP[0] >= 2*refP[1]){ // first triangle
				triagIndex = 0; v1 = Point(1.0,0.0); v2 = Point(1.0,0.5);
			}else if ((refP[0] <= 2*refP[1]) && (refP[0] >= refP[1])){// second triangle
				triagIndex = 1; v1 = Point(1.0,0.5); v2 = Point(1.0,1.0);
			}else if ((refP[0] <= refP[1]) && (2*refP[0] >= refP[1])){// third triangle
				triagIndex = 2; v1 = Point(1.0,1.0); v2 = Point(0.5,1.0);
			}else{ // fourth triangle
				triagIndex = 3;v1 = Point(0.5,1.0); v2 = Point(0.0,1.0);
			}
			refInv = case3[triagIndex];
			break;
		}
		default:{
			TEST_FOR_EXCEPTION( true, RuntimeError,"get3DRealCoordsOnSurf , too many triangles " << nrTriag);
		}
	}
	// get the three points of the triangle
	Point p0 = edgeIntersectPoint[triangles[3*triagIndex]];
	Point p1 = edgeIntersectPoint[triangles[3*triagIndex+1]];
	Point p2 = edgeIntersectPoint[triangles[3*triagIndex+2]];
	SUNDANCE_MSG3(verb, tabs << "get3DRealCoordsOnSurf p0="<<p0<<" , p1="<<p1<<" , p2="<<p2);
	SUNDANCE_MSG3(verb, tabs << "get3DRealCoordsOnSurf v1=" << v1 << " , v2=" << v2 << " , triagIndex=" << triagIndex);

	// transfor to the reference 2D coordinates
	double ref_x = refInv[0] * refP[0] + refInv[1]*refP[1];
	double ref_y = refInv[2] * refP[0] + refInv[3]*refP[1];
	SUNDANCE_MSG3(verb, tabs << " after traf to ref_x ="<<ref_x<<" , ref_y="<<ref_y );

	// transform from 2D ref to 3D triangle real coordinate
	realPoint = p0 + ref_x*(p1-p0) + ref_y*(p2-p0);
	SUNDANCE_MSG3(verb, tabs << " end get3DRealCoordsOnSurf , realPoint="<<realPoint );
	SUNDANCE_MSG3(verb, tabs << " end get3DRealCoordsOnSurf cellPoints[0]="<<cellPoints[0]<<" , cellPoints[7]="<<cellPoints[7] );

	// transform the point in real coordinates back to reference 3D coordinates
	// this works only for cubes (structured)
	realPoint[0] = (realPoint[0] - cellPoints[0][0]) / (cellPoints[7][0] - cellPoints[0][0]);
	realPoint[1] = (realPoint[1] - cellPoints[0][1]) / (cellPoints[7][1] - cellPoints[0][1]);
	realPoint[2] = (realPoint[2] - cellPoints[0][2]) / (cellPoints[7][2] - cellPoints[0][2]);
	SUNDANCE_MSG3(verb, tabs << " end get3DRealCoordsOnSurf triagIndex="<<triagIndex<<" , realPoint="<<realPoint );
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
