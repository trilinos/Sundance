/*
 * SundanceRefinementBase.cpp
 *
 *  Created on: Apr 29, 2010
 *      Author: benk
 */

#include "SundanceRefinementBase.hpp"

using namespace Sundance;

int GeometryRefinement::refine(const int cellLevel ,
		            const Point& cellPos ,
		            const Point& cellDimameter) const{

	//int verb = 5;
	// no refinement if the level is above the specified
	if (cellLevel >= level_) { return 0; }

	switch (curve_.getCurveDim()){
	case 1:{
		Point p0 = cellPos - 0.5*cellDimameter;
		Point p1 = Point( p0[0] , p0[1] + cellDimameter[1]);
		Point p2 = Point( p0[0] + cellDimameter[0] , p0[1]);
		Point p3 = cellPos + 0.5*cellDimameter;

		//SUNDANCE_OUT(verb > 3, " p0:" << p0 << " , p1:" << p1 << " , p2:" << p2 <<" , p3:" << p3);

		double v0 = curve_.curveEquation( p0 );
		double v1 = curve_.curveEquation( p1 );
		double tmp = v0;
		if (v1*tmp < 0.0) { return 1; } tmp = (v1 > 0.0?1.0:-1.0);
		double v2 = curve_.curveEquation( p2 );
		if (v2*tmp < 0.0) { return 1; } tmp = (v2 > 0.0?1.0:-1.0);
		double v3 = curve_.curveEquation( p3 );
		if (v3*tmp < 0.0) { return 1; }

		// no refinement
		return 0;
	} break;
	case 2: {
		Point p0 = cellPos - 0.5*cellDimameter;
		Point p1 = Point( p0[0] + cellDimameter[0] , p0[1] , p0[2]);
		Point p2 = Point( p0[0] , p0[1] + cellDimameter[1] , p0[2]);
		Point p3 = Point( p0[0]  + cellDimameter[0], p0[1] + cellDimameter[1] , p0[2]);
		Point p4 = Point( p0[0] , p0[1] , p0[2] + cellDimameter[2]);
		Point p5 = Point( p0[0] + cellDimameter[0] , p0[1] , p0[2] + cellDimameter[2]);
		Point p6 = Point( p0[0] , p0[1] + cellDimameter[1] , p0[2] + cellDimameter[2]);
		Point p7 = cellPos + 0.5*cellDimameter;
        // evaluate the curves on the points
		//SUNDANCE_OUT(4 > 3, "GeometryRefinement::refine p0:" << p0 << " , p1:" << p1 << " , p2:" << p2 <<" , p3:" << p3
		//		<< " p4:" << p4 << " , p5:" << p5 << " , p6:" << p6 <<" , p7:" << p7 );

		double v0 = curve_.curveEquation( p0 ) , v1 = curve_.curveEquation( p1 );
		double tmp = (v0 > 0.0?1.0:-1.0);
		if (v1*tmp < 0.0) { return 1; } tmp = (v1 > 0.0?1.0:-1.0);
		double v2 = curve_.curveEquation( p2 );
		if (v2*tmp < 0.0) { return 1; } tmp = (v2 > 0.0?1.0:-1.0);
		double v3 = curve_.curveEquation( p3 );
		if (v3*tmp < 0.0) { return 1; } tmp = (v3 > 0.0?1.0:-1.0);
		double v4 = curve_.curveEquation( p4 );
		if (v4*tmp < 0.0) { return 1; } tmp = (v4 > 0.0?1.0:-1.0);
		double v5 = curve_.curveEquation( p5 );
		if (v5*tmp < 0.0) { return 1; } tmp = (v5 > 0.0?1.0:-1.0);
		double v6 = curve_.curveEquation( p6 );
		if (v6*tmp < 0.0) { return 1; } tmp = (v6 > 0.0?1.0:-1.0);
		double v7 = curve_.curveEquation( p7 );
		if (v7*tmp < 0.0) { return 1; }
		// no refinement
		return 0;
	 } break;
	}
	// by default no refinement
	return 0;
}

