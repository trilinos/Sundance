/* @HEADER@ */
// ************************************************************************
//
//                              Sundance
//                 Copyright (2005) Sandia Corporation
//
// Copyright (year first published) Sandia Corporation.  Under the terms
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
// retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Kevin Long (krlong@sandia.gov),
// Sandia National Laboratories, Livermore, California, USA
//
// ************************************************************************
/* @HEADER@ */

#include "SundanceBox2D.hpp"
#include "SundancePoint.hpp"
#include "SundanceDefs.hpp"
#include "SundancePolygon2D.hpp"

using namespace Sundance;

Box2D::Box2D(double px, double py, double ox, double oy, double a1, double a2,  bool flipD) :
	CurveBase(1, a1, a2, flipD), px_(px), py_(py), ox_(ox), oy_(oy)
{
}

Box2D::~Box2D()
{
}

Expr Box2D::getParams() const
{
	// return the parameters of the box
	return Expr(List(px_, py_, ox_ , oy_));
}

double Box2D::curveEquation(const Point& evalPoint) const
{
	int verb = 0;
	TEST_FOR_EXCEPTION(evalPoint.dim() != 2, std::runtime_error ,
			"Box2D::curveEquation() evaluation point dimension must be 2");

	// calculate the distance compared to the middle point
	double distX =  fabs(px_ + 0.5*ox_ - evalPoint[0]) - 0.5*ox_;
	double distY =  fabs(py_ + 0.5*oy_ - evalPoint[1]) - 0.5*oy_;
	distX = (distX > distY) ? distX : distY ;

	SUNDANCE_OUT(verb > 3, " Box2D::curveEquation for:" << evalPoint << " is: " << distX);

	return flipDomains_*distX;
}

void Box2D::returnIntersect(const Point& start, const Point& end, int& nrPoints,
		Array<double>& result) const
{

	Array<Point> t;
	returnIntersectPoints(start, end, nrPoints, t);

	result.resize(nrPoints);

	// Return coordinates instead of t values
	for (int i = 0; i < nrPoints; i++)
	{
		Point tmp( end[0]-t[i][0]/(end[0]-start[0]) , end[1]-t[i][1]/(end[1]-start[1]) );
		result[i] =  sqrt( tmp*tmp );
	}
}

void Box2D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
    int verb = 0;
    double num_zero = 1e-14;
    // first implementation
	//TEST_FOR_EXCEPTION( true , std::runtime_error , "Box2D::returnIntersect() not implemented yet");
    nrPoints = 0;
	// calc cut in X direction
	if ( fabs(start[0] - end[0]) < num_zero ){
        //
	} else {
	   SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints , 1 ")
	   // y = a*x + b
       double a = (end[1] - start[1]) / (end[0] - start[0]);
       double b = start[1] - a*start[0];
       double y1 = a*px_ + b;
       double y2 = a*(px_ + ox_) + b;
       Point p1( px_ , y1 );
       Point p2( px_ + ox_ , y2 );
       SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  px:" << px_ << " px+ox:" << px_+ox_ );
       SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  py:" << py_ << " py+oy:" << py_+oy_ );
       SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  p1:" << p1 << " p2:" << p2 );
       SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  start:" << start << " end:" << end );
       // todo: compare to numerical zero
       if ( (px_ <= p1[0] ) && ( px_ + ox_ >= p1[0] ) && (start[0] <= p1[0] ) && ( end[0]  >= p1[0] ) &&
    		(py_ <= p1[1] ) && ( py_ + oy_ >= p1[1] ) && (start[1] <= p1[1] ) && ( end[1]  >= p1[1] )){
    	   result.resize(nrPoints+1);
    	   result[nrPoints] = p1; //(p1[1] - end[1])*(start[1] - end[1]) / (start[0] - end[0]);
    	   nrPoints++;
    	   SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  ADD p1:" << p1 << " nrPoints:" << nrPoints);
       }
       if ( (px_ <= p2[0] ) && ( px_ + ox_ >= p2[0] ) && (start[0] <= p2[0] ) && ( end[0]  >= p2[0] ) &&
    		(py_ <= p2[1] ) && ( py_ + oy_ >= p2[1] ) && (start[1] <= p2[1] ) && ( end[1]  >= p2[1] )){
    	   result.resize(nrPoints+1);
    	   result[nrPoints] = p2; //(p2[1] - end[1])*(start[1] - end[1]) / (start[0] - end[0]);
    	   nrPoints++;
    	   SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  ADD p2:" << p2 << " nrPoints:" << nrPoints);
       }
	}

	// calc cut in Y direction
	if ( fabs(start[1] - end[1]) < num_zero ){
        //
	} else {
	   SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints , 2 ")
       // x = a*y + b
       double a = (end[0] - start[0]) / (end[1] - start[1]);
       double b = start[0] - a*start[1];
       double x1 = a*py_ + b;
       double x2 = a*(py_ + oy_) + b;
       Point p1( x1 , py_ );
       Point p2( x2 , py_ + oy_ );
       SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  px:" << px_ << " px+ox:" << px_+ox_ );
       SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  py:" << py_ << " py+oy:" << py_+oy_ );
       SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  p1:" << p1 << " p2:" << p2 );
       SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  start:" << start << " end:" << end );
       // todo: compare to numerical zero
       if ( (px_ <= p1[0] ) && ( px_ + ox_ >= p1[0] ) && (start[0] <= p1[0] ) && ( end[0]  >= p1[0] ) &&
    		(py_ <= p1[1] ) && ( py_ + oy_ >= p1[1] ) && (start[1] <= p1[1] ) && ( end[1]  >= p1[1] ) && (nrPoints < 2)){
    	   result.resize(nrPoints+1);
    	   result[nrPoints] = p1; //(p1[1] - end[1])*(start[1] - end[1]) / (start[0] - end[0]);
    	   nrPoints++;
    	   SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  ADD p1:" << p1 << " nrPoints:" << nrPoints);
       }
       if ( (px_ <= p2[0] ) && ( px_ + ox_ >= p2[0] ) && (start[0] <= p2[0] ) && ( end[0]  >= p2[0] ) &&
    		(py_ <= p2[1] ) && ( py_ + oy_ >= p2[1] ) && (start[1] <= p2[1] ) && ( end[1]  >= p2[1] ) && (nrPoints < 2)){
    	   result.resize(nrPoints+1);
    	   result[nrPoints] = p2; //(p2[1] - end[1])*(start[1] - end[1]) / (start[0] - end[0]);
    	   nrPoints++;
    	   SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints  ADD p2:" << p2 << " nrPoints:" << nrPoints);
       }
	}
	SUNDANCE_OUT(verb > 3, " Box2D::returnIntersectPoints END , nrPoints:" << nrPoints
			<< " , start:" << start << " , end:" << end);
}


const RCP<CurveBase> Box2D::getPolygon(const Mesh& mesh ,double resolution) const
{

   int nrXP = ::ceil(ox_/resolution);
   double stepX = ox_ / (double) nrXP;
   int nrYP = ::ceil(oy_/resolution);
   double stepY = oy_ / (double) nrYP;

   int verb = 0;

   // start bottom left
   Array<Point> points(2*(nrXP) + 2*nrYP + 1);
   int actualPoint = 0;

   // bottom
   for (int i = 0 ; i < nrXP ; i++){
	   Point tmp( px_ + stepX*(double)i , py_);
	   SUNDANCE_MSG3( verb , " Box2D::getPolygon add index " << actualPoint << " tmp=" << tmp );
	   points[actualPoint++] = tmp;
   }
   // right
   for (int i = 0 ; i < nrYP ; i++){
	   Point tmp( px_ + ox_ , py_ +  stepY*(double)i );
	   SUNDANCE_MSG3( verb , " Box2D::getPolygon add index " << actualPoint << " tmp=" << tmp );
	   points[actualPoint++] = tmp;
   }

   // add top right corner
   points[actualPoint++] = Point(px_ + ox_ , py_ + oy_ );
   SUNDANCE_MSG3( verb , " Box2D::getPolygon add index " << actualPoint << " tmp=" << points[actualPoint-1] );

   // top
   for (int i = nrXP-1 ; i >= 0 ; i--){
	   Point tmp( px_ + stepX*(double)i , py_ + oy_ );
	   SUNDANCE_MSG3( verb , " Box2D::getPolygon add index " << actualPoint << " tmp=" << tmp );
	   points[actualPoint++] = tmp;
   }
   //left
   for (int i = nrYP-1 ; i > 0 ; i--){
	   Point tmp( px_ , py_ + stepY*(double)i);
	   SUNDANCE_MSG3( verb , " Box2D::getPolygon add index " << actualPoint << " tmp=" << tmp );
	   points[actualPoint++] = tmp;
   }
   points.resize(actualPoint);

   // return the polygon
   return rcp(new Polygon2D( mesh , points , _alpha1 , _alpha2 , (flipDomains_ < 0) ));
}
