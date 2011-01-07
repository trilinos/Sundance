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

#include "SundancePolygon2D.hpp"
#include "SundancePoint.hpp"
#include "SundanceMesh.hpp"
#include "SundanceDefs.hpp"

#include <iostream>
#include <fstream>

using namespace Sundance;

int Polygon2D::intersectionEdge_ = -1;

Polygon2D::Polygon2D(const Mesh& mesh , const Array<Point>& points , double a1, double a2) :
	CurveBase(1, a1, a2), mesh_(mesh)
{
	int verb = 0;
	// just store the input points
	polyPoints_.resize(points.size());
	SUNDANCE_MSG3( verb , "Polygon2D() Ctor nrPoints=" << points.size() );
	for (int j = 0 ; j < points.size() ; j++ ){
		polyPoints_[j] = points[j];
		SUNDANCE_MSG3( verb , " point[" << j << "] = " << polyPoints_[j]);
	}
    // get the maxCellsLID for each point
	computeMaxCellLIDs();
}

Polygon2D::Polygon2D(const Mesh& mesh , const std::string& filename , double a1, double a2) :
CurveBase(1, a1, a2) , mesh_(mesh)
{
	std::string str_tmp;
	std::ifstream myfile;
	std::vector<std::string> elems;
	int verb = 0;

	char *pEnd;
	double d1 , d2;
	int index_tmp = 0;

	SUNDANCE_MSG3( verb , "Polygon2D() read from file ");

	myfile.open( filename.c_str() , std::ios::in );
	polyPoints_.resize(index_tmp);

	// each line of the file should contain a point
	while (!myfile.eof())
	{
		getline( myfile , str_tmp );
		d1 = strtod (str_tmp.c_str() , &pEnd);
		d2 = strtod (pEnd , NULL);
		// we add the point to the polygon
		index_tmp = polyPoints_.size();
		polyPoints_.resize(index_tmp+1);
		Point tmpP(d1,d2);
		polyPoints_[index_tmp] = tmpP;
		SUNDANCE_MSG3( verb , " point[" << index_tmp << "] = " << polyPoints_[index_tmp]);
    }
	// close the file
	myfile.close();
	// compute the maxCells in which are the points
	computeMaxCellLIDs();
}

void Polygon2D::computeMaxCellLIDs(){

	int meshDim = mesh_.spatialDim();
	int nrLID = mesh_.numCells(meshDim);
	int verb = 0;

	//SUNDANCE_MSG3( verb , " Polygon2D::computeMaxCellLIDs " << polyPoints_.size() );

	// initially set everything to -1, so in parallel case we'll know which do not belong to this processor
	pointsMaxCellLID_.resize(polyPoints_.size());
	for (int j = 0 ; j < polyPoints_.size() ; j++){ pointsMaxCellLID_[j] = -1; }

	for (int cellLID = 0 ; cellLID < nrLID ; cellLID++){
		// set the LID first for all points to -1
		// run through each cell, and in each cell look for each point which is contained in that cell
		//
		switch (mesh_.cellType(meshDim)){
			case QuadCell:{
				// in the case of Quadcells we always have structured cells (Quadcells)
				int tmp;
				Point p0 = mesh_.nodePosition( mesh_.facetLID(meshDim,cellLID,0,0,tmp) );
				//Point p1 = mesh_.nodePosition( mesh_.facetLID(meshDim,cellLID,0,1,tmp) );
				//Point p2 = mesh_.nodePosition( mesh_.facetLID(meshDim,cellLID,0,2,tmp) );
				Point p3 = mesh_.nodePosition( mesh_.facetLID(meshDim,cellLID,0,3,tmp) );
				SUNDANCE_MSG3( verb , " Polygon2D::computeMaxCellLIDs cellLID =" << cellLID << " ,p0:" << p0 << " ,p3:" << p3 );
				for (int j = 0 ; j < polyPoints_.size() ; j++)
				{
	               //SUNDANCE_MSG3( verb , " test j = " << j );
                   Point& tmp = polyPoints_[j];
                   if ( (tmp[0] >= p0[0]) &&  (tmp[0] <= p3[0]) &&
                		(tmp[1] >= p0[1]) &&  (tmp[1] <= p3[1]))
                   {
                	   // store the maxCellLID
                	   //SUNDANCE_MSG3( verb , " Polygon2D::computeMaxCellLIDs j=" << j << " maxCellLID=" << polyPoints_.size() );
                	   pointsMaxCellLID_[j] = cellLID;
                   }
				}
				break;
			}
			case TriangleCell:{
				// todo:
				break;
			}
			default: { /* throw error */}
		}
	} // - from the for loop
}

Expr Polygon2D::getParams() const
{
	//todo: return the points of the polygon, only later for the optimization
	// if necessary at all ...
	return Expr(List(0.0));
}

double Polygon2D::curveEquation(const Point& evalPoint) const
{
	//int verb = 0;
	double dist = 1e+20;
	double dist_tmp = 0.0;
	double sign_tmp = 0.0;
	Point p0 , p1;
	double eps = 1e-12;

    for (int ii = 0 ; ii < polyPoints_.size() ; ii++)
    {
    	if (ii < polyPoints_.size() - 1)
    	{
        	p0 = polyPoints_[ii];
        	p1 = polyPoints_[ii+1];
    	}
    	else
    	{   // the last line from size-1 to 0
        	p0 = polyPoints_[ii];
        	p1 = polyPoints_[0];
    	}
    	Point intP(0.0,0.0);

    	// this is the formula from the equation V*U = 0 and det|P-P0,P-P1| = 0 (from paper)
    	double a[4] = { p1[0]-p0[0] , p1[1]-p0[1] ,
    			        p0[1]-p1[1] , p1[0]-p0[0] };
    	double b[2] = { evalPoint[0]*(p1[0] - p0[0]) + evalPoint[1]*(p1[1] - p0[1]) ,
    			        p0[1]*(p1[0] - p0[0]) - p0[0]*(p1[1] - p0[1]) };

    	// simple Gauss elimination
    	if (::fabs(a[0]) < eps ){
        	double fakt = -a[0]/a[2];
        	intP[1] = (b[0] + fakt*b[1])/(a[1] + fakt*a[3]);
        	intP[0] = (b[1] - a[3]*intP[1]) / a[2];
    	}
    	else{
        	double fakt = -a[2]/a[0];
        	intP[1] = (b[1] + fakt*b[0])/(a[3] + fakt*a[1]);
        	intP[0] = (b[0] - a[1]*intP[1]) / a[0];
    	}
    	// now we have the intersection point
    	//SUNDANCE_MSG3( verb , " intP= " << intP );

    	if  (  (    ((intP[0] >= p0[0]-eps) && (intP[0] <= p1[0]+eps))
    		     || ((intP[0] <= p0[0]+eps) && (intP[0] >= p1[0]-eps)) )
    		&& (    ((intP[1] >= p0[1]-eps) && (intP[1] <= p1[1]+eps))
        		 || ((intP[1] <= p0[1]+eps) && (intP[1] >= p1[1]-eps)) )
			)
    	{
    		dist_tmp = 0.0;
    	}
    	else
    	{
    		double d1 = ::sqrt((intP-p0)*(intP-p0));
    		double d2 = ::sqrt((intP-p1)*(intP-p1));
    		dist_tmp = ( d1 < d2 )? d1 : d2;
    	}
   		//SUNDANCE_MSG3( verb , "curveEquation() evalPoint=" << evalPoint << " intP=" << intP <<
   		//		       " is between points" << p0 << " and " << p1 << " dist_tmp=" << dist_tmp );
       	dist_tmp = ::sqrt(dist_tmp*dist_tmp + (intP-evalPoint)*(intP-evalPoint) );
       	Point v1 = p1 - p0;
       	Point v2 = evalPoint - p0;
       	sign_tmp = v1[0]*v2[1] - v1[1]*v2[0];
       	// if the vector product is negative then the point is inside
       	if ( sign_tmp > 0.0){
      		dist_tmp = -dist_tmp;
       	}
       	//SUNDANCE_MSG3( verb , "curveEquation() dist_tmp=" << dist_tmp );
       	// we store the absolute minimal distance
       	if (::fabs(dist) > ::fabs(dist_tmp)){
       		//SUNDANCE_MSG3( verb , "curveEquation() store distance " << dist_tmp << " prev dist: " << dist);
       		dist = dist_tmp;
       	}

    }

    //SUNDANCE_MSG3( verb , "curveEquation() valPoint=" << evalPoint << " return distance = " << dist);
	return dist;
}

void Polygon2D::returnIntersectPoints(const Point& start, const Point& end, int& nrPoints,
		Array<Point>& result) const
{
	int verb = 0;
	nrPoints = 0;
	result.resize(nrPoints);
	Point p0 , p1;
	double eps = 1e-12;

	// the edge index which will be intersecting as last
	intersectionEdge_ = -1;

	SUNDANCE_MSG3( verb , "returnIntersectPoints() start= " << start << " end="<< end );
    for (int ii = 0 ; ii < polyPoints_.size() ; ii++)
    {
    	if (ii < polyPoints_.size() - 1)
    	{
        	p0 = polyPoints_[ii];
        	p1 = polyPoints_[ii+1];
    	}
    	else
    	{   // the last line from size-1 to 0
        	p0 = polyPoints_[ii];
        	p1 = polyPoints_[0];
    	}
    	Point intP(0.0,0.0);
    	// this is the formula where twice the determinant should be zero
    	double a[4] = { p0[1]-p1[1] , p1[0]-p0[0] ,
    			        start[1]-end[1] , end[0]-start[0]};
    	double b[2] = { p0[1]*(p1[0]-p0[0]) - p0[0]*(p1[1]-p0[1]),
    			        start[1]*(end[0] - start[0]) - start[0]*(end[1] - start[1]) };
    	if ( ::fabs( a[0]*a[3]-a[1]*a[2]) < 1e-10){
    		// matrix is singular det|A|~0, no solution
    		continue;
    	}

    	// simple Gauss elimination
       	if (::fabs(a[0]) < eps ){
           	double fakt = -a[0]/a[2];
           	intP[1] = (b[0] + fakt*b[1])/(a[1] + fakt*a[3]);
           	intP[0] = (b[1] - a[3]*intP[1]) / a[2];
        }
        else
        {
           	double fakt = -a[2]/a[0];
           	intP[1] = (b[1] + fakt*b[0])/(a[3] + fakt*a[1]);
           	intP[0] = (b[0] - a[1]*intP[1]) / a[0];
        }
        // now we have the intersection point
        //SUNDANCE_MSG3( verb , "returnIntersectPoints() intP= " << intP << " p0= " << p0 << " p1= " << p1);

       	if  (      (    ((intP[0] >= p0[0]-eps) && (intP[0] <= p1[0]+eps))
        		     || ((intP[0] <= p0[0]+eps) && (intP[0] >= p1[0]-eps)) )
        		&& (    ((intP[1] >= p0[1]-eps) && (intP[1] <= p1[1]+eps))
            		 || ((intP[1] <= p0[1]+eps) && (intP[1] >= p1[1]-eps)) )
            	&& (    ((intP[0] >= start[0]-eps) && (intP[0] <= end[0]+eps))
           		     || ((intP[0] <= start[0]+eps) && (intP[0] >= end[0]-eps)) )
           		&& (    ((intP[1] >= start[1]-eps) && (intP[1] <= end[1]+eps))
               		 || ((intP[1] <= start[1]+eps) && (intP[1] >= end[1]-eps)) )
    		)
        {
       		    // the the intersection point to the results since it is between the two end points
        		//SUNDANCE_MSG3( verb , "returnIntersectPoints() found intP=" << intP << " is between points");
        		result.resize(nrPoints+1);
        		result[nrPoints] = intP;
               	nrPoints++;
               	// store the line index which
               	intersectionEdge_ = ii;
        }
    }// from the for loop
}

void Polygon2D::returnIntersect(const Point& start, const Point& end, int& nrPoints, Array<
		double>& result) const
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


void Polygon2D::writeToVTK(const std::string& filename) const
{
     std::ofstream myfile;
     myfile.open(filename.c_str());
     myfile << "# vtk DataFile Version 2.0 \n";
     myfile << "Generated by Sundance::Polygon2D Author: Janos Benk \n";
     myfile << "ASCII\n";
     myfile << "\n";
     myfile << "DATASET UNSTRUCTURED_GRID\n";
     myfile << "POINTS "<< polyPoints_.size() << " float \n";
     myfile << "\n";
     for (int ii = 0 ; ii < polyPoints_.size() ; ii++){
    	 myfile << polyPoints_[ii][0] << " " << polyPoints_[ii][1] << " 0\n";
     }
     myfile << "\n \n";
     myfile << "CELLS " << polyPoints_.size()<< " " << polyPoints_.size()*3 << "\n";
     myfile << "\n";
     for (int ii = 0 ; ii < polyPoints_.size()-1 ; ii++){
    	 myfile << "2 "<< ii << " " << ii +1 << "\n";
     }
     myfile << "2 "<< polyPoints_.size()-1 << " 0\n";
     myfile << "\n \n";
     myfile << "CELL_TYPES " << polyPoints_.size() << "\n";
     myfile << "\n";
     for (int ii = 0 ; ii < polyPoints_.size() ; ii++){
    	 myfile << "3\n";
     }
     myfile << "\n";

	 // later eventually there could be data on the polygon
     myfile.close();
}


RCP<CurveBase> Polygon2D::unite(ParametrizedCurve& c1 , ParametrizedCurve& c2)
{
   const CurveBase* pb1 = c1.ptr().get();
   const CurveBase* pb2 = c2.ptr().get();

   const Polygon2D* pol1 = dynamic_cast<const Polygon2D*>(pb1);
   const Polygon2D* pol2 = dynamic_cast<const Polygon2D*>(pb2);

   if ( (pol1 == 0) || (pol2 == 0)){
	   TEST_FOR_EXCEPTION( true , RuntimeError, "Polygon2D::unite one of the inputs is not a polygon 2D ");
   }

   // the array where the resulting points will be stored
   Array<Point> allPoints(0);

   int nrP1 = pol1->polyPoints_.size();
   int nrP2 = pol2->polyPoints_.size();
   int p1Index = 0;
   int p2Index = 0;
   int pNewIndex = 0;
   int polyIndex = 0;
   bool doLoop = true;
   int verb = 0;

   SUNDANCE_MSG3( verb , "unite() starting ... " );
   // 1) start with the frst polygon, iterate as long one point will be outside
   while (doLoop)
   {
	   SUNDANCE_MSG3( verb , "unite() testing p1Index=" << p1Index );
	  if ((pol2->curveEquation(pol1->polyPoints_[p1Index+1]) > 0) || (p1Index >= nrP1 - 1)){
		  doLoop = false;
		  // when already the first point is inside then add to the new polygon the first point
		  if (pol2->curveEquation(pol1->polyPoints_[p1Index]) > 0){
			   allPoints.resize(pNewIndex+1);
			   allPoints[ pNewIndex ] = pol1->polyPoints_[p1Index];
			   pNewIndex = pNewIndex + 1;
		  }
	  }
	  else
	  {
		  p1Index = p1Index + 1;
	  }
   }
   SUNDANCE_MSG3( verb , "unite() for starting found p1Index=" << p1Index << " pNewIndex=" << pNewIndex );

   // 2) iterate as long it is outside the other, with the current polygon
   // 3) if we have switch , then determine intersection point, and continue with the other polygon
   // 4) stop if we run out of points in the case of one polygon
   while ( (p1Index < nrP1) && (p2Index < nrP2)){
	   SUNDANCE_MSG3( verb , "unite() next point polyIndex = " << polyIndex);
	   if (polyIndex == 0)
	   { // first polygon
		   if ( (pol2->curveEquation(pol1->polyPoints_[p1Index]) > 0) &&
				(pol2->curveEquation(pol1->polyPoints_[p1Index+1]) > 0) ) {
			   // the next point is also inside and
			   allPoints.resize(pNewIndex+1);
			   allPoints[ pNewIndex ] = pol1->polyPoints_[p1Index];
			   SUNDANCE_MSG3( verb , "unite() step forward with polygon1 " << p1Index << " pNewIndex="
					    << pNewIndex << " allPoints[ pNewIndex ]=" << allPoints[ pNewIndex ]);
			   pNewIndex = pNewIndex + 1;
			   p1Index = p1Index + 1;
		   } else {
			   // calculate intersection point
			   Array<Point> intesectPoints;
			   int nrIntersectPoints;
			   pol2->returnIntersectPoints( pol1->polyPoints_[p1Index] , pol1->polyPoints_[p1Index+1] , nrIntersectPoints , intesectPoints );
			   SUNDANCE_MSG3( verb , "unite() calculate intersection of p2 with p1Index = " << p1Index << " nrP=" << nrIntersectPoints
					   << " , " << pol1->polyPoints_[p1Index] << " , " << pol1->polyPoints_[p1Index+1]
					  << " , " << intersectionEdge_ << " p2I=" << p2Index);
			   if (nrIntersectPoints > 0)
			   {
				   if (intersectionEdge_ < p2Index) {
					   p2Index = nrP2; // exit the whole while loop
					   break;
				   }
				   // add the first point and switch if we are not at the beginning
				   if (pNewIndex > 1){
					   polyIndex = 1;
					   allPoints.resize(pNewIndex+1);
					   allPoints[ pNewIndex ] = pol1->polyPoints_[p1Index];
					   pNewIndex = pNewIndex + 1;
					   p2Index = intersectionEdge_ + 1;
				   }
				   else{
					   p1Index = p1Index + 1;
				   }
				   allPoints.resize(pNewIndex+1);
				   allPoints[ pNewIndex ] = intesectPoints[nrIntersectPoints-1];
				   SUNDANCE_MSG3( verb , "unite() switch from p1 to p2 p2Index=" << p2Index << " pNewIndex=" <<
						   pNewIndex << " allPoints[ pNewIndex ]=" << allPoints[ pNewIndex ] );
				   pNewIndex = pNewIndex + 1;
				   // switch the polygon , only if this is the first point
			   }
			   else
			   {   // error
				   TEST_FOR_EXCEPTION( true , RuntimeError, "Polygon2D::unite nrIntersectPoints == 0");
			   }
		   }
	   }
	   else
	   { // second polygon outside
		   if ( (pol1->curveEquation(pol2->polyPoints_[p2Index]) > 0) &&
				(pol1->curveEquation(pol2->polyPoints_[p2Index+1]) > 0) ) {
			   // the next point is also inside and
			   allPoints.resize(pNewIndex+1);
			   allPoints[ pNewIndex ] = pol2->polyPoints_[p2Index];
			   SUNDANCE_MSG3( verb , "unite() step forward with polygon2 " << p2Index << " pNewIndex="
					    << pNewIndex << " allPoints[ pNewIndex ]=" << allPoints[ pNewIndex ]);
			   pNewIndex = pNewIndex + 1;
			   p2Index = p2Index + 1;
		   } else {
			   // calculate intersection point
			   Array<Point> intesectPoints;
			   int nrIntersectPoints;
			   pol2->returnIntersectPoints( pol2->polyPoints_[p2Index] , pol2->polyPoints_[p2Index+1] , nrIntersectPoints , intesectPoints );
			   SUNDANCE_MSG3( verb , "unite() calculate intersection of p1 with p2Index = " << p2Index << " nrP=" << nrIntersectPoints
					   << " , " << pol2->polyPoints_[p2Index] << " , " << pol2->polyPoints_[p2Index+1]
					   << " , " << intersectionEdge_ << " p1I=" << p1Index);
			   // todo: for circle some reason no intersection point is found
			   if ( (nrIntersectPoints > 0) || (intersectionEdge_ < p1Index) )
			   {
				   if (intersectionEdge_ < p1Index) {
					   p1Index = nrP1; // exit the whole while loop
					   break;
				   }
				   p1Index = intersectionEdge_ + 1;
				   allPoints.resize(pNewIndex+1);
				   allPoints[ pNewIndex ] = pol2->polyPoints_[p2Index];
				   pNewIndex = pNewIndex + 1;
				   allPoints.resize(pNewIndex+1);
				   allPoints[ pNewIndex ] = intesectPoints[nrIntersectPoints-1];
				   SUNDANCE_MSG3( verb , "unite() switch from p2 to p1 p1Index=" << p1Index << " pNewIndex=" <<
						   pNewIndex << " allPoints[ pNewIndex ]=" << allPoints[ pNewIndex ] );
				   pNewIndex = pNewIndex + 1;
				   // switch the polygon
				   polyIndex = 0;
			   }
			   else
			   {   // error
				   TEST_FOR_EXCEPTION( true , RuntimeError, "Polygon2D::unite nrINtersectPoints == 0");
			   }
		   }
	   }
   }

   // now create one polygon
   return rcp( new Polygon2D(pol1->mesh_ , allPoints , pol1->_alpha1 , pol2->_alpha2) );
}
