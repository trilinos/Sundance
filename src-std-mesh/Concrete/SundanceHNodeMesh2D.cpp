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
/*
 * SundanceHNodeMesh2D.cpp
 *
 *  Created on: Sep 8, 2009
 *      Author: benk
 */

#include "SundanceHNodeMesh2D.hpp"

#include "SundanceMeshType.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceMaximalCofacetBatch.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceDebug.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceCollectiveExceptionCheck.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace std;

#define P1 0.333333333333
#define P2 0.666666666666

#define P12 0.4444444444444
#define P22 0.5555555555555

#define DM_NR_POINTS 28
#define DM_NR_EDGE   52
#define DM_NR_CELL   19

#define DM_NR_EDGE_R   47
#define DM_NR_CELL_R   17

Point HNodeMesh2D::returnPoint(0.0 , 0.0);

/** This is defined on the unit square */ //                                              9   10                   //14
double HNodeMesh2D::_point_x_coords[29] = { 0.0, 1.0, 0.0, 1.0, P1 , 0.0, P1 , P2 , P2 , 1.0, P2 , 1.0, P1 , 0.0 , P1 ,
		                                    P2 , P12 , P1 , P12 , P22 , /*20*/ P22 , P2 , P22 , P2 , P12 , P1 , P12 , P22  };

/** */                                    //                                              9   10                   //14
double HNodeMesh2D::_point_y_coords[29] = { 0.0, 0.0, 1.0, 1.0, 0.0 , P1, P1 , 0.0 , P1 , P1, P2 , P2 , P2 , P2 ,  1.0 ,
                                            1.0 , P1 , P12 , P12 , P1 , /*20*/ P12 , P12 , P22 , P22 , P22 , P22 , P2 , P2  };

/** */
int HNodeMesh2D::_edgePoints[52][2] = { {0,1} , {0,2} , {1,3} , {2,3} , {0,4} , {0,5} , {4,6} , {5,6} , {4,7} , {7,8} , {6,8},
		                                {7,1} , {1,9} , {8,9} , {8,10}, {9,11},{10,11}, {6,12},{12,10}, {5,13},{13,12},{13,2},
		                               {12,14}, {2,14},{10,15},{14,15}, {11,3}, {15,3}, {6,16}, {6,17},{16,18},{17,18},{16,19},
                                       {19,20},{18,20},{19,8} , {8,21},{20,21},{20,22},{21,23},{22,23},{18,24},{24,22},{17,25},
                                       {25,24},{25,12},{24,26},{12,26},{22,27},{26,27},{23,10},{27,10}};

/** */
int HNodeMesh2D::_cellPoints[19][4] = { {0,1,2,3} , {0,4,5,6}, {4,7,6,8}, {7,1,8,9}, {8,9,10,11}, {6,8,12,10}, {5,6,13,12},
		                               {13,12,2,14}, {12,10,14,15}, {10,11,15,3} , {6,16,17,18}, {16,19,18,20}, {19,8,20,21} ,
		                               {20,21,22,23}, {18,20,24,22}, {17,18,25,24}, {25,24,12,26}, {24,22,26,27}, {22,23,27,10} };

/** */
int HNodeMesh2D::_cellEdges[19][4] = { {0,1,2,3} , {4,5,6,7}, {8,6,9,10}, {11,9,12,13}, {13,14,15,16}, {10,17,14,18}, {7,19,17,20},
                                       {20,21,22,23}, {18,22,24,25}, {16,24,26,27} ,{28,29,30,31}, {32,30,33,34}, {35,33,36,37} ,
                                       {37,38,39,40}, {34,41,38,42}, {31,43,41,44}, {44,45,46,47}, {42,46,48,49}, {40,48,50,51} };
/** */
int HNodeMesh2D::_isCellLeaf[19] = {  0 ,   1, 1, 1, 1, 0, 1, 1, 1, 1   ,   1, 1, 1, 1, 1, 1, 1, 1, 1};

/** The edge is not leaf when is not in at least one leaf cell(2D)*/
int HNodeMesh2D::_isEdgeLeaf[52] = {  0,0,0,0,1,1,1,1,1,1  ,  1,1,1,1,1,1,1,1,1,1  ,  1,1,1,1,1,1,1,1,1,1  ,  1,1,1,1,1,1,1,1,1,1
		                           ,  1,1,1,1,1,1,1,1,1,1  ,  1,1  };

int HNodeMesh2D::_Cellindex[19] = { 1,2,3,4,6,7,8,9,10,11   ,   12, 13, 14, 15, 16,17,18,19};

int HNodeMesh2D::_Edgeindex[52] = {  4,5,6,7,8,9,10,11,12,13,14,15  ,  16,17,18,19,20,21,22,23,24,25  ,
		                             26,27,28,29,30,31,32,33,34,35  ,  36,37,38,39,40,41,42,43,44,45  ,  46,47,48,49,50,51  };

int HNodeMesh2D::_CellReindex[19] = { -1 ,   0, 1, 2, 3, -1, 4, 5, 6, 7   ,   8, 9, 10, 11, 12, 13, 14, 15, 16 };

int HNodeMesh2D::_EdgeReindex[52] = {  -1,-1,-1,-1,0,1,2,3,4,5  ,  6,7,8,9,10,11,12,13,14,15  ,  16,17,18,19,20,21,22,23,24,25  ,
		                               26,27,28,29,30,31,32,33,34,35  ,  36,37,38,39,40,41,42,43,44,45  ,  46,47  };

/** */
int HNodeMesh2D::_parentCellIndex[19] = {0 , 0,0,0,0,0,0,0,0,0  , 5,5,5,5,5,5,5,5,5};

/** */
//int HNodeMesh2D::_IndexInParentCell[19] = {0 , 0,1,2,3,4,5,6,7,8  , 0,1,2,3,4,5,6,7,8};
int HNodeMesh2D::_IndexInParentCell[19] = {0,1,2,3,5,6,7,8  , 0,1,2,3,4,5,6,7,8};

/** */
int HNodeMesh2D::_cellLevel[19] = {1 , 2,2,2,2,2,2,2,2,2  ,  3,3,3,3,3,3,3,3,3 };

/** */
int HNodeMesh2D::_pointMaxCoFacet[29][4] = {{ 1,-1,-1,-1} , {-1, 3,-1,-1} , {-1,-1, 7,-1} , {-1,-1,-1, 9} , { 1, 2,-1,-1},
		                                    { 6,-1, 1,-1} , {5 ,6 ,2 , 1} , { 3, 2,-1,-1} , {4 ,5 ,3 , 2} , {-1, 4,-1, 3},
		                                    { 9, 8, 4, 5} , {-1,9 ,-1, 4} , { 8, 7, 5, 6} , {7 ,-1, 6,-1} , {-1,-1, 8, 7},
		                                    {-1,-1, 9, 8} , {11,12,-1,-1} , {13,-1,12,-1} , {14,13,11,12} , {10,11,-1,-1},
		                                    {15,14,10,11} , {-1,15,-1,10} , {16,17,15,14} , {-1,16,-1,15} , {17,18,14,13},
		                                    {18,-1,13,-1} , {-1,-1,17,18} , {-1,-1,16,17} , {-1,-1,-1,-1} };

/** The hanging edges have on both side cells (this is the convention what I fixed )*/
int HNodeMesh2D::_edgeMaxCoFacet[52][2] ={  { 0,-1},{ 0,-1},{-1, 0},{-1, 0},{ 1,-1},{ 1,-1},{ 2, 1},{ 6, 1},{2 ,-1},{ 3, 2},
		                                    { 5, 2},{ 3,-1},{-1, 3},{ 4, 3},{ 4, 5},{-1, 4},{ 9, 4},{ 5, 6},{ 8, 5},{ 6,-1},
		                                    { 7, 6},{ 7,-1},{ 8, 7},{-1, 7},{ 9, 8},{-1, 8},{-1, 9},{-1, 9},{12, 4},{12, 6},
		                                    {11,12},{13,12},{11, 2},{10,11},{14,11},{10, 2},{ 4,10},{15,10},{15,14},{ 4,15},
		                                    {16,15},{14,13},{17,14},{13, 6},{18,13},{18, 6},{17,18},{ 8,18},{16,17},{ 8,17},
		                                    { 8,16},{ 8,16}  };

int HNodeMesh2D::_pointIsHanging[29] = { 0,0,0,0,0,0,0,0,0,0  ,  0,0,0,0,0,0,1,1,0,1 , 0,1,0,1,0,1,1,1,0,  };

int HNodeMesh2D::_edgeIsHanging[52] = { 0,0,0,0,0,0,0,0,0,0  ,  0,0,0,0,0,0,0,0,0,0 , 0,0,0,0,0,0,0,0,1,1
		                              , 0,0,1,0,0,1,1,0,0,1  ,  0,0,0,1,0,1,0,1,0,1 , 1,1  };

double HNodeMesh2D::_divFact[4] = { 1 , 3 , 9 , 27 };

double HNodeMesh2D::_returnDoubleVect[2];

int HNodeMesh2D::_returnIntVect[4];


HNodeMesh2D::HNodeMesh2D(int dim, const MPIComm& comm ,
	    const MeshEntityOrder& order)
: MeshBase(dim, comm , order),_dimension(dim), _comm(comm)
{
}

void HNodeMesh2D::createMesh(
                      double position_x,
			          double position_y,
			          double offset_x,
			          double offset_y,
			          double resolution_x,
			          double resolution_y
){
	_pos_x = position_x;
	_pos_y = position_y;
	_ofs_x = offset_x;
	_ofs_y = offset_y;
	_res_x = resolution_x;
	_res_y = resolution_y;
}


HNodeMesh2D::~HNodeMesh2D() {
}


int HNodeMesh2D::numCells(int dim) const  {
	//printf("HNodeMesh2D::numCells(int dim):   dim:%d \n" , dim );
	switch (dim){
	case 0: return DM_NR_POINTS;
	case 1: return DM_NR_EDGE_R;
	case 2: return DM_NR_CELL_R;
	}
	return 0;
}

Point HNodeMesh2D::nodePosition(int i) const {
	//SUNDANCE_VERB_HIGH("nodePosition(int i)");
	//printf("HNodeMesh2D::nodePosition(int i)   i:%d \n", i);
	HNodeMesh2D::returnPoint[0] = _pos_x + _ofs_x * _point_x_coords[i];
	HNodeMesh2D::returnPoint[1] = _pos_y + _ofs_y * _point_y_coords[i];
	return HNodeMesh2D::returnPoint;
}

const double* HNodeMesh2D::nodePositionView(int i) const {
	//printf("HNodeMesh2D::nodePositionView(int i)   i:%d \n", i);
	//SUNDANCE_VERB_HIGH("nodePosition(int i)");
	_returnDoubleVect[0] = _pos_x + _ofs_x * _point_x_coords[i];
	_returnDoubleVect[1] = _pos_y + _ofs_y * _point_y_coords[i];
	return (const double*)(&_returnDoubleVect);
}

void HNodeMesh2D::getJacobians(int cellDim, const Array<int>& cellLID,
                          CellJacobianBatch& jBatch) const
{
	  //printf("HNodeMesh2D::getJacobians  cellDim:%d  _x:%f , _y:%f \n",cellDim, _ofs_x , _ofs_y );
	  SUNDANCE_VERB_HIGH("getJacobians()");
	  TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), InternalError,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	  int nCells = cellLID.size();
	  int LID;
	  Point pnt(0.0,0.0);

	  jBatch.resize(cellLID.size(), spatialDim(), cellDim);
	  if (cellDim < spatialDim()) // they need the Jacobian of a lower dinemsional element
	  {
		   for (int i=0; i<nCells; i++)
		    {
		      double* detJ = jBatch.detJ(i);
		      switch(cellDim)
		      {
		        case 0: *detJ = 1.0;
		          break;
		        case 1:
				  LID = _Edgeindex[cellLID[i]];
			      pnt[0] = (_point_x_coords[_edgePoints[LID][1]] - _point_x_coords[_edgePoints[LID][0]])* _ofs_x;
			      pnt[1] = (_point_y_coords[_edgePoints[LID][1]] - _point_y_coords[_edgePoints[LID][0]])* _ofs_y;
		          *detJ = sqrt(pnt * pnt); // the length of the edge
		        break;
		        default:
		          TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
		            "cellDim=" << cellDim << " in HNodeMesh2D::getJacobians()");
		      }
		    }
	  }else{ // they request the complete Jacoby matrix for this bunch of elements
		    //Array<double> J(cellDim*cellDim);
		    SUNDANCE_VERB_HIGH("cellDim == spatialDim()");
		    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
		    {
		      double* J = jBatch.jVals(i);
		      switch(cellDim)
		      {
		        case 2:
				  LID = _Cellindex[cellLID[i]];
		          J[0] =  (_point_x_coords[_cellPoints[LID][0]] - _point_x_coords[_cellPoints[LID][1]])*_ofs_x; //_ofs_x / _divFact[_cellLevel[LID]];
		          J[1] = 0.0;  J[2] = 0.0;   //
		          J[3] =  (_point_y_coords[_cellPoints[LID][0]] - _point_y_coords[_cellPoints[LID][2]])*_ofs_y; //_point_x_coords[]     //_ofs_y / _divFact[_cellLevel[LID]];
		        break;
		        default:
		          TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
		            "cellDim=" << cellDim
		            << " in HNodeMesh2D::getJacobians()");
		      }
		    }
	  }
}

void HNodeMesh2D::getCellDiameters(int cellDim, const Array<int>& cellLID,
                              Array<double>& cellDiameters) const {
	 TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), InternalError,
	    "cellDim=" << cellDim << " is not in expected range [0, " << spatialDim() << "]");
	 SUNDANCE_VERB_HIGH("getCellDiameters()");
	  cellDiameters.resize(cellLID.size());
	  Point pnt(0.0,0.0);
	  int LID;
	  if (cellDim < spatialDim())
	  {
		  //printf("HNodeMesh2D::getCellDiameters(), cellDim < spatialDim() \n ");
	    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	    {
	      switch(cellDim)
	      {
	        case 0:
	             cellDiameters[i] = 1.0;
	          break;
	        case 1:  //length of the edge
			  LID = _Edgeindex[cellLID[i]];
		      pnt[0] = (_point_x_coords[_edgePoints[LID][1]] - _point_x_coords[_edgePoints[LID][0]])* _ofs_x;
		      pnt[1] = (_point_y_coords[_edgePoints[LID][1]] - _point_y_coords[_edgePoints[LID][0]])* _ofs_y;
		      cellDiameters[i] = sqrt(pnt * pnt); // the length of the edge
	        break;
	        default:
	          TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
	            "cellDim=" << cellDim << " in HNodeMesh2D::getCellDiameters()");
	      }
	    }
	  }
	  else
	  {
		  //printf("HNodeMesh2D::getCellDiameters(), cellDim == spatialDim() \n ");
	    for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	    {
	      switch(cellDim)
	      {
	        case 2:
			  LID = _Cellindex[cellLID[i]];
		      pnt[0] = (_point_x_coords[_cellPoints[LID][3]] - _point_x_coords[_cellPoints[LID][0]])* _ofs_x;
		      pnt[1] = (_point_y_coords[_cellPoints[LID][3]] - _point_y_coords[_cellPoints[LID][0]])* _ofs_y;
	          cellDiameters[i] = sqrt(pnt * pnt);
	        break;
	        default:
	          TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
	            "cellDim=" << cellDim
	            << " in HNodeMesh2D::getCellDiameters()");
	      }
	    }
	  }
}

void HNodeMesh2D::pushForward(int cellDim, const Array<int>& cellLID,
                         const Array<Point>& refQuadPts,
                         Array<Point>& physQuadPts) const {

	//printf("HNodeMesh2D::pushForward cellDim:%d\n",cellDim);
	  TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), InternalError,
	    "cellDim=" << cellDim
	    << " is not in expected range [0, " << spatialDim()
	    << "]");

	  int nQuad = refQuadPts.size();
	  Array<double> J(cellDim*cellDim);
	  Point pnt( 0.0 , 0.0 );
	  Point pnt1( 0.0 , 0.0 );

	  if (physQuadPts.size() > 0) physQuadPts.resize(0);
	  physQuadPts.reserve(cellLID.size() * refQuadPts.size());
	  for (unsigned int i=0; i<(unsigned int)cellLID.size(); i++)
	  {
	    switch(cellDim)
	    {
	      case 0: // integrate one point
	         physQuadPts.append(pnt);
	        break;
	      case 1:{ // integrate on one line
			     int LID = _Edgeindex[cellLID[i]];
		    	 pnt[0] = _pos_x + _point_x_coords[_edgePoints[LID][0]] * _ofs_x ;
		    	 pnt[1] = _pos_y + _point_y_coords[_edgePoints[LID][0]] * _ofs_y ;
		    	 pnt1[0] = (_point_x_coords[_edgePoints[LID][1]] - _point_x_coords[_edgePoints[LID][0]]) * _ofs_x ;
		    	 pnt1[1] = (_point_y_coords[_edgePoints[LID][1]] - _point_y_coords[_edgePoints[LID][0]]) * _ofs_y ;
	             for (int q=0; q<nQuad; q++) {
	               physQuadPts.append(pnt + (pnt1)*refQuadPts[q][0]);
	        	}
	      break;}
	      case 2:{
			     int LID = _Cellindex[cellLID[i]];
		    	 pnt[0] = _pos_x + _point_x_coords[_cellPoints[LID][0]] * _ofs_x;
		    	 pnt[1] = _pos_y + _point_y_coords[_cellPoints[LID][0]] * _ofs_y;
		    	 pnt1[0] = (_point_x_coords[_cellPoints[LID][3]] - _point_x_coords[_cellPoints[LID][0]] ) * _ofs_x ;
		    	 pnt1[1] = (_point_y_coords[_cellPoints[LID][3]] - _point_y_coords[_cellPoints[LID][0]] ) * _ofs_y ;
		         for (int q=0; q<nQuad; q++) {
		          	  physQuadPts.append( pnt
		           	    + Point( refQuadPts[q][0] * pnt1[0] ,
		           	    		refQuadPts[q][1] * pnt1[1] ) );
		         }
	      break;}
	      default:
	        TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
	          "in HNodeMesh2D::getJacobians()");
	    }
	  }
}

int HNodeMesh2D::ownerProcID(int cellDim, int cellLID) const  {
	 SUNDANCE_VERB_HIGH("ownerProcID()");
	 return 0;
}


int HNodeMesh2D::numFacets(int cellDim, int cellLID,
                      int facetDim) const  {
	SUNDANCE_VERB_HIGH("numFacets()");
	if (cellDim==1) { // 1 dimension
         return 2; //one line has 2 points
    }
    else if (cellDim==2) { // 2 dimensions
         return 4; //one quad has 4 edges and 4 points
    }
	return -1;
}

int HNodeMesh2D::facetLID(int cellDim, int cellLID,
                     int facetDim, int facetIndex,
                     int& facetOrientation) const  {
	facetOrientation = 0;
	//printf("HNodeMesh2D::facetLID  cellDim: %d , cellLID: %d , facetDim %d , facetIndex:%d  \n" , cellDim , cellLID , facetDim , facetIndex );
	int rnt = -1;
	if (facetDim == 0){ // return the Number/ID of a Vertex
		if (cellDim == 2 ){
		    int LID = _Cellindex[cellLID];
		    rnt = _cellPoints[LID][facetIndex];
	    }
	    else if ((cellDim==1)){
	        int LID = _Edgeindex[cellLID];
	        rnt = _edgePoints[LID][facetIndex];
	    }
	} else if (facetDim == 1){
	        int LID = _Cellindex[cellLID];
	        rnt = _cellEdges[LID][facetIndex];
			rnt = _EdgeReindex[rnt];
	       }
	//printf(" RET = %d \n" ,rnt );
	return rnt;
}


void HNodeMesh2D::getFacetLIDs(int cellDim,
                          const Array<int>& cellLID,
                          int facetDim,
                          Array<int>& facetLID,
                          Array<int>& facetSign) const {
	SUNDANCE_VERB_HIGH("getFacetLIDs()");
	//printf("PeanoMesh2D::getFacetLIDs()  cellDim:%d  cellLID.size():%d  facetDim:%d\n" , cellDim, (int)cellLID.size() , facetDim);
    int LID = 0 , cLID , facetOrientation ;
    int ptr = 0;

    int nf = numFacets(cellDim, cellLID[0], facetDim);
    facetLID.resize(cellLID.size() * nf);
    facetSign.resize(cellLID.size() * nf);
    // At this moment we just use the previous function
	for (unsigned int i = 0 ; i < (unsigned int)cellLID.size() ; i++){
		  cLID = cellLID[i];
	      for (int f=0; f<nf; f++, ptr++) {
	    	  // we use this->facetLID caz facetLID is already used as variable
			  LID = this->facetLID( cellDim, cLID, facetDim, f , facetOrientation);
			  //printf("LID:%d , cellDim:%d , cLID:%d , facetDim:%d , f:%d , facetOrientation:%d \n"
			  //	  ,LID , cellDim, cLID, facetDim, f , facetOrientation );
	          facetLID[ptr] = LID;
	          facetSign[ptr] = facetOrientation;
	      }
	}
	//printf("PeanoMesh2D::getFacetLIDs()  DONE. \n");
}


const int* HNodeMesh2D::elemZeroFacetView(int cellLID) const {
    int LID = _Cellindex[cellLID];
	_returnIntVect[0] = _cellPoints[LID][0];
	_returnIntVect[1] = _cellPoints[LID][1];
	_returnIntVect[2] = _cellPoints[LID][2];
	_returnIntVect[3] = _cellPoints[LID][3];
	return (const int*)(&_returnIntVect);
}


int HNodeMesh2D::numMaxCofacets(int cellDim, int cellLID) const  {
	//SUNDANCE_VERB_HIGH("numMaxCofacets()");
	//printf("HNodeMesh2D::numMaxCofacets():  cellDim:%d cellLID:%d \n",cellDim, cellLID );
	int rnt = -1;
	if (cellDim==0) { // 1 dimension
        int sum = 0;
        for (int i = 0 ; i < 4 ; i++)
        	if (_pointMaxCoFacet[cellLID][i] >= 0) sum++;
        // return the value, how many cells has this point, on the leaf level
        rnt = sum;
    }
    else if (cellDim==1) { // 2 dimensions
        int LID = _Edgeindex[cellLID];
		if ((_edgeMaxCoFacet[LID][0] >= 0) && (_edgeMaxCoFacet[LID][1] >= 0))
			rnt = 2;
		else
			rnt = 1;
    }
	//printf(" RET = %d \n" ,rnt );
	return rnt;
}


int HNodeMesh2D::maxCofacetLID(int cellDim, int cellLID,
                       int cofacetIndex,
                       int& facetIndex) const  {

	//printf("HNodeMesh2D::maxCofacetLID() cellDim:%d,  cellLID:%d, cofacetIndex:%d , facetIndex:%d\n",
	//		  cellDim,  cellLID, cofacetIndex , facetIndex);
	int rnt =-1;
	if (cellDim==0) { // 1 dimension
		//facetIndex = cofacetIndex;
		int actCoFacetIndex = 0;
        int LID = _Cellindex[cellLID];
		for (int ii = 0 ; ii < 4 ; ii++){
			if ( _pointMaxCoFacet[LID][ii] >= 0 ){
				if ( actCoFacetIndex < cofacetIndex ){
					actCoFacetIndex++;
				}else{
					facetIndex = ii;
					rnt = _pointMaxCoFacet[LID][ii];
					break;
				}
			}
		}
    }
    else if (cellDim==1) { // 2 dimensions
    	int orientation = 0;
    	int addFakt = 0;
    	int maxCoFacet = 0;
        int LID = _Edgeindex[cellLID];
		double tmp_val = (_point_x_coords[_edgePoints[LID][1]] - _point_x_coords[_edgePoints[LID][0]]);

		if (fabs(tmp_val) > 0.00001 ){
			orientation = 0;   addFakt = 2;
		}else{
			orientation = 1;   addFakt = 1;
		}
		//printf("HNodeMesh2D::maxCofacetLID() , orientation:%d , addFakt:%d \n" ,orientation , addFakt );
        // return the index in the vector, which later will be corrected later
		int actCoFacetIndex = 0;
		for (int ii = 0 ; ii < 2 ; ii++){
			// todo: "2-ii" needs to be corrected
			if ( _edgeMaxCoFacet[LID][1-ii] >= 0 ){
				if ( actCoFacetIndex < cofacetIndex ){
					actCoFacetIndex++;
				}else{
					facetIndex = 1-ii;
					maxCoFacet = _edgeMaxCoFacet[LID][1-ii];
					break;
				}
			}
		}
		//printf("HNodeMesh2D::maxCofacetLID() , facetIndex:%d , _edgeMaxCoFacet[0]:%d ,"
		//		"_edgeMaxCoFacet[1]:%d\n" ,facetIndex , _edgeMaxCoFacet[LID][0] , _edgeMaxCoFacet[LID][1]);
		// calculate the correct facetIndex
		if ( orientation == 0 ){
			facetIndex = facetIndex + facetIndex*addFakt; // this should be eighter 0 or 3
		} else {
			facetIndex = facetIndex + addFakt; // this should be eighter 1 or 2
		}
		//printf(" maxCoFacet = %d \n" ,maxCoFacet );
		rnt = ( maxCoFacet );
    }
	// transform back to leaf indexing
	rnt = _CellReindex[rnt];
	//printf(" RET = %d  facetIndex:%d\n" ,rnt ,facetIndex);
	return rnt;
}

void HNodeMesh2D::getCofacets(int cellDim, int cellLID,
                 int cofacetDim, Array<int>& cofacetLIDs) const {
    TEST_FOR_EXCEPTION(true, InternalError," HNodeMesh2D::getCofacets() not implemented yet");
}


void HNodeMesh2D::getMaxCofacetLIDs(const Array<int>& cellLIDs,
  MaximalCofacetBatch& cofacets) const {
     TEST_FOR_EXCEPTION(true, InternalError," HNodeMesh2D::getMaxCofacetLIDs() not implemented yet");
}


int HNodeMesh2D::mapGIDToLID(int cellDim, int globalIndex) const  {
	SUNDANCE_VERB_HIGH("mapGIDToLID()");
	return globalIndex;
}


bool HNodeMesh2D::hasGID(int cellDim, int globalIndex) const {
	SUNDANCE_VERB_HIGH("hasGID()");
	return true;
}


int HNodeMesh2D::mapLIDToGID(int cellDim, int localIndex) const  {
	SUNDANCE_VERB_HIGH("mapLIDToGID()");
	return localIndex;
}


CellType HNodeMesh2D::cellType(int cellDim) const  {
	 switch(cellDim)
	  {
	    case 0:  return PointCell;
	    case 1:  return LineCell;
	    case 2:  return QuadCell;
	    case 3:  return BrickCell;
	    default:
	      return NullCell; // -Wall
	  }
}


int HNodeMesh2D::label(int cellDim, int cellLID) const {
   TEST_FOR_EXCEPTION(true, InternalError," HNodeMesh2D::label() not implemented yet");
   return 0;
}


void HNodeMesh2D::getLabels(int cellDim, const Array<int>& cellLID,
		Array<int>& labels) const {
   TEST_FOR_EXCEPTION(true, InternalError," HNodeMesh2D::getLabels() not implemented yet");
}

Set<int> HNodeMesh2D::getAllLabelsForDimension(int cellDim) const {
   Set<int>                 rtn;
   TEST_FOR_EXCEPTION(true, InternalError," HNodeMesh2D::getAllLabelsForDimension() not implemented yet");
   return rtn;
}

void HNodeMesh2D::getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const {
   TEST_FOR_EXCEPTION(true, InternalError," HNodeMesh2D::getLIDsForLabel() not implemented yet");
}

void HNodeMesh2D::setLabel(int cellDim, int cellLID, int label) {
   TEST_FOR_EXCEPTION(true, InternalError," HNodeMesh2D::setLabel() not implemented yet");
}


void HNodeMesh2D::assignIntermediateCellGIDs(int cellDim) {
	SUNDANCE_VERB_HIGH("assignIntermediateCellGIDs()");
}


bool HNodeMesh2D::hasIntermediateGIDs(int dim) const {
	SUNDANCE_VERB_HIGH("hasIntermediateGIDs()");
	return true; // true means they have been synchronized ... not used now
}

// =============================== HANGING NODE FUNCTIONS ==========================
bool HNodeMesh2D::isElementHangingNode(int cellDim , int cellLID) const {
	//printf("HNodeMesh2D::isElementHangingNode  cellDim:%d LID:%d \n",cellDim, cellLID);
	if (cellDim==0) { // 1 dimension
        return (_pointIsHanging[cellLID] == 1);
    }
    else if (cellDim==1) { // 2 dimensions
    	int LID = _Edgeindex[cellLID];
        return (_edgeIsHanging[LID] == 1);
    } else {
	return false; //Wall
    }
}

int HNodeMesh2D::indexInParent(int maxCellLID) const {
	return _IndexInParentCell[maxCellLID];
}

void HNodeMesh2D::returnParentFacets(  int childCellLID , int dimFacets ,
		                         Array<int> &facetsLIDs , int &parentCellLIDs ) const {
	parentCellLIDs = _parentCellIndex[_Cellindex[childCellLID]];
	facetsLIDs.resize(4);
	//printf("HNodeMesh2D::returnParentFacets  childCellLID:%d dimFacets:%d  parentCellLIDs:%d\n",childCellLID, dimFacets , parentCellLIDs);
	// this is the same for edges and for points
	facetsLIDs[0] = facetLID_tree( 2 , parentCellLIDs ,  dimFacets , 0 );
	facetsLIDs[1] = facetLID_tree( 2 , parentCellLIDs ,  dimFacets , 1 );
	facetsLIDs[2] = facetLID_tree( 2 , parentCellLIDs ,  dimFacets , 2 );
	facetsLIDs[3] = facetLID_tree( 2 , parentCellLIDs ,  dimFacets , 3 );
}

/** only used in determining the parents*/
int HNodeMesh2D::facetLID_tree(int cellDim, int cellLID,
                     int facetDim, int facetIndex) const{
    int rnt = -1;
	if (facetDim == 0){ // return the Number/ID of a Vertex
		     rnt = _cellPoints[cellLID][facetIndex];
	} else if (facetDim == 1){
	      if (_isEdgeLeaf[_cellEdges[cellLID][facetIndex]] == 1 ){
	    	 rnt = _cellEdges[cellLID][facetIndex];
	    	 rnt = _EdgeReindex[rnt];
	      }
	}
	//printf("HNodeMesh2D::facetLID_tree cellDim:%d, cellLID:%d, facetDim:%d, facetIndex:%d RET = %d \n"
	//		, cellDim , cellLID , facetDim , facetIndex , rnt );
	return rnt;
}
