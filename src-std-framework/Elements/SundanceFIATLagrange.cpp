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

#include "SundanceFIATLagrange.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace Teuchos;

#ifdef HAVE_FIAT

static int line_sdvert_to_fvert[] = {0,1};
static int line_sdline_to_fline[] = {0};
static int *line_sd_to_fiat[] = {line_sdvert_to_fvert,line_sdline_to_fline};

static int tri_sdvert_to_fvert[] = {0,1,2};
//static int tri_sdline_to_fline[] = {2,0,1};
static int tri_sdline_to_fline[] = {2,0,1};
static int tri_sdtri_to_ftri[] = {0};

static int *tri_sd_to_fiat[] = {tri_sdvert_to_fvert,
				tri_sdline_to_fline,
				tri_sdtri_to_ftri};

static int tet_sdvert_to_fvert[] = {0,1,2,3};
static int tet_sdline_to_fline[] = {2,0,1,3,4,5};
static int tet_sdtri_to_ftri[] = {0,1,2,3};
static int tet_sdtet_to_ftet[] = {0};

static int *tet_sd_to_fiat[] = {tet_sdvert_to_fvert,
				tet_sdline_to_fline,
				tet_sdtri_to_ftri,
				tet_sdtet_to_ftet};

static int **sd_to_fiat[] = {NULL,
			     line_sd_to_fiat,
			     tri_sd_to_fiat,
			     tet_sd_to_fiat};


static CellType sdim_to_cellType[] = { PointCell ,
                                       LineCell ,
                                       TriangleCell ,
                                       TetCell };


// spatialDim is the topological dimension of the mesh
// cellType is the kind of facet on which we're counting the
// degrees of freedom.  We count the degrees of freedom
// associated with that facet and the lower-dimensional facets
// covering it.  For Lagrange elements, this happens to be
// the same number of nodes as belong to the Lagrange element
// of the lower topological dimension.

FIATLagrange::FIATLagrange( int order ) 
  : order_( order )
{
  TEST_FOR_EXCEPTION(order < 0, RuntimeError,
    "invalid polynomial orde r=" << order
    << " in  FIATLagrange ctor");
  fiatElems_[ LineCell ] = new FIAT::LagrangeElement( 1 , order );
  fiatElems_[ TriangleCell ] = new FIAT::LagrangeElement( 2 , order );
  fiatElems_[ TetCell ] = new FIAT::LagrangeElement( 3 , order );
}


bool FIATLagrange::supportsCellTypePair(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maxCellType)
  {
    case LineCell:
      switch(cellType)
      {
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    case TriangleCell:
      switch(cellType)
      {
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    case TetCell:
      switch(cellType)
      {
        case TetCell:
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    default:
      return false;
  }
}


int FIATLagrange::nReferenceDOFs(
  const CellType& maximalCellType, 
  const CellType& cellType
  ) const 
{
  if (PointCell == cellType) return 1;
  else return fiatElems_[ cellType ]->getPolynomialSet()->size();
}

void FIATLagrange::getReferenceDOFs(
  const CellType& maximalCellType, 
  const CellType& cellType,
  Array<Array<Array<int> > > &dofs ) const
{
  if (PointCell == cellType) {
    dofs.resize(1);
    dofs[0] = tuple(tuple(0));
  }
  else {
    //    int sdim = dimension( cellType );
    //    int **transfer = sd_to_fiat[sdim];
    vector<vector<vector<int> > > eids = 
      fiatElems_[ cellType ]->getDualBasis()->getEntityIds();
    dofs.resize(eids.size());
    for (int d=0;d<(int)eids.size();d++) {
      dofs[d].resize(eids[d].size());
      //int *transfer_cur = transfer[d];
      for (int e=0;e<(int)eids[d].size();e++) {
        dofs[d][e].resize( eids[d][e].size() );
        for (int n=0;n<(int)eids[d][e].size();n++) {
          dofs[d][e][n] = eids[d][e][n];
          //dofs[d][transfer_cur[e]][n] = eids[d][e][n];
        }
      }
    }
  }

  return;
}

void FIATLagrange::refEval( 
  const CellType& maximalCellType, 
  const CellType& cellType,
  const Array<Point>& pts ,
  const MultiIndex& deriv ,
  Array<Array<Array<double> > >& result ) const
{
  result.resize(1);
  int nbf = nNodes( spatialDim , cellType );
  result[0].resize(pts.length());
  for (int i=0;i<pts.length();i++) {
    result[0][i].resize(nbf);
  }

  if (PointCell == cellType) {
    result[0] = tuple(tuple(1.0));
    return;
  }
  else {
    // convert points into array for fiat
    // this involves:
    // i) conversion from (0,1) --> (-1,1) reference domains
    // ii) embedding points on cellType into spatialDim (this
    //   is just appending enough -1's to the point after I transform
    // iii) getting it into the right data structure.
    blitz::Array<double,2> fiat_pts( pts.length() , spatialDim );
    int pd = dimension( cellType );
    for (int i=0;i<pts.length();i++) {
      for (int j=0;j<pd;j++) {
        fiat_pts(i,j) = 2.0 * ( pts[i][j] - 0.5 );
      }
      for (int j=pd;j<spatialDim;j++) {
        fiat_pts(i,j) = -1.0;
      }
    }

    // convert deriv multi index into array for fiat.
    blitz::Array<int,1> alpha(spatialDim);
    for (int i=0;i<spatialDim;i++) {
      alpha(i) = deriv[i];
    }

    // chain rule for coordinate change;
    double factor = pow( 2.0 , deriv.order() );

    CellType ct = sdim_to_cellType[spatialDim];

    FIAT::RCP<FIAT::ScalarPolynomialSet> phis
      = fiatElems_[ct]->getPolynomialSet()->multi_deriv_all( alpha );

    blitz::Array<double,2> fiat_results1 = phis->tabulate( fiat_pts );
    blitz::Array<double,2> fiat_results(fiat_results1.rows(),fiat_results1.columns());
    fiat_results = factor * fiat_results1;
    
    //    std::cout << "tabulated FIAT results" << endl << fiat_results << endl;

    const vector<vector<vector<int> > >&eids 
      = fiatElems_[ct]->getDualBasis()->getEntityIds();

    int cur=0;
    int **sd_to_fiat_spd = sd_to_fiat[spatialDim];
    for (int d=0;d<=pd;d++) {
      int *sd_to_fiat_spd_d = sd_to_fiat_spd[d];
      for (int e=0;e<FIAT::Shapes::num_entities(pd,d);e++) {
	int fiat_e = sd_to_fiat_spd_d[e];
	for (int n=0;n<(int)eids[d][e].size();n++) {
	  for (int p=0;p<(int)pts.length();p++) {
	    result[0][p][cur] = fiat_results( eids[d][fiat_e][n] , p );
	  }
	  cur++;
	}
      }
    }
    return;
  }
}


void FIATLagrange::print(std::ostream& os) const
{
  os << "FIATLagrange(" << order() << ")";
}

#endif
