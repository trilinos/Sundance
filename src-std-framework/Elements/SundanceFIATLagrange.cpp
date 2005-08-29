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
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

//#ifdef BLARF


static CellType sdim_to_cellType[] = { PointCell ,
                                       LineCell ,
                                       TriangleCell ,
                                       TetCell };

static int dim_to_facet_id[] = {0,2,0,0};

// spatialDim is the topological dimension of the mesh
// cellType is the kind of facet on which we're counting the
// degrees of freedom.  We count the degrees of freedom
// associated with that facet and the lower-dimensional facets
// covering it.  For Lagrange elements, this happens to be
// the same number of nodes as belong to the Lagrange element
// of the lower topological dimension.

FIATLagrange::FIATLagrange( int order ) :
  order_( order )
{
  fiatElems_[ LineCell ] = new FIAT::LagrangeElement( 1 , order );
  fiatElems_[ TriangleCell ] = new FIAT::LagrangeElement( 2 , order );
  fiatElems_[ TetCell ] = new FIAT::LagrangeElement( 3 , order );
}

int FIATLagrange::nNodes( int spatialDim ,
                          const CellType& cellType ) const
{
  if (PointCell == cellType) return 1;
  else return fiatElems_[ cellType ]->getPolynomialSet()->size();
}

void FIATLagrange::getLocalDOFs( const CellType &cellType ,
                                 Array<Array<Array<int> > > &dofs ) const
{
  if (PointCell == cellType) {
    dofs.resize(1);
    dofs[0] = tuple(tuple(0));
  }
  else {
    vector<vector<vector<int> > > eids = 
      fiatElems_[ cellType ]->getDualBasis()->getEntityIds();
    dofs.resize(eids.size());
    for (int d=0;d<eids.size();d++) {
      dofs[d].resize(eids[d].size());
      for (int e=0;e<eids[d].size();e++) {
        dofs[d][e].resize( eids[d][e].size() );
        for (int n=0;n<eids[d][e].size();n++) {
          dofs[d][e][n] = eids[d][e][n];
        }
      }
    }
  }

  return;
}

void FIATLagrange::refEval( int spatialDim ,
                            const CellType& cellType ,
                            const Array<Point>& pts ,
                            const MultiIndex& deriv ,
                            Array<Array<double> >& result ) const
{
  int nbf = nNodes( spatialDim , cellType );
  result.resize(pts.length());
  for (int i=0;i<pts.length();i++) {
    result[i].resize(nbf);
  }

  if (PointCell == cellType) {
    result = tuple(tuple(1.0));
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

    
    std::cout << "points in fiat" << endl << fiat_pts << endl;

    // convert deriv multi index into array for fiat.
    blitz::Array<int,1> alpha(spatialDim);
    for (int i=0;i<spatialDim;i++) {
      alpha(i) = deriv[i];
    }

    // chain rule for coordinate change;
    double factor = pow( 0.5 , deriv.order() );

    std::cout << "derivative multiindex in FIAT:" << alpha << endl;

    CellType ct = sdim_to_cellType[spatialDim];

    FIAT::RCP<FIAT::ScalarPolynomialSet> phis
      = fiatElems_[ct]->getPolynomialSet()->multi_deriv_all( alpha );

    blitz::Array<double,2> fiat_results 
      = phis->tabulate( fiat_pts );
    
    std::cout << "tabulated FIAT results" << endl << fiat_results << endl;

    const vector<vector<vector<int> > >&eids 
      = fiatElems_[ct]->getDualBasis()->getEntityIds();

    // extract pertinent basis functions.
    // ones from vertex 0 and 1 go in, plus from edge 2
    // and maybe triangle 0 and tet 0.


    if (spatialDim == 1) {
      int cur = 0;

      // get vertex dof

      for (int v=0;v<2;v++) {
        for (int n=0;n<eids[0][v].size();n++) {
          for (int p=0;p<pts.length();p++) {
            result[p][cur] = factor * fiat_results( eids[0][v][n] , p );
          }
          cur++;
        }
      }
      // get internal dof
      for (int e=0;e<1;e++) {
        for (int n=0;n<eids[1][e].size();n++) {
          for (int p=0;p<pts.length();p++) {
            result[p][cur] = factor * fiat_results( eids[1][e][n] , p );
          }
          cur++;
        }
      }
    }
    else if (spatialDim == 2) {
      if (LineCell==cellType) {   // line of a triangle
        int cur = 0;
        for (int v=0;v<2;v++) {  // grab two vertices
          for (int n=0;n<eids[0][v].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[0][v][n] , p );
            }
            cur++;
          }
        }
        // grab line 2
        for (int e=2;e<3;e++) {
          for (int n=0;n<eids[1][e].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[1][e][n] , p );
            }
            cur++;
          }
        } 
      }
      else { // triangle
        int cur = 0;
        for (int v=0;v<3;v++) {  // grab vertices
          for (int n=0;n<eids[0][v].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[0][v][n] , p );
            }
            cur++;
          }
        }
        for (int e=0;e<3;e++) { // Lines swizzled to match Sundance
          for (int n=0;n<eids[1][e].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[1][(e+2)%3][n] , p );
            }
            cur++;
          }
        }
        for (int f=0;f<1;f++) { // interior
          for (int n=0;n<eids[2][f].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[2][f][n] , p );
            }
            cur++;
          }
        }
      }
    }
    else if (spatialDim == 3) {
      if (LineCell==cellType) {   // line of a tet
        int cur = 0;
        for (int v=0;v<2;v++) {  // grab two vertices
          for (int n=0;n<eids[0][v].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[0][v][n] , p );
            }
            cur++;
          }
        }
        // grab line 2
        for (int e=2;e<3;e++) {
          for (int n=0;n<eids[1][e].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[1][e][n] , p );
            }
            cur++;
          }
        } 
      }
      else if (TriangleCell==cellType) {
        int cur = 0;
        for (int v=0;v<2;v++) {  // grab two vertices
          for (int n=0;n<eids[0][v].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[0][v][n] , p );
            }
            cur++;
          }
        }
        // grab line 2
        for (int e=2;e<3;e++) {
          for (int n=0;n<eids[1][e].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[1][e][n] , p );
            }
            cur++;
          }
        }
        // grab face 0
        for (int f=0;f<1;f++) {
          for (int n=0;n<eids[2][0].size();n++) {
            for (int p=0;p<pts.length();p++) {
              result[p][cur] = factor * fiat_results( eids[2][f][n] , p );
            }
            cur++;
          }
        }

      }
      else {
        for (int i=0;i<pts.length();i++) {
          for (int j=0;j<nbf;j++) { 
            result[i][j] = factor * fiat_results( j , i );
          }
        }
      }
      
    }
    return;
  }
}


void FIATLagrange::print(ostream& os) const
{
  os << "FIATLagrange(" << order() << ")";
}

//#endif
