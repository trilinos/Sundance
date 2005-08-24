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
#include "SundanceFIATExpansion.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

#ifdef IFDEF_OUT_OLD_CODE

// <------------ Added array bounds

static double V[][9] = { {0.333333333333 , -0.5 , -0.166666666667 , 0.333333333333 , 0.5 , -0.166666666667 , 0.333333333333 , 0.0 , 0.333333333333} };
static double D[][2][9] = { { {-0.5 , -0.5 , -0.5 , 0.5 , 0.5 , 0.5 , 0.0 , 0.0 , 0.0} } , { {-0.5 , -0.5 , -0.5 , 1.38777878078e-17 , 1.38777878078e-17 , 1.38777878078e-17 , 0.5 , 0.5 , 0.5} } };


// <------- CHANGED V_ and D_ to VDM_ and derivMats_

FIATLagrange::FIATLagrange(int order)
  : ScalarBasis(), order_(order), 
    VDM_((order+1)*(order+2)/2,Array<double>((order+1)*(order+2)/2)), 
    derivMats_(2,Array<Array<double> >((order+1)*(order+2)/2,
                                       Array<double>((order+1)*(order+2)/2)))
{
  int i;
  int dim = (order+1)*(order+2)/2;
  if (order > MAXDEGREE_) {
    SUNDANCE_ERROR( "Maximal degree exceeded" );
  }

  SundanceStdFwk::Internal::doublesIntoArray( dim,dim,V[order-1],VDM_ );
  for (i=0;i<2;i++) {
    SundanceStdFwk::Internal::doublesIntoArray(dim,dim,D[order-1][i],derivMats_[i]);
  }

  static bool first = true;
  if (first)
    {
      cerr << endl 
           << " ----------------- FIAT internal tables ----------------------" 
           << endl;
      cerr << "VDM_ = " << VDM_ << endl << endl;
      cerr << "derivMats_ = " << derivMats_ << endl << endl << endl;
      first = false;
    }
}

void FIATLagrange::print(ostream& os) const 
{
  os << "FIATLagrange(" << order_ << ")";
}

int FIATLagrange::nNodes(int /*spatialDim*/,
                         const CellType& cellType) const
{
  switch(cellType)
    {
    case TriangleCell:
      {
	return (order_+1)*(order_+2)/2;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in FIATLagrange basis");
      return -1; // -Wall
      }
    }
}  // <--------------------- ADDED A BRACE

void FIATLagrange::getLocalDOFs(const CellType& cellType,
                                Array<Array<Array<int> > >& dofs) const 
{
  switch(cellType)
    {
    case TriangleCell:
      {
        int n = order()-1;
        dofs.resize(3);
        dofs[0] = tuple(tuple(0), tuple(1), tuple(2));
        dofs[1] = tuple(makeRange(3,2+n), 
                        makeRange(3+n, 2+2*n),
                        makeRange(3+2*n, 2+3*n));
                                
        dofs[2] = tuple(makeRange(3, order()));
        return;
      }
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "Cell type "
                         << cellType << " not implemented in FIATLagrange basis");
    }
}


Array<int> FIATLagrange::makeRange(int low, int high)
{
  if (high < low) return Array<int>();

  Array<int> rtn(high-low+1);
  for (int i=0; i<rtn.length(); i++) rtn[i] = low+i;
  return rtn;
}

/* result: rows are points, columns are bf */
void FIATLagrange::refEval(int /*spatialDim*/,
                           const CellType& cellType,
                           const Array<Point>& pts,
                           const MultiIndex& deriv,
                           Array<Array<double> >& result) const
{
  result.resize(pts.length());
  int i,j;
  int dim = (order_+1)*(order_+2)/2;
  // <----------- CHANGE FROM REF TO VALUE 
  Array<Array<double> > tmp1(dim,Array<double>(pts.length()));
  Array<Array<double> > tmp2(dim,Array<double>(pts.length()));

  /* Transform input points from KRL's {{0,0}, {1,0}, {0,1}} frame to
  * RCK's {{-1,-1}, {1, -1}, {-1, 1}} frame */
  Array<Point> evalPts(pts.length());
  for (i=0; i<evalPts.length(); i++) 
    {
      evalPts[i] = /* Point(-1.0, -1.0) + 2.0* */pts[i];
    }
    

  switch(cellType)
    {
    case TriangleCell:
      /* evaluate orthogonal basis on triangles */
      SundanceStdFwk::Internal::phis(order_,evalPts,tmp1);
      
      /* convert to nodal basis */
      SundanceStdFwk::Internal::matmul( VDM_ , tmp1 , tmp2 );

      /* apply derivatives */
      for (i=0;i<2*deriv.order();i++) 
        {
          for (j=0;j<2;j++) 
            {
              SundanceStdFwk::Internal::matmul( derivMats_[i] , tmp2 , tmp1 );
              SundanceStdFwk::Internal::matcopy( tmp1 , tmp2 );
            }
      }
      
      /* copy result, transposing */
      for (int k=0; k<pts.length(); k++) result[k].resize(dim);
      for (i=0;i<dim;i++) {
        for (j=0;j<pts.length();j++) {
          result[j][i] = tmp2[i][j];
        }
      }
      return;
    default:
      SUNDANCE_ERROR("FIATLagrange::refEval() unimplemented for cell type ");

    }
}

#endif

