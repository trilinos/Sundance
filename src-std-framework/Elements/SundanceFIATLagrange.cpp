/* @HEADER@ */
/* @HEADER@ */

#include "SundanceFIATLagrange.hpp"
#include "SundanceFIATExpansion.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;


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

  SundanceStdFwk::Internal::doublesIntoArray( dim,dim,V[order],VDM_ );
  for (i=0;i<2;i++) {
    SundanceStdFwk::Internal::doublesIntoArray(dim,dim,D[order][i],derivMats_[i]);
  }
}

void FIATLagrange::print(ostream& os) const 
{
  os << "FIATLagrange(" << order_ << ")";
}

int FIATLagrange::nNodes(const CellType& cellType) const
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
void FIATLagrange::refEval(const CellType& cellType,
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

  switch(cellType)
    {
    case TriangleCell:
      /* evaluate orthogonal basis on triangles */
      SundanceStdFwk::Internal::phis(order_,pts,tmp1);
      
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
      for (int k=0; i<pts.length(); i++) result[i].resize(dim);
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

