/* @HEADER@ */
/* @HEADER@ */


#ifndef SUNDANCE_FIATEXPANSION_H
#define SUNDANCE_FIATEXPANSION_H

#include "SundanceDefs.hpp"
#include "SundancePoint.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceExceptions.hpp"



namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;

  namespace Internal
  {
    double eval_jacobi( double a, double b, int n, double x );
  
    double psitilde_a( int i , double z );

    double psitilde_b( int i , int j , double z ) ;

    Point xi_to_eta( const Point& xi );

    void phi( int p , int q , const Array<Point>& pts ,
              Array<double>& result );

    void phis( int degree , const Array<Point>& pts ,
               Array<Array<double> >& result );

    void matmul( const Array<Array<double> >& A ,
                 const Array<Array<double> >& B ,
                 Array<Array<double> >& C ) ;

    void matcopy( const Array<Array<double> >& Mat1, 
                  Array<Array<double> >& Mat2 );

    /* assumes vals is column-major (rows are contiguous memory) */
    /* I'll use this to create the correct matrices from FIAT dums */
    void doublesIntoArray( int m, int n, const double* vals, 
                           Array<Array<double> >& Mat );
  }

}



#endif

