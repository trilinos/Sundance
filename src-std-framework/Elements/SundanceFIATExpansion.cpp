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

#include "SundanceFIATExpansion.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace Teuchos;

namespace SundanceStdFwk
{ 
  namespace Internal
  {

    double eval_jacobi( double a, double b, int n, double x )
    {
      int k;

      if ( 0 == n ) {
        return 1.0;
      }
      else if ( 1 == n ) {
        return 0.5 * ( a - b + ( a + b + 2.0 ) * x );
      }
      else {
        double a1,a2,a3,a4;
        double apb = a + b;
        double pn2 = 1.0;
        double pn1 = 0.5 * ( a - b + ( apb + 2.0 ) * x );
        double p = 0;
        for ( k = 2 ; k <= n ; k++ ) {
          a1 = 2.0 * k * ( k + apb ) * ( 2.0 * k + apb - 2.0 );
          a2 = ( 2.0 * k + apb - 1.0 ) * ( a * a - b * b );
          a3 = ( 2.0 * k + apb - 2.0 )
            * ( 2.0 * k + apb - 1.0 )
            * ( 2.0 * k + apb );
          a4 = 2.0 * ( k + a - 1.0 ) * ( k + b - 1.0 )
            * ( 2.0 * k + apb );
          a2 = a2 / a1;
          a3 = a3 / a1;
          a4 = a4 / a1;
          p = ( a2 + a3 * x ) * pn1 - a4 * pn2;
          pn2 = pn1;
          pn1 = p;
        }
        return p;
      }
    }
  
    double psitilde_a( int i , double z )
    {
      return eval_jacobi( 0 , 0 , i , z );
    }

    double psitilde_b( int i , int j , double z ) 
    {
      return pow( 0.5 * ( 1.0 - z ) , i ) 
        * eval_jacobi( 2.0 * i + 1.0 , 0 , j , z );
    }

    Point xi_to_eta( const Point& xi )
    {
      double eta0;
      if ( xi[1] == 1.0 ) {
        eta0 = -1.0;
      }
      else {
        eta0 = 2.0 * ( 1.0 * xi[0] ) / ( 1.0 - xi[1] ) - 1.0;
      }
      return Point(eta0,xi[1]);
    }

    void phi( int p , int q , const Array<Point>& pts ,
              Array<double>& result )
    {
      int i;
      for (i=0;i<pts.length();i++) {
        Point eta = xi_to_eta( pts[i] );
        result[i] = psitilde_a( p , eta[0] ) * psitilde_b( p , q , eta[1] );
      }
      return;
    }

    void phis( int degree , const Array<Point>& pts ,
               Array<Array<double> >& result ) {
      Array<double> tmp(pts.length());
      int i,k,j,cur=0;
      for (k=0;k<=degree;k++) {
        for (i=0;i<=k;i++) {
          phi(k-i,i,pts,tmp);
          for (j=0;j<pts.length();j++) {
            result[cur][j] = tmp[j];
          }
          cur++;
        }
      }
      return;
    }

    void matmul( const Array<Array<double> >& A ,
                 const Array<Array<double> >& B ,
                 Array<Array<double> >& C ) 
    {
      int i,j,k;
      int n_rows_a = A.length();
      int n_cols_a = A[0].length();
      int n_rows_b = B.length();
      int n_cols_b = B[0].length();
      int n_rows_c = C.length();
      int n_cols_c = C[0].length();
      /* bounds checking */
      //bvbw temporary hack
#ifndef TRILINOS_7
      if (n_cols_a != n_rows_b) {
        SUNDANCE_ERROR("Nonconforming matrices");
      }
      if (n_rows_a != n_rows_c || n_cols_b != n_cols_c ) {
        SUNDANCE_ERROR("Nonconforming matrices");
      }
#else
      if (n_cols_a != n_rows_b) {
        SUNDANCE_ERROR7("Nonconforming matrices");
      }
      if (n_rows_a != n_rows_c || n_cols_b != n_cols_c ) {
        SUNDANCE_ERROR7("Nonconforming matrices");
      }
#endif
      for (i=0;i<n_rows_a;i++) {
        for (j=0;j<n_cols_b;j++) {
          C[i][j] = 0.0;
          for (k=0;k<n_cols_a;k++) {
            C[i][j] += A[i][k] * B[k][j];
          }
        }
      }

      return;
    }

    void matcopy( const Array<Array<double> >& Mat1, 
                  Array<Array<double> >& Mat2 )
    {
      for (int i=0;i<Mat1.length();i++) {
        for (int j=0;j<Mat1[0].length();j++) {
          Mat2[i][j] = Mat1[i][j];
        }
      }
    }

    /* assumes vals is column-major (rows are contiguous memory) */
    /* I'll use this to create the correct matrices from FIAT dums */
    void doublesIntoArray( int m, int n, const double* vals, 
                           Array<Array<double> >& Mat )
    {
      int i,j,cur=0; 
      for (i=0;i<m;i++) { 
        for (j=0;j<n;j++) { 
          Mat[i][j] = vals[cur++];
        } 
      } 
      return;
    } 

  }
}
