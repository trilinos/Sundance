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

