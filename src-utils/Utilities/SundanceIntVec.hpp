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

#ifndef SUNDANCE_INTVEC_H
#define SUNDANCE_INTVEC_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "PlayaPrintable.hpp"

namespace Sundance
{
using Teuchos::Array;
using Teuchos::tuple;

/** 
 * An integer vector class for use in the generalized chain rule
 * of Constantine and Savits (1996). See the paper by CS for 
 * definitions of the various operations.
 */
class IntVec : public Playa::Printable
{
public:
  /** */
  IntVec(){}
  /** */
  IntVec(const Array<int>& d) : data_(d) {}
  /** */
  IntVec(int n);

  /** */
  int size() const {return data_.size();}

  /** */
  int operator[](int i) const {return data_[i];}

  /** */
  int& operator[](int i) {return data_[i];}

  /** */
  IntVec operator+(const IntVec& other) const ;

  /** */
  IntVec operator*(int alpha) const ;

  /** Return the factorial of this vector, as defined by CS */
  int factorial() const ;

  /** */
  int pow(const IntVec& other) const ;

  /** Return the sum of elements in this vector */
  int abs() const ;

  /** Return the infinity norm of this vector */
  int norm() const ;

  /** */
  bool operator==(const IntVec& other) const ;

  /** */
  bool operator<(const IntVec& other) const ;

  /** */
  void print(std::ostream& os) const ;

  /** Get the length-M partitions of this vector. These are all 
   * list of exactly M vectors \f$ v_i\f$ such that
   * \f[ \sum_{i=1}^M v_i = this. \f]
   */
  void getPartitions(int M, Array<Array<IntVec> >& parts) const ;

private:
  Array<int> data_;
};



/** \relates IntVec */
inline IntVec operator*(int a, const IntVec& x)
{
  return x * a;
}

/** \relates IntVec*/
inline std::ostream& operator<<(std::ostream& os, const IntVec& v)
{
  v.print(os);
  return os;
}

/** */
inline IntVec intVec(int a)
{
  Array<int> dat = tuple(a);
  return dat;
}

/** */
inline IntVec intVec(int a, int b)
{
  Array<int> dat = tuple(a, b);
  return dat;
}

/** */
inline IntVec intVec(int a, int b, int c)
{
  Array<int> dat = tuple(a, b, c);
  return dat;
}

/** */
inline IntVec intVec(int a, int b, int c, int d)
{
  Array<int> dat = tuple(a, b, c, d);
  return dat;
}

/** */
inline IntVec intVec(int a, int b, int c, int d, int e)
{
  Array<int> dat = tuple(a, b, c, d, e);
  return dat;
}



}




#endif
