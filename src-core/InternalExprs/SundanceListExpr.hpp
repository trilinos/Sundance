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

#ifndef SUNDANCE_LISTEXPR_H
#define SUNDANCE_LISTEXPR_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "Teuchos_Array.hpp"

namespace SundanceCore
{
using namespace SundanceUtils;
using namespace Teuchos;
using std::string;


/** */
class ListExpr : virtual public ExprBase
{
public:
  /** */
  ListExpr();

  /** */
  ListExpr(const Array<Expr>& elements);

  /** */
  virtual ~ListExpr() {;}

  /** */
  const Expr& element(int i) const {return elements_[i];}

  /** */
  void append(const Expr& expr);

  /** */
  Expr flatten() const ;

  /** */
  Expr join(const Expr& other) const ;

  /** */
  unsigned int size() const ;

  /** */
  unsigned int totalSize() const ;

  /** Write a simple text description suitable 
   * for output to a terminal */
  virtual ostream& toText(ostream& os, bool paren) const ;

  /** Write in a form suitable for LaTeX formatting */
  virtual ostream& toLatex(ostream& os, bool paren) const ;


  /** Write in XML */
  virtual XMLObject toXML() const ;

  /** */
  virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

private:
  Array<Expr> elements_;
};
}

#endif
