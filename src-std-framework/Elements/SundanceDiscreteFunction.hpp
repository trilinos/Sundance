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

#ifndef SUNDANCE_DISCRETEFUNCTION_H
#define SUNDANCE_DISCRETEFUNCTION_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceDiscreteFunctionData.hpp"
#include "SundanceFuncWithBasis.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "PlayaVectorDecl.hpp"

namespace Sundance
{
using namespace Teuchos;
  

/** 
 * DiscreteFunction represents a function that is discretized
 * on a finite-element space.
 */
class DiscreteFunction : public DiscreteFunctionStub,
                         public FuncWithBasis
{
public:
  /** */
  DiscreteFunction(const DiscreteSpace& space, const std::string& name="");

  /** */
  DiscreteFunction(const DiscreteSpace& space, const Vector<double>& vec, 
    const std::string& name="");

  /** */
  DiscreteFunction(const DiscreteSpace& space, const double& constantValue,
    const std::string& name="");
  /** */
  DiscreteFunction(const DiscreteSpace& space, const Array<string>& names);

  /** */
  DiscreteFunction(const DiscreteSpace& space, const Vector<double>& vec, 
    const Array<string>& names);

  /** */
  DiscreteFunction(const DiscreteSpace& space, const double& constantValue,
    const Array<string>& name);

  /** */
  static const DiscreteFunction* discFunc(const Expr& expr);


  /** */
  static DiscreteFunction* discFunc(Expr& expr);

  /** */
  void updateGhosts() const ;

  /** */
  void setVector(const Vector<double>& vec);

  /** */
  const Vector<double>& getVector() const 
    {return data_->getVector();}

  /** */
  const DiscreteSpace& discreteSpace() const 
    {return data_->discreteSpace();}

  /** */
  const Mesh& mesh() const {return discreteSpace().mesh();}

  /** */
  const RCP<DOFMapBase>& map() const {return discreteSpace().map();}


  RCP<GhostView<double> >  ghostView() const 
    {return data_->ghostView();}

  const DiscreteFunctionData* data() const {return data_.get();}


  /** virtual destructor */
  virtual ~DiscreteFunction() {;}

  /* boilerplate */
  GET_RCP(ExprBase);


  /** */
  RCP<const MapStructure> getLocalValues(int cellDim, 
    const Array<int>& cellLID,
    Array<Array<double> >& localValues) const ;


private:
  /** */
  RCP<DiscreteFuncDataStub> getRCP(DiscreteFunctionData* ptr);

  RCP<DiscreteFunctionData> data_;

};


/** \relates DiscreteFunction
 * Replace the vector in oldVals with the vector from newVals.
 */
void updateDiscreteFunction(const Expr& newVals, Expr oldVals);


/** \relates DiscreteFunction
 * Make a copy of the discrete function u0. The copy will have a shallow
 * copy of u0's space, and a deep copy of u0's vector. 
 */
Expr copyDiscreteFunction(const Expr& u0, const string& name = "");


/** \relates DiscreteFunction
 * Add a vector v to the vector underlying the discrete function u.
 */
void addVecToDiscreteFunction(Expr u, const Vector<double>& v);

/** \relates DiscreteFunction
 * Get a shallow copy of the vector underlying a discrete function 
 */
Vector<double> getDiscreteFunctionVector(const Expr& u);


/** \relates DiscreteFunction
 * Set the vector underlying a discrete function 
 */
void setDiscreteFunctionVector(Expr u, const Vector<double>& v);


/** \relates DiscreteFunction
 * Get the mesh underlying a discrete function 
 */
Mesh getDiscreteFunctionMesh(const Expr& u);


}



#endif
