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

#ifndef SUNDANCE_FUNCTIONAL_H
#define SUNDANCE_FUNCTIONAL_H

#include "SundanceDefs.hpp"
#include "SundanceLinearProblem.hpp"
#include "SundanceNonlinearProblem.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "TSFNonlinearOperator.hpp"
#include "TSFLinearSolverDecl.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFVectorType.hpp"

namespace SundanceStdFwk
{
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace Teuchos;

/**
 *
 */
class Functional
  : public SundanceUtils::ParameterControlledObjectWithVerbosity<Functional>
{
public:
  /** */
  Functional(){;}

  /** */
  Functional(
    const Mesh& mesh, 
    const Expr& integral, 
    const TSFExtended::VectorType<double>& vecType,
    const ParameterList& verbParams = *defaultVerbParams());

  /** */
  Functional(
    const Mesh& mesh, 
    const Expr& integral, 
    const Expr& essentialBC,
    const TSFExtended::VectorType<double>& vecType,
    const ParameterList& verbParams = *defaultVerbParams());

  /** */
  LinearProblem linearVariationalProb(const Expr& var,
    const Expr& varEvalPts,
    const Expr& unk,
    const Expr& fixed,
    const Expr& fixedEvalPts) const ;

    
  /** */
  NonlinearProblem
  nonlinearVariationalProb(const Expr& var,
    const Expr& varEvalPts,
    const Expr& unk,
    const Expr& unkEvalPts,
    const Expr& fixed,
    const Expr& fixedEvalPts) const ;


  /** */
  FunctionalEvaluator evaluator(const Expr& var,
    const Expr& varEvalPts,
    const Expr& fixed,
    const Expr& fixedEvalPts) const ;


  /** */
  FunctionalEvaluator evaluator(const Expr& var,
    const Expr& varEvalPts) const ;

  /** */
  const Mesh& mesh() const {return mesh_;}
    
  /** */
  static RefCountPtr<ParameterList> defaultVerbParams();
                                  

private:
  Mesh mesh_;

  Expr integral_;

  Expr bc_;

  TSFExtended::VectorType<double> vecType_;
    
};
}


#endif
