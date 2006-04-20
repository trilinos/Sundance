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

#include "SundanceFunctionalEvaluator.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;



namespace SundanceStdFwk
{
  double evaluateIntegral(const Mesh& mesh, const Expr& expr)
  {
    FunctionalEvaluator eval(mesh, expr);
    return eval.evaluate();
  }
}

FunctionalEvaluator::FunctionalEvaluator()
  : assembler_(),
    varValues_(),
    vecType_(),
    gradient_()
{}

FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh,
                                         const Expr& integral)
  : assembler_(),
    varValues_(),
    vecType_(),
    gradient_()
{
  Array<Expr> fields;
  Expr bcs;
  
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(integral, bcs, fields, fields));
  
  
  assembler_ = rcp(new Assembler(mesh, eqnSet));
}


FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh,
                                         const Expr& integral,
                                         const Expr& bcs,
                                         const Expr& var,
                                         const Expr& varValues,
                                         const VectorType<double>& vectorType)
  : assembler_(),
    varValues_(varValues),
    vecType_(vectorType),
    gradient_()
{
  Array<Expr> v = tuple(var.flatten());
  Array<Expr> v0 = tuple(varValues.flatten());
  Array<Expr> fixed;
  
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(integral, bcs, v, v0, fixed, fixed));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vectorType), tuple(vectorType)));
}


FunctionalEvaluator::FunctionalEvaluator(const Mesh& mesh,
                                         const Expr& integral,
                                         const Expr& bcs,
                                         const Expr& vars,
                                         const Expr& varEvalPts,
                                         const Expr& fields,
                                         const Expr& fieldValues,
                                         const VectorType<double>& vectorType)
  : assembler_(),
    varValues_(varEvalPts),
    vecType_(vectorType),
    gradient_()
{
  Array<Expr> f = tuple(fields.flatten());
  Array<Expr> f0 = tuple(fieldValues.flatten());
  Array<Expr> v = tuple(vars.flatten());
  Array<Expr> v0 = tuple(varEvalPts.flatten());
  
  RefCountPtr<EquationSet> eqnSet 
    = rcp(new EquationSet(integral, bcs, v, v0, f, f0));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vectorType), tuple(vectorType)));
}

double FunctionalEvaluator::evaluate() const
{
  double value;
  assembler_->evaluate(value);
  return value;
}


Vector<double> FunctionalEvaluator::evalGradientVector(double& value) const 
{
  assembler_->evaluate(value, gradient_);

  return gradient_;
}

Expr FunctionalEvaluator::evalGradient(double& value) const 
{
  Vector<double> g = evalGradientVector(value);

  Array<Expr> rtn(assembler_->rowSpace().size());
  for (unsigned int i=0; i<rtn.size(); i++)
    {
      string name = "gradient";
      if (rtn.size() > 1) name += "[" + Teuchos::toString(i) + "]";
      rtn[i] = new DiscreteFunction(*(assembler_->rowSpace()[i]),
                                    g.getBlock(i), name);
    }
  if ((int)rtn.size()==1)
    {
      return rtn[0];
    }
  else
    {
      return new ListExpr(rtn);
    }
}


double FunctionalEvaluator::fdGradientCheck(double h) const
{
  bool showAll = false;
  Tabs tabs;
  double f0, fPlus, fMinus;
  Expr gradF0 = evalGradient(f0);
  


  DiscreteFunction* df = DiscreteFunction::discFunc(varValues_);
  DiscreteFunction* dg = DiscreteFunction::discFunc(gradF0);
  Vector<double> x = df->getVector();

  RefCountPtr<GhostView<double> > xView = df->ghostView();
  RefCountPtr<GhostView<double> > gradF0View = dg->ghostView();


  TEST_FOR_EXCEPTION(xView.get() == 0, RuntimeError, 
                     "bad pointer in FunctionalEvaluator::fdGradientCheck");
  TEST_FOR_EXCEPTION(gradF0View.get() == 0, RuntimeError, 
                     "bad pointer in FunctionalEvaluator::fdGradientCheck");

  int nTot = x.space().dim();
  int n = x.space().numLocalElements();
  int lowestIndex = x.space().lowestLocallyOwnedIndex();

  cerr << tabs << "doing fd check:  h=" << h << endl;
  Array<double> df_dx(n);

  int localIndex = 0;
  for (int globalIndex=0; globalIndex<nTot; globalIndex++)
    {
      double tmp=0.0;
      bool isLocal = globalIndex >= lowestIndex 
        && globalIndex < (lowestIndex+n);
      if (isLocal)
        {
          tmp = xView->getElement(globalIndex);
          x.setElement(globalIndex, tmp + h);
        }
      df->setVector(x);
      fPlus = evaluate();
      if (isLocal)
        {
          x.setElement(globalIndex, tmp - h);
        }
      df->setVector(x);
      fMinus = evaluate();
      
      if (isLocal)
        {
          df_dx[localIndex] = (fPlus - fMinus)/2.0/h;
          if (showAll)
            {
              cerr << "i " << globalIndex << " x_i=" << tmp 
                   << " f(x)=" << f0 
                   << " f(x+h)=" << fPlus 
                   << " f(x-h)=" << fMinus << endl;
            }
          x.setElement(globalIndex, tmp);
          localIndex++;
        }
      df->setVector(x);
    }
  
  double localMaxErr = 0.0;

  showAll = true;
  for (int i=0; i<n; i++)
    {
      int globalIndex = lowestIndex + i;
      double num =  fabs(df_dx[i]-gradF0View->getElement(globalIndex));
      double den = fabs(df_dx[i]) + fabs(gradF0View->getElement(globalIndex));
      double r = 0.0;
      if (fabs(den) > 1.0e-10) r = num/den;
      else r = 1.0;
      if (showAll)
        {
          cerr << "i " << i << " g=" << globalIndex << " FD=" << df_dx[i] 
               << " grad=" << gradF0View->getElement(globalIndex)
               << " r=" << r << endl;
        }
      if (localMaxErr < r) localMaxErr = r;
    }
  cout << "local max error = " << localMaxErr << endl;
  
  double maxErr = localMaxErr;
  MPIComm::world().allReduce((void*) &localMaxErr, (void*) &maxErr, 1, 
                             MPIComm::DOUBLE, MPIComm::MAX);
  cerr << tabs << "fd check: max error = " << maxErr << endl;

  return maxErr;
}
