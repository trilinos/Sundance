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


#include "SundanceEvaluationTester.hpp"
#include "SundanceOut.hpp"
#include "SundanceExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceMultiIndex.hpp"

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceTesting;
using namespace TSFExtended;
using namespace Teuchos;
using namespace std;



EvaluationTester::EvaluationTester(const Expr& e)
  : e_(e),
    rqc_(),
    context_(),
    mgr_(),
    mediator_(),
    tem_(),
    ev_(dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get())),
    sparsity_(),
    unkIDToDiscreteIDMap_()
{
  verbosity() = classVerbosity();
  Tabs tabs;
  int maxDiffOrder = 2;

  SUNDANCE_VERB_LOW(tabs << "creating tester for expression " << e.toString());

  rqc_ = RegionQuadCombo(rcp(new CellFilterStub()), 
                         rcp(new QuadratureFamilyStub(1)));


  


  Array<Expr> unks;
  RefCountPtr<Set<int> > unkID = rcp(new Set<int>());

  SUNDANCE_VERB_LOW(tabs << "finding unknowns...");

  const UnknownFuncElement* ue 
    = dynamic_cast<const UnknownFuncElement*>(e[0].ptr().get());
  if (ue!=0)
    {
      unkID->put(ue->funcID());
      unks.append(e[0]);
    }
  else
    {
      ev_->getUnknowns(*unkID, unks);
    }

  context_ = EvalContext(rqc_, maxDiffOrder, EvalContext::nextID());  
  mgr_.setRegion(context_);


  Array<Expr> u0;

  for (unsigned int i=0; i<unks.size(); i++)
    {
      Tabs tabs1;
      const UnknownFuncElement* fe 
        = dynamic_cast<const UnknownFuncElement*>(unks[i].ptr().get());
      TEST_FOR_EXCEPTION(fe==0, InternalError,
                         "unk " << unks[i] << " is not an UnknownFunction");
      const TestUnknownFunction* tu 
        = dynamic_cast<const TestUnknownFunction*>(fe->master());

      TEST_FOR_EXCEPTION(tu==0, InternalError,
                         "unk " << unks[i] << " is not a TestUnknownFunction");

      Expr discFunc = tu->createDiscreteFunction();

      const DiscreteFuncElement* df 
        = dynamic_cast<const DiscreteFuncElement*>(discFunc[0].ptr().get());
      TEST_FOR_EXCEPTION(df==0, InternalError,
                         "df " << discFunc 
                         << " is not a DiscreteFuncElement");

      SUNDANCE_VERB_LOW(tabs1 << i << " unk=" << unks[i] 
                        << " disc=" << discFunc);

      u0.append(discFunc);
      unkIDToDiscreteIDMap_.put(fe->funcID(), df->funcID());
    }

  Expr unkList = new ListExpr(unks);
  Expr discList = new ListExpr(u0);

  SUNDANCE_VERB_LOW(tabs << "creating test eval mediator...");
  
  tem_ = new TestEvalMediator(discList);
  mediator_ = rcp(tem_);

  mgr_.setMediator(mediator_);

  SUNDANCE_VERB_LOW(tabs << "setting up evaluation...");

  DerivSet d = SymbPreprocessor::setupExpr(e[0], 
                                           unkList,
                                           discList,
                                           context_);

  sparsity_ = ev_->sparsitySuperset(context_);
}



double EvaluationTester::evaluate(Array<double>& firstDerivs, 
                                  Array<Array<double> >& secondDerivs) const 
{
  Tabs tabs;

  firstDerivs.resize(tem_->numFields());
  secondDerivs.resize(tem_->numFields());
  for (int i=0; i<tem_->numFields(); i++) 
    {
      firstDerivs[i] = 0.0;
      secondDerivs[i].resize(tem_->numFields());
      for (int j=0; j<tem_->numFields(); j++)
        {
          secondDerivs[i][j] = 0.0;
        }
    }

  Array<double> constantResults;
  Array<RefCountPtr<EvalVector> > vectorResults;

  
  SUNDANCE_VERB_MEDIUM(tabs << endl << tabs << "evaluating...");

  tem_->setEvalPoint(ADField::evalPoint());
  ev_->evaluator(context_)->resetNumCalls();
  ev_->evaluate(mgr_, constantResults, vectorResults);

  if (verbosity() > VerbLow)
    {
      ev_->sparsitySuperset(context_)->print(cerr, vectorResults,
                                             constantResults);
    }

  int vectorCount=0;
  int constantCount=0;
  double rtn = 0.0;

  for (int i=0; i<sparsity_->numDerivs(); i++)
    {
      const MultipleDeriv& md = sparsity_->deriv(i);
      
      Array<int> fieldIndex;
      Array<MultiIndex> mi;

      SUNDANCE_VERB_MEDIUM("md=" << md);

      for (MultipleDeriv::const_iterator 
             iter=md.begin(); iter != md.end(); iter++)
        {
          const Deriv& d = *iter;
          TEST_FOR_EXCEPTION(d.isCoordDeriv(), InternalError,
                             "coordinate deriv found in TestEvalMediator::"
                             "sumFunctionalChainRule");
          const FunctionalDeriv* f = d.funcDeriv();
          int uid = f->funcID();
          SUNDANCE_VERB_EXTREME("deriv=" << d << " uid=" << uid);
          TEST_FOR_EXCEPTION(!unkIDToDiscreteIDMap_.containsKey(uid),
                             InternalError,
                             "uid " << uid << " not found in map " 
                             <<unkIDToDiscreteIDMap_ );
          int fid = unkIDToDiscreteIDMap_.get(uid);
          int m = tem_->funcIdToFieldNumberMap().get(fid);
          fieldIndex.append(m);
          mi.append(f->multiIndex());
        }

      SUNDANCE_VERB_MEDIUM("field indices are " << fieldIndex
                           << ", multiindices are=" <<mi.toString() );


      int resultIndex;
      if (sparsity_->state(i)==ConstantDeriv)
        {
          resultIndex = constantCount++;
        }
      else
        {
          resultIndex = vectorCount++;
        }
      
      
      double fieldDerivs = 1.0;
      for (unsigned int k=0; k<fieldIndex.size(); k++)
        {
          fieldDerivs *= tem_->evalDummyBasis(fieldIndex[k], mi[k]);
        }

      double coeff;
      if (sparsity_->state(i)==ConstantDeriv)
        {
          coeff = constantResults[resultIndex];
        }
      else
        {
          coeff = vectorResults[resultIndex]->start()[0];
        }
      SUNDANCE_VERB_HIGH("field deriv " << md.toString() 
                           << " value=" << fieldDerivs << " coeff=" << coeff);

      if (fieldIndex.size() == 0)
        {
          SUNDANCE_VERB_HIGH("adding " << coeff * fieldDerivs << " to result");
          rtn += coeff * fieldDerivs;
        }
      else if (fieldIndex.size() == 1)
        {
          SUNDANCE_VERB_HIGH("adding " << coeff * fieldDerivs << " to first deriv");
          firstDerivs[fieldIndex[0]] += coeff * fieldDerivs;
        }
      else
        {
          int multiplicity = 1;
          if (fieldIndex[0] != fieldIndex[1] || !(mi[0] == mi[1]))
            {
              multiplicity = 2;
            }
          SUNDANCE_VERB_HIGH("adding " << coeff * fieldDerivs << " to second deriv with multiplicity " << multiplicity);
          secondDerivs[fieldIndex[0]][fieldIndex[1]] += multiplicity * coeff * fieldDerivs;
      //     if (fieldIndex[0] != fieldIndex[1])
//             {
//               secondDerivs[fieldIndex[1]][fieldIndex[0]] += multiplicity * coeff * fieldDerivs;
//             }

        }
    }

  Array<Array<double> > tmp = secondDerivs;
  /* symmetrize the second derivs */
  for (int i=0; i<tem_->numFields(); i++) 
    {
      for (int j=0; j<tem_->numFields(); j++)
        {
          secondDerivs[i][j] = (tmp[i][j] + tmp[j][i])/2.0;
        }
    }
  return rtn;
}


double EvaluationTester
::fdEvaluate(const double& step, const double& tol, const double& tol2,
             bool& isOK)
{
  Array<double> afdFirst(tem_->numFields());
  Array<Array<double> > afdSecond(tem_->numFields()); 

  Array<double> fdFirst(tem_->numFields());
  Array<Array<double> > fdSecond(tem_->numFields());

  Array<double> firstPlus(tem_->numFields());
  Array<double> firstMinus(tem_->numFields());
  Array<Array<double> > tmpSecond(tem_->numFields());
  
  for (int i=0; i<tem_->numFields(); i++) 
    {
      afdSecond[i].resize(tem_->numFields());
      fdSecond[i].resize(tem_->numFields());
      tmpSecond[i].resize(tem_->numFields());
    }

  cerr << setprecision(14);
  double f0 = evaluate(afdFirst, afdSecond);

  for (int i=0; i<tem_->numFields(); i++)
    {
      Tabs tab2;
      double A0 = tem_->fieldCoeff(i);
      tem_->setFieldCoeff(i, A0 - step);
      double fMinus = evaluate(firstMinus, tmpSecond);
      tem_->setFieldCoeff(i, A0 + step);
      double fPlus = evaluate(firstPlus, tmpSecond);
      fdFirst[i] = (fPlus - fMinus)/2.0/step;
      tem_->setFieldCoeff(i, A0);
      double error1 = fabs(fdFirst[i] - afdFirst[i])/(step + fabs(afdFirst[i]));
      cerr << tab2 << endl;
      SUNDANCE_VERB_MEDIUM(tab2 << "f(A_i+h)=" << fPlus << "      f(A_i-h)=" << fMinus);

      cerr << tab2 << i << " exact first deriv=" << afdFirst[i]
           << "    fd=" << fdFirst[i] << "    |exact - fd|=" 
           << error1 << endl;
      if (error1 > tol) 
        {
          isOK = false;
          cerr << tab2 << "first deriv wrt field "
               << i << " calculation FAILED" << endl;
        }
      cerr << tab2 << endl;

      /* second deriv wrt this field by finite differences 
      * on the AFD first derivs*/
      fdSecond[i][i] = (firstPlus[i] - firstMinus[i])/2.0/step;
      double error2 = fabs(fdSecond[i][i] - afdSecond[i][i])/(step + fabs(afdSecond[i][i]));
      cerr << tab2 << i << " exact second deriv=" << afdSecond[i][i]
           << "    fd=" << fdSecond[i][i] << "    |exact - fd|=" 
           << error2 << endl;
      if (error2 > tol2) 
        {
          isOK = false;
          cerr << tab2 << "second deriv calculation wrt field "
               << i << " FAILED" << endl;
          Tabs tab3;
          cerr << tab3 << "f'(A_i+h) = " << firstPlus[i] << endl;
          cerr << tab3 << "f'(A_i-h) = " << firstMinus[i] << endl;
          cerr << tab3 << "step=" << step << endl;
        }
      
      /* mixed partials by finite differences on the AFD first derivs*/
      for (int j=0; j<i; j++)
        {
         //  /* -1,-1 node */
//           double B0 = tem_->fieldCoeff(j);
//           tem_->setFieldCoeff(i, A0 - step);
//           tem_->setFieldCoeff(j, B0 - step);
//           double fMM = evaluate(tmpFirst, tmpSecond);
//           /* -1,+1 node */
//           tem_->setFieldCoeff(i, A0 - step);
//           tem_->setFieldCoeff(j, B0 + step);
//           double fMP = evaluate(tmpFirst, tmpSecond);
//           /* +1,+1 node */
//           tem_->setFieldCoeff(i, A0 + step);
//           tem_->setFieldCoeff(j, B0 + step);
//           double fPP = evaluate(tmpFirst, tmpSecond);
//           /* +1,-1 node */
//           tem_->setFieldCoeff(i, A0 + step);
//           tem_->setFieldCoeff(j, B0 - step);
//           double fPM = evaluate(tmpFirst, tmpSecond);
//           fdSecond[i][j] = (fPP + fMM - fPM - fMP)/4.0/step/step;
//           fdSecond[j][i] = fdSecond[i][j];
          fdSecond[i][j] = (firstPlus[j] - firstMinus[j])/2.0/step;
          error2 = fabs(fdSecond[i][j] - afdSecond[i][j])/(step + fabs(afdSecond[i][j]));
          cerr << tab2 << "(" << i << ", " << j 
               << ") exact mixed deriv=" << afdSecond[i][j]
               << "    fd=" << fdSecond[i][j] << "    |exact - fd|=" 
           << error2 << endl;
          if (error2 > tol2) 
            {
              isOK = false;
              cerr << tab2 << "mixed partial deriv calculation wrt fields "
                   << i << ", " << j << " FAILED" << endl;
              Tabs tab3;
              cerr << tab3 << "f'(A_i+h) = " << firstPlus[j] << endl;
              cerr << tab3 << "f'(A_i-h) = " << firstMinus[j] << endl;
              cerr << tab3 << "step=" << step << endl;
            }
          tem_->setFieldCoeff(i, A0);
          //   tem_->setFieldCoeff(j, B0);
        }
    }

  

  return f0;
}


