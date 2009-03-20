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

#include "SundanceGrouperBase.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceFuncWithBasis.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceTestFunction.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;



void GrouperBase::extractWeakForm(const EquationSet& eqn,
                                  const MultipleDeriv& functionalDeriv,
                                  BasisFamily& varBasis, 
                                  BasisFamily& unkBasis,
                                  MultiIndex& miVar, MultiIndex& miUnk,
                                  int& rawVarID, int& rawUnkID,  
                                  int& reducedVarID, int& reducedUnkID,  
                                  int& testBlock, int& unkBlock, 
                                  bool& isOneForm) const
{
  Tabs tab0;

  MultipleDeriv::const_iterator iter;

  isOneForm = false;  

  if (functionalDeriv.size()==0) return;
  TEST_FOR_EXCEPTION(functionalDeriv.size() > 2, InternalError,
                     "WeakFormBatch::extractWeakForm detected a functional "
                     "derivative of order > 2: " 
                     << functionalDeriv.toString());

  bool foundUnk = false;
  bool foundVar = false;

  SUNDANCE_LEVEL2("extract weak form",
    tab0 << "extracting weak form for functional derivative " 
    << functionalDeriv);

  for (iter = functionalDeriv.begin(); iter != functionalDeriv.end(); iter++)
    {
      Tabs tab;
      const Deriv& d = *iter;
      
      TEST_FOR_EXCEPTION(!d.isFunctionalDeriv(), InternalError,
                         "WeakFormBatch::extractWeakForm "
                         "detected a non-functional derivative: "
                         << functionalDeriv.toString());
      
      const FunctionalDeriv* f = d.funcDeriv();
      
      const SymbolicFuncElement* s 
        = dynamic_cast<const SymbolicFuncElement*>(f->func());

      TEST_FOR_EXCEPTION(s==0, InternalError, 
                         "WeakFormBatch::extractWeakForm failed to cast "
                         "function to SymbolicFuncElement");
      

      int funcID = f->funcComponentID();
      int myIndex = s->myIndex();

      if (!foundUnk && eqn.hasUnkID(funcID))
        {
          const UnknownFuncElement* u
            = dynamic_cast<const UnknownFuncElement*>(s);
          TEST_FOR_EXCEPTION(u==0, InternalError, 
                             "WeakFormBatch::extractWeakForm could not cast "
                             "unknown function to UnknownFuncElement");
          foundUnk = true;
          reducedUnkID = eqn.reducedUnkID(funcID);
          rawUnkID = funcID;
          unkBlock = eqn.blockForUnkID(funcID);

          SUNDANCE_LEVEL2("extract weak form",
            tab << "found reducedUnkID=" << reducedUnkID);

          unkBasis = UnknownFunctionData::getData(u)->basis()[myIndex];
          SUNDANCE_LEVEL2("extract weak form",
            tab << "found unkBasis=" << unkBasis);

          miUnk = f->multiIndex();
          SUNDANCE_LEVEL2("extract weak form",
                       tab << "found unk multi index=" << miUnk.toString());
        }
      else
        {
          foundVar = true;
          reducedVarID = eqn.reducedVarID(funcID);
          rawVarID = funcID;
          testBlock = eqn.blockForVarID(funcID);

          SUNDANCE_LEVEL2("extract weak form",
            tab << "found varID=" << reducedVarID);

          const UnknownFuncElement* u
            = dynamic_cast<const UnknownFuncElement*>(s);

          const TestFuncElement* t
            = dynamic_cast<const TestFuncElement*>(s);

          TEST_FOR_EXCEPTION(u==0 && t==0, InternalError, 
                             "WeakFormBatch::extractWeakForm could not cast "
                             "variational function to either an "
                             "UnknownFuncElement or a TestFuncElement");

          if (t != 0) 
            {
              varBasis = TestFunctionData::getData(t)->basis()[myIndex];
            }
          else
            {
              varBasis = UnknownFunctionData::getData(u)->basis()[myIndex];
            }
          SUNDANCE_LEVEL2("extract weak form", 
            tab << "found varBasis=" << varBasis);

          miVar = f->multiIndex();
          SUNDANCE_LEVEL2("extract weak form", 
            tab << "found var multi index=" << miVar.toString());
        }
    }

  if (!foundUnk) isOneForm = true;
}
