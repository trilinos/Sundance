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

#include "SundanceSparsitySuperset.hpp"
#include "SundanceSparsitySubset.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceTabs.hpp"



using namespace SundanceCore;
using namespace SundanceUtils;
using namespace TSFExtended;
using namespace SundanceCore::Internal;
using namespace Teuchos;


SparsitySubset::SparsitySubset(SparsitySuperset* master)
  : master_(master),
    derivs_(),
    derivToIndexMap_(),
    ptr_()
{}

SparsitySuperset* SparsitySubset::master() 
{
  TEST_FOR_EXCEPTION(master_==0, RuntimeError,
                     "null pointer to master sparsity pattern detected "
                     "in SparsitySubset::master()");
  return master_;
}

const SparsitySuperset* SparsitySubset::master() const 
{
  TEST_FOR_EXCEPTION(master_==0, RuntimeError,
                     "null pointer to master sparsity pattern detected "
                     "in SparsitySubset::master()");
  return master_;
}

void SparsitySubset::addDeriv(const MultipleDeriv& d, 
                              const DerivState& state)
{
  master()->addDeriv(d, state);
  if (!containsDeriv(d))
    {
      int index = derivs_.size();
      derivs_.append(d);
      ptr_.append(master()->getIndex(d));
      derivToIndexMap_.put(d, index);
    }
}

void SparsitySubset::addDeriv(const Deriv& d, 
                              const DerivState& state)
{
  MultipleDeriv md;
  md.put(d);
  addDeriv(md, state);
}



bool SparsitySubset::containsDeriv(const MultipleDeriv& d) const
{
  return derivToIndexMap_.containsKey(d);
}

int SparsitySubset::getIndex(const MultipleDeriv& d) const
{
  if (containsDeriv(d))
    {
      return derivToIndexMap_.get(d);
    }
  return -1;
}

int SparsitySubset::numDerivs() const 
{
  return derivs_.size();
}

bool SparsitySubset::isConstant(int i) const
{
  return master()->isConstant(ptr_[i]);
}

bool SparsitySubset::isSpatialDeriv(int i) const
{
  return master()->isSpatialDeriv(ptr_[i]);
}

const MultiIndex& SparsitySubset::multiIndex(int i) const 
{
  return master()->multiIndex(ptr_[i]);
}

const MultipleDeriv& SparsitySubset::deriv(int i) const 
{
  return master()->deriv(ptr_[i]);
}


const DerivState& SparsitySubset::state(int i) const 
{
  return master()->state(ptr_[i]);
}




void SparsitySubset::print(ostream& os) const 
{
  Tabs tabs;

  /* find the maximum size of the string reps of the derivatives.
   * We'll use this to set the field width for printing derivatives. */
  int maxlen = 25;
  for (int i=0; i<numDerivs(); i++)
    {
      int s = deriv(i).toString().length();
      if (s > maxlen) maxlen = s;
    }


  os << tabs << "SparsitySubset" << endl;
  for (int i=0; i<numDerivs(); i++)
    {
      os << tabs << i << "\tderiv=\t";
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << deriv(i).toString() << "\tstate=\t" ;
      switch(state(i))
        {
        case ZeroDeriv:
          os  << "Zero" << endl;
          break;
        case ConstantDeriv:
          os << "Constant" << endl;
          break;
        case VectorDeriv:
          os << "Vector" << endl;
          break;
        }
    }
}

void SparsitySubset::print(ostream& os, 
                           const Array<RefCountPtr<EvalVector> >& vecResults,
                           const Array<double>& constantResults) const 
{
  Tabs tabs;

  /* find the maximum size of the string reps of the derivatives.
   * We'll use this to set the field width for printing derivatives. */
  int maxlen = 25;
  for (int i=0; i<numDerivs(); i++)
    {
      int s = deriv(i).toString().length();
      if (s > maxlen) maxlen = s;
    }


  int vecIndex=0;
  int constIndex = 0;
  os << tabs << "Results Subset" << endl;
  for (int i=0; i<numDerivs(); i++)
    {
      os << tabs << i << "\t\t";
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << deriv(i).toString() << "\t\t";
      switch(state(i))
        {
        case ZeroDeriv:
          break;
        case ConstantDeriv:
          os << constantResults[constIndex++] << endl;
          break;
        case VectorDeriv:
          if (vecResults[vecIndex].get()==0)
            {
              os << "{Null}";
            }
          else
            {
              vecResults[vecIndex]->print(os);
            }
          vecIndex++;
          os << endl;
          break;
        }
    }
}




string SparsitySubset::toString() const 
{
	ostringstream ss;
	print(ss);
	return ss.str();
}

DerivSet SparsitySubset::derivSet() const
{
  DerivSet rtn;
  for (int i=0; i<numDerivs(); i++) 
    {
      if (state(i) != ZeroDeriv) rtn.put(deriv(i));
    }
  return rtn;
}
