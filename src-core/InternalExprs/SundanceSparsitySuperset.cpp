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
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"



using namespace SundanceCore;
using namespace SundanceUtils;
using namespace TSFExtended;
using namespace SundanceCore::Internal;
using namespace Teuchos;




SparsitySuperset::SparsitySuperset()
  : maxOrder_(0),
    derivToIndexMap_(),
    derivs_(),
    states_(),
    multiIndex_(),
    subsets_(),
    allMultiIndices_(),
    allFuncIDs_(),
    numConstantDerivs_(0),
    numVectorDerivs_(0)
{}

void SparsitySuperset::addSubset(const Set<MultiIndex>& multiIndices,
                                 const Set<MultiSet<int> >& funcIDs)
{
  if (hasSubset(multiIndices, funcIDs)) return;
  
  allMultiIndices_.merge(multiIndices);
  allFuncIDs_.merge(funcIDs);
  RefCountPtr<SparsitySubset> s = rcp(new SparsitySubset(this));
  subsets_.put(keyPair(multiIndices, funcIDs), s);
}

bool SparsitySuperset::hasSubset(const Set<MultiIndex>& multiIndices,
                                 const Set<MultiSet<int> >& funcIDs) const 
{
  return subsets_.containsKey(keyPair(multiIndices, funcIDs));
}

const RefCountPtr<SparsitySubset>& 
SparsitySuperset::subset(const Set<MultiIndex>& multiIndices,
                         const Set<MultiSet<int> >& funcIDs) const
{
  return subsets_.get(keyPair(multiIndices, funcIDs));
}

RefCountPtr<SparsitySubset> 
SparsitySuperset::subset(const Set<MultiIndex>& multiIndices,
                         const Set<MultiSet<int> >& funcIDs) 
{
  return subsets_.get(keyPair(multiIndices, funcIDs));
}

RefCountPtr<SparsitySubset> 
SparsitySuperset::findSubset(const Set<MultiIndex>& multiIndices,
                             const Set<MultiSet<int> >& funcIDs) const
{
  /* first see if the subset exists */
  if (hasSubset(multiIndices, funcIDs))
    {
      return subsets_.get(keyPair(multiIndices, funcIDs));
    }

  /* next, try to find it as the union of two or more existing subsets */
  assembleSubsetUnions();
  
  TEST_FOR_EXCEPTION(!hasSubset(multiIndices, funcIDs),
                     InternalError,
                     "sparsity subset for MI=" << multiIndices.toString()
                     << ", funcID=" << funcIDs.toString()
                     << " not found in superset " 
                     << *this);

  return subsets_.get(keyPair(multiIndices, funcIDs));
}

void SparsitySuperset::assembleSubsetUnions() const
{
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "assembling subset unions");
  typedef Map<keyPair, RefCountPtr<SparsitySubset> >::const_iterator iter;

  int oldSize = subsets_.size();

  for (iter i=subsets_.begin(); i != subsets_.end(); i++)
    {
      const keyPair& key1 = i->first;
      const RefCountPtr<SparsitySubset>& sub1 = i->second;
      for (iter j=subsets_.begin(); j != i; j++)
        {
          Tabs tab1;
          const keyPair& key2 = j->first;
          const RefCountPtr<SparsitySubset>& sub2 = j->second;
          Set<MultiIndex> newMiSet = key1.first();
          newMiSet.merge(key2.first());
          Set<MultiSet<int> > newFuncIDSet = key1.second();
          newFuncIDSet.merge(key2.second());
          keyPair newKey(newMiSet, newFuncIDSet);
          /* if the union is not a new subset, try again */
          if (subsets_.containsKey(newKey)) continue;
          /* add the new subset */
          SparsitySuperset* me = const_cast<SparsitySuperset*>(this);
          RefCountPtr<SparsitySubset> s = rcp(new SparsitySubset(me));
          subsets_.put(newKey, s);
          for (int k=0; k<sub1->numDerivs(); k++)
            {
              s->addDeriv(sub1->deriv(k), sub1->state(k));
            }
          
          for (int k=0; k<sub2->numDerivs(); k++)
            {
              s->addDeriv(sub2->deriv(k), sub2->state(k));
            }
          SUNDANCE_VERB_HIGH(tab1 << "added new subset " 
                             << endl << tab1 << *s
                             << endl << tab1 << "with key " 
                             "miSet=" << newMiSet.toString()
                             << ", funcIDSet=" << newFuncIDSet.toString());
        }
    }

  /* if we have created new subsets, continue the process by recursion
   * until all subset combinations are enumerated */
  if (oldSize != subsets_.size())
    {
      SUNDANCE_VERB_HIGH(tab << "recursing");
      assembleSubsetUnions();
    }
}



void SparsitySuperset::addDeriv(const MultipleDeriv& d, 
                               const DerivState& state)
{
  maxOrder_ = max(d.order(), maxOrder_);

  if (containsDeriv(d))
    {
      const DerivState& oldState = states_[getIndex(d)];
      if (state > oldState) 
        {
          states_[getIndex(d)] = state;
          numConstantDerivs_--;
          numVectorDerivs_++;
        }
    }
  else
    {
      int index = derivs_.size();
      derivs_.append(d);
      states_.append(state);
      derivToIndexMap_.put(d, index);
      MultiIndex mi;
      for (MultipleDeriv::const_iterator i=d.begin(); 
           i != d.end(); i++)
        {
          if (i->isCoordDeriv())
            {
              MultiIndex m;
              int dir = i->coordDeriv()->dir();
              m[dir] = 1;
              mi = mi + m;
            }

        }
      multiIndex_.append(mi);
      if (state==ConstantDeriv) numConstantDerivs_++;
      else numVectorDerivs_++;
    }
}

void SparsitySuperset::addDeriv(const Deriv& d, 
                               const DerivState& state)
{
  MultipleDeriv md;
  md.put(d);
  addDeriv(md, state);
}



bool SparsitySuperset::containsDeriv(const MultipleDeriv& d) const
{
  return derivToIndexMap_.containsKey(d);
}

int SparsitySuperset::getIndex(const MultipleDeriv& d) const
{
  if (!containsDeriv(d)) return -1;
  return derivToIndexMap_.get(d);
}



void SparsitySuperset::print(ostream& os,
                             const Array<RefCountPtr<EvalVector> >& vecResults,
                             const Array<double>& constantResults) const 
{
  Tabs tabs;

  /* find the maximum size of the string reps of the derivatives.
   * We'll use this to set the field width for printing derivatives. */
  int maxlen = 25;
  for (unsigned int i=0; i<derivs_.size(); i++)
    {
      int s = derivs_[i].toString().length();
      if (s > maxlen) maxlen = s;
    }


  int vecIndex=0;
  int constIndex = 0;
  os << tabs << "Results Superset" << endl;
  for (unsigned int i=0; i<derivs_.size(); i++)
    {
      os << tabs << i << "\t\t";
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << derivs_[i].toString() << "\t\t" ;
      switch(states_[i])
        {
        case ZeroDeriv:
          os  << "Zero" << endl;
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



void SparsitySuperset::displayAll(ostream& os) const 
{
  Tabs tabs;

  /* find the maximum size of the string reps of the derivatives.
   * We'll use this to set the field width for printing derivatives. */
  int maxlen = 25;
  for (unsigned int i=0; i<derivs_.size(); i++)
    {
      int s = derivs_[i].toString().length();
      if (s > maxlen) maxlen = s;
    }


  os << tabs << "SparsitySuperset" << endl;
  for (unsigned int i=0; i<derivs_.size(); i++)
    {
      os << tabs << i << "\tderiv=\t";
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << derivs_[i].toString() << "\tstate=\t" ;
      switch(states_[i])
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

  os << tabs << "Subsets" << endl;
  for (Map<keyPair, RefCountPtr<SparsitySubset> >::const_iterator
         i=subsets_.begin(); i != subsets_.end(); i++)
    {
      Tabs tabs1;
      os << tabs1 << "Subset: " << i->first << endl;
      i->second->print(os);
    }
}

void SparsitySuperset::print(ostream& os) const 
{
  Tabs tabs;

  /* find the maximum size of the string reps of the derivatives.
   * We'll use this to set the field width for printing derivatives. */
  int maxlen = 25;
  for (unsigned int i=0; i<derivs_.size(); i++)
    {
      int s = derivs_[i].toString().length();
      if (s > maxlen) maxlen = s;
    }


  os << tabs << "SparsitySuperset" << endl;
  for (unsigned int i=0; i<derivs_.size(); i++)
    {
      os << tabs << i << "\tderiv=\t";
      os.width(maxlen);
      os.setf(ios_base::left, ios_base::adjustfield);
      os << derivs_[i].toString() << "\tstate=\t" ;
      switch(states_[i])
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

  os << tabs << "Subsets" << endl;
  for (Map<keyPair, RefCountPtr<SparsitySubset> >::const_iterator
         i=subsets_.begin(); i != subsets_.end(); i++)
    {
      Tabs tabs1;
      os << tabs1 << "Subset: " << i->first << endl;
      i->second->print(os);
    }

}

string SparsitySuperset::toString() const 
{
	ostringstream ss;
	print(ss);
	return ss.str();
}

DerivSet SparsitySuperset::derivSet() const
{
  DerivSet rtn;
  for (int i=0; i<numDerivs(); i++) 
    {
      if (state(i) != ZeroDeriv) rtn.put(deriv(i));
    }
  return rtn;
}
