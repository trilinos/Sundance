/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSparsitySuperset.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceTabs.hpp"



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
    numConstantDerivs_(0),
    numVectorDerivs_(0)
{}

void SparsitySuperset::addSubset(const Set<MultiIndex>& multiIndices)
{
  if (hasSubset(multiIndices)) return;
  
  allMultiIndices_.merge(multiIndices);
  RefCountPtr<SparsitySubset> s = rcp(new SparsitySubset(this));
  subsets_.put(multiIndices, s);
}

bool SparsitySuperset::hasSubset(const Set<MultiIndex>& multiIndices) const 
{
  return subsets_.containsKey(multiIndices);
}

const RefCountPtr<SparsitySubset>& SparsitySuperset::subset(const Set<MultiIndex>& multiIndices) const
{
  return subsets_.get(multiIndices);
}

RefCountPtr<SparsitySubset> SparsitySuperset::subset(const Set<MultiIndex>& multiIndices) 
{
  return subsets_.get(multiIndices);
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
  for (Map<Set<MultiIndex>, RefCountPtr<SparsitySubset> >::const_iterator
         i=subsets_.begin(); i != subsets_.end(); i++)
    {
      Tabs tabs1;
      os << tabs1 << "Subset: " << i->first.toString() << endl;
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
