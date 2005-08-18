#include "SundanceIntHashSet.hpp"
#include "SundanceExceptions.hpp"


using namespace SundanceUtils;
using namespace Teuchos;

IntHashSet::IntHashSet()
  : capacity_(0),
    data_(),
    size_(0)
{;}

void IntHashSet::setCapacity(unsigned int capacity)
{
  capacity_ = capacity;
  data_.resize(capacity_);
}
    



bool IntHashSet::contains(int x) const 
{
  int ptr = hashFunc(x);
  const std::list<int>& d = data_[ptr];
  for (std::list<int>::const_iterator i=d.begin(); i != d.end(); i++)
    {
      if (x == *i) return true;
    }
  return false;
}

void IntHashSet::fillArray(int* a) const
{
  int k = 0;
  for (unsigned int i=0; i<data_.size(); i++)
    {
      for (std::list<int>::const_iterator 
             j=data_[i].begin(); j!=data_[i].end(); j++, k++)
        {
          a[k] = *j;
        }
    }
}



