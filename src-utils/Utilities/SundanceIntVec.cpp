#include "SundanceIntVec.hpp"
#include "SundanceDebug.hpp"
#include "SundanceCombinatorialUtils.hpp"

namespace Sundance
{
using Teuchos::Array;

IntVec::IntVec(int n)
  : data_(n)
{
  for (int i=0; i<n; i++) data_[i] = 0;
}

IntVec IntVec::operator+(const IntVec& other) const
{
  TEST_FOR_EXCEPT(size() != other.size());

  IntVec rtn(size());
  for (int i=0; i<size(); i++) rtn[i] = data_[i] + other[i];
  
  return rtn;
}


IntVec IntVec::operator*(int a) const
{
  IntVec rtn(size());
  for (int i=0; i<size(); i++) rtn[i] = a*data_[i];
  
  return rtn;
}

int IntVec::factorial() const
{
  int rtn=1;

  for (int i=0; i<size(); i++)
  {
    int n_i = data_[i];
    for (int j=1; j<=n_i; j++) rtn *= j;
  }
  return rtn;
}


int IntVec::pow(const IntVec& other) const
{
  TEST_FOR_EXCEPT(size() != other.size());
  int rtn=1;

  for (int i=0; i<size(); i++)
  {
    int n_i = data_[i];
    int p_i = other[i];
    for (int j=1; j<=p_i; j++) rtn *= n_i;
  }
  return rtn;
}

int IntVec::abs() const
{
  int rtn = 0;
  for (int i=0; i<size(); i++)
  {
    rtn += ::abs(data_[i]);
  }
  return rtn;
}

int IntVec::norm() const
{
  int rtn = 0;
  for (int i=0; i<size(); i++)
  {
    if (::abs(data_[i]) > rtn) rtn = ::abs(data_[i]);
  }
  return rtn;
}

bool IntVec::operator==(const IntVec& other) const
{
  if (size() != other.size()) return false;
  for (int i=0; i<size(); i++)
  {
    if (data_[i] != other.data_[i]) return false;
  }
  return true;
}

bool IntVec::operator<(const IntVec& other) const
{
  if (size() < other.size()) return true;
  if (size() > other.size()) return false;
  if (abs() < other.abs()) return true;
  if (abs() > other.abs()) return false;

  for (int i=0; i<size(); i++)
  {
    if (data_[i] < other.data_[i]) return true;
    if (data_[i] > other.data_[i]) return false;
  }

  return false;
}

void IntVec::print(std::ostream& os) const 
{
  os << "IVec[" << data_ << "]";
}


void IntVec::getPartitions(int M, Array<Array<IntVec> >& parts) const
{
  Array<Array<Array<int> > > rComp(size());
  Array<int> radix(size());
  
  for (int i=0; i<size(); i++)
  {
    restrictedCompositions(data_[i], M, rComp[i]);
    radix[i] = rComp[i].size();
  }

  Array<int> pick(size(), 0);
  bool workLeft = true;
  while (workLeft)
  {
    Array<IntVec> part(M);
    for (int j=0; j<M; j++) part[j] = IntVec(size());
    for (int i=0; i<size(); i++)
    {
      TEST_FOR_EXCEPT(pick[i] >= rComp[i].size());
      const Array<int>& p = rComp[i][pick[i]];
      TEST_FOR_EXCEPT(p.size() != M);
      for (int j=0; j<M; j++) part[j][i] = p[j];
    }
    bool isNonzero = true;
    for (int i=0; i<part.size(); i++) 
    {
      if (part[i].abs()==0) {isNonzero = false; break;}
    }
    if (isNonzero) parts.append(part);
    workLeft = nextNum(pick, radix);
  }
}

}



	





