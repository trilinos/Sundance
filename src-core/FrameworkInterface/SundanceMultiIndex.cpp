/* @HEADER@ */
/* @HEADER@ */

#include "SundanceMultiIndex.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

MultiIndex::MultiIndex()
	: m_(maxDim(), 0)
{;}

MultiIndex::MultiIndex(int x, int y, int z)
	: m_(maxDim(), 0)
{
	m_[0] = x;
	m_[1] = y;
	m_[2] = z;
}

MultiIndex MultiIndex::operator+(const MultiIndex& other) const 
{
	MultiIndex rtn;

	for (int i=0; i<maxDim(); i++)
		{
			rtn.m_[i] = m_[i] + other[i];
		}
	return rtn;
}

string MultiIndex::toString() const
{
	return "(" + Teuchos::toString(m_[0]) + ","
		+ Teuchos::toString(m_[1]) + ","
		+ Teuchos::toString(m_[2]) + ")";
} 

XMLObject MultiIndex::toXML() const
{
	XMLObject rtn("MultiIndex");
	rtn.addAttribute("indices", toString());
	return rtn;
}

bool MultiIndex::operator==(const MultiIndex& m) const 
{
	for (int i=0; i<maxDim(); i++)
		{
			if (m_[i] != m[i]) return false;
		}
	return true;
}

bool MultiIndex::operator<(const MultiIndex& m) const 
{
	for (int i=0; i<maxDim(); i++)
		{
			if (m_[i] > m.m_[i]) return false;
      if (m_[i] < m.m_[i]) return true;
		}
	return false;
}



int MultiIndex::order() const 
{
  int h = 0;
	for (int i=0; i<maxDim(); i++)
		{
			h += m_[i];
		}
	return h;
}

int MultiIndex::firstOrderDirection() const 
{
  TEST_FOR_EXCEPTION(order() != 1, InternalError,
                     "bad order in MultiIndex::firstOrderDirection() const");
	for (int i=0; i<maxDim(); i++)
		{
			if (m_[i] == 1) return i;
		}
  return -1;
}



string MultiIndex::coordForm() const
{
  string rtn;
  
  for (int i=0; i<m_[0]; i++)
		{
      rtn += "x";
		}
  for (int i=0; i<m_[1]; i++)
		{
      rtn += "y";
		}
  for (int i=0; i<m_[2]; i++)
		{
      rtn += "z";
		}
	return rtn;
}


