/* @HEADER@ */
/* @HEADER@ */

#include "SundanceTabs.hpp"

using namespace SundanceUtils;


Tabs::Tabs(char c)
  : c_(c)
{
  tabLevel()++;
}

Tabs::~Tabs()
{
  tabLevel()--;
}

void Tabs::print(ostream& os) const
{
  if (showDepth()) 
    {
      os << "[" << tabLevel() << "]";
    }
  int n = tabLevel() * tabSize();
  for (int i=0; i<n; i++) os << c_;
}



