/* @HEADER@ */
/* @HEADER@ */

#include "SundanceTabs.hpp"

using namespace SundanceUtils;


Tabs::Tabs()
{
  tabLevel()++;
}

Tabs::~Tabs()
{
  tabLevel()--;
}

void Tabs::print(ostream& os) const
{
  os << "[" << tabLevel() << "]";
  for (int i=0; i<tabLevel(); i++) os << "==";
  os << "> ";
}



