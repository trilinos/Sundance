/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaTabs.hpp"

using namespace Playa;


Tabs::Tabs(bool jump)
  : jump_(jump), myLevel_(0)
{
  if (jump) tabLevel()++;
  myLevel_ = tabLevel();
}

Tabs::~Tabs()
{
  if (jump_) tabLevel()--;
}

void Tabs::print(std::ostream& os) const
{
  if (showDepth()) 
    {
      os << "[" << tabLevel() << "]";
    }
  int n = myLevel_ * tabSize();
  for (int i=0; i<n; i++) os << " ";
}



