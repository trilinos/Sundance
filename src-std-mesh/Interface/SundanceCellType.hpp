/* @HEADER@ */
/* @HEADER@ */

#ifndef CELLTOPOLOGYCODE_H
#define CELLTOPOLOGYCODE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"

namespace Sundance
{
  /** */
  enum CellTopologyCode {NullCell, PointCell, LineCell, TriangleCell, 
                         TetCell, QuadCell, BrickCell, PrismCell};

  /** \relates CellTopologyCode */
  string toString(const CellTopologyCode& c) ;

  /** \relates CellTopologyCode */
  int dimension(const CellTopologyCode& c) ;


  /** \relates CellTopologyCode */
  inline ostream& operator<<(ostream& os, const CellTopologyCode& c)
  {
    os << toString(c);
    return os;
  }

  
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */   

#endif


