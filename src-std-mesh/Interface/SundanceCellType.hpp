/* @HEADER@ */
/* @HEADER@ */

#ifndef CELLTOPOLOGYCODE_H
#define CELLTOPOLOGYCODE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"

namespace SundanceStdMesh
{
  /** */
  enum CellType {NullCell, PointCell, LineCell, TriangleCell, 
                         TetCell, QuadCell, BrickCell, PrismCell};

  /** \relates CellType */
  string toString(const CellType& c) ;

  /** \relates CellType */
  int dimension(const CellType& c) ;


  /** \relates CellType */
  inline ostream& operator<<(ostream& os, const CellType& c)
  {
    os << toString(c);
    return os;
  }

  
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */   

#endif


