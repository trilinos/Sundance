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

#ifndef SUNDANCE_CELLSET_H
#define SUNDANCE_CELLSET_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellSetBase.hpp"
#include "SundanceCellPredicate.hpp"
#include "TSFHandle.hpp"


namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;
  using namespace TSFExtended;

  namespace Internal
  {

    /** 
     * CellSet is, you guessed it, a set of cells in a mesh. Cells are 
     * represented by their LID relative to the mesh. 
     * 
     * 
     * @see CellFilter, CellIterator
     **/
    class CellSet : public TSFExtended::Handle<CellSetBase>
    {
    public:
      /* handle boilerplate */
      HANDLE_CTORS(CellSet, CellSetBase);

      /** The ID number of the mesh in which these cells exist */
      int meshID() const {return ptr()->meshID();}
      
      /** The mesh in which these cells exist */
      const Mesh& mesh() const {return ptr()->mesh();}

      /** The dimension of the cells contained in this set */
      int dimension() const {return ptr()->dimension();}

      /** The type of the cells contained in this set */
      const CellType& cellType() const {return ptr()->cellType();}

      /** An iterator pointing to the beginning of the set */
      CellIterator begin() const {return ptr()->begin();}

      /** An iterator pointing to the end of the set */
      CellIterator end() const {return ptr()->end();}

      /** Return a cell set that is the union of this set and another set */
      CellSet setUnion(const CellSet& other) const ;

      /** Return a cell set that is the intersection
       *  of this set and another set */
      CellSet setIntersection(const CellSet& other) const ;

      /** Return a cell set that is the difference
       *  of this set and another set */
      CellSet setDifference(const CellSet& other) const ;

      /** */
      CellSet subset(const RefCountPtr<CellPredicate>& test) const ;
      
    private:
      void checkCompatibility(const string& op, const CellSet& other) const ;
    };


  }

  STREAM_OUT(Internal::CellSet);
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
