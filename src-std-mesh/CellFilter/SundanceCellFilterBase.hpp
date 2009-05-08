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

#ifndef SUNDANCE_CELLFILTERBASE_H
#define SUNDANCE_CELLFILTERBASE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceMap.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore;
  using namespace Teuchos;
  
  namespace Internal
  {
    /** 
     * Base class for CellFilter objects.
     *
     * <h4> Notes for subclass implementors </h4>
     * 
     * Derived classes must implement the methods
     * <ul>
     * <li> internalGetCells() -- returns the set of cells that 
     * pass through this filter
     * <li> dimension() -- returns the dimension of the cells that
     * will pass through this filter
     * <li> toXML() -- writes an XML description of the filter
     * <li> lessThan() -- compares to another cell filter. Used to store
     * cell filters in STL containers. 
     * <li> typeName() -- returns the name of the subclass. Used in ordering.
     * </ul>
     */
    class CellFilterBase : public CellFilterStub
    {
    public:
      /** Empty ctor */
      CellFilterBase();

      /** virtual dtor */
      virtual ~CellFilterBase(){;}

      /** Find the cells passing this filter on the given mesh. This
       * method will cache the cell sets it computes for each mesh  */
      CellSet getCells(const Mesh& mesh) const ;

      /** Return the dimension of the cells that will be identified
       * by this filter when acting on the given mesh */
      virtual int dimension(const Mesh& mesh) const = 0 ;

      /** */
      void registerSubset(const CellFilter& sub) const ;

      /** */
      void registerLabeledSubset(int label, const CellFilter& sub) const ;

      /** */
      void registerDisjoint(const CellFilter& sub) const ;

      /** */
      const Set<CellFilter>& knownSubsets() const {return subsets_;}

      /** */
      const Set<CellFilter>& knownDisjoints() const {return disjoints_;}

      /** */
      virtual string toString() const {return name_;}

      /** Print to a stream */
      virtual string description() const 
      {return toString();}

      /** */
      void setName(const string& name) {name_ = name;}

    protected:

      /** */
      virtual CellSet internalGetCells(const Mesh& mesh) const = 0 ;

    private:
      /** cache of previously computed cell sets */
      mutable CellSet cellSetCache_;

      /** */
      mutable Set<CellFilter> subsets_;

      /** */
      mutable Set<CellFilter> disjoints_;

      /** */
      mutable SundanceUtils::Map<int, CellFilter> labeledSubsets_;

      /** */
      string name_;

    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
