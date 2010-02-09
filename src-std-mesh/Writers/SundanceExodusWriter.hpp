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

#ifndef SUNDANCE_EXODUSWRITER_H
#define SUNDANCE_EXODUSWRITER_H


#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"
#include "SundanceCellType.hpp"

namespace SundanceStdMesh
{
  /**
   * ExodusWriter writes a mesh or fields to an ExodusII file
   */
  class ExodusWriter : public FieldWriterBase
  {
  public:
    /** */
    ExodusWriter(const string& filename) 
      : FieldWriterBase(filename) {;}
    
    /** virtual dtor */
    virtual ~ExodusWriter(){;}

    /** */
    virtual void write() const ;

    /** Return a ref count pointer to self */
    virtual RefCountPtr<FieldWriterBase> getRcp() {return rcp(this);}

    /** */
    void writeParallelInfo(const string& filename) const ;


  private:    
    /** */
    void getCharpp(const Array<std::string>& s, Array<const char*>& p) const ;

    /** */
    void findNodeSets(
      Array<CellFilter>& nodesetFilters,
      Array<int>& omnipresentFuncs,
      Array<RefCountPtr<Array<int> > >& funcsForNodeset,
      Array<RefCountPtr<Array<int> > >& nodesForNodeset,
      Array<int>& nsID,
      Array<int>& nNodesPerSet,
      Array<int>& nsNodePtr,
      RefCountPtr<Array<int> > allNodes
      ) const ;

    /** */
    void findBlocks(
      Array<CellFilter>& blockFilters,
      Array<int>& omnipresentFuncs,
      Array<RefCountPtr<Array<int> > >& funcsForBlock,
      Array<RefCountPtr<Array<int> > >& elemsForBlock,
      Array<int>& elemIDs,
      Array<int>& nElemsPerBlock,
      Array<int>& blockElemPtr,
      RefCountPtr<Array<int> > allElems
      ) const ;



    /** */
    void offset(Array<int>& x) const ;

    /** */
    std::string elemType(const CellType& type) const ;

    /** */
    void writeMesh(int exoID, 
      const Array<CellFilter>& nodesetFilters,
      const Array<int>& nsID,
      const Array<int>& nNodesPerSet,
      const Array<int>& nsNodePtr,
      const RefCountPtr<Array<int> >& allNodes) const ;

    /** */
    void writeFields(int exoID, 
      const Array<CellFilter>& nodesetFilters,
      const Array<int>& omnipresentNodalFuncs,
      const Array<int>& omnipresentElemFuncs,
      const Array<RefCountPtr<Array<int> > >& funcsForNodeset,
      const Array<RefCountPtr<Array<int> > >& nodesForNodeset,
      const Array<int>& nsID) const ;
    
    
  };
}


#endif
