/* @HEADER@ */
/* @HEADER@ */


#ifndef SUNDANCE_CELLPREDICATEBASE_H
#define SUNDANCE_CELLPREDICATEBASE_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "TSFHandleable.hpp"
#include "Teuchos_XMLObject.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include <typeinfo>

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;
  
  namespace Internal
  {
    /** 
     * CellPredicateBase is the base class for predicate objects
     * that test cells
     * against some condition. A simulation developer needing
     * some specialized method for identifying cells might implement
     * a custom cell predicate by extending this function. However,
     * the most common cases, selection by cell label or cell position,
     * have already been implemented 
     * in LabelCellPredicate and PositionalCellPredicate.
     */
    class CellPredicateBase 
      : public TSFExtended::Handleable<CellPredicateBase>,
        public Noncopyable,
        public TSFExtended::ObjectWithVerbosity<CellPredicateBase>
    {
    public:
      /** Empty ctor */
      CellPredicateBase();

      /** virtual dtor */
      virtual ~CellPredicateBase(){;}
      
      /** Test whether the cell with the given LID satisfies the condition */
      virtual bool test(int cellLID) const = 0 ;

      
      /** Set the current mesh and dimension on which cells are to be tested */
      virtual void setMesh(const Mesh& mesh, int cellDim) const 
      {mesh_ = mesh; cellDim_ = cellDim;}

      /** Write to XML */
      virtual XMLObject toXML() const = 0 ;

      /** */
      virtual bool lessThan(const CellPredicateBase* other) const = 0 ;


      /** */
      virtual string typeName() const {return typeid(*this).name();}
    protected:

      /** */
      const Mesh& mesh() const {return mesh_;}


      /** */
      int cellDim() const {return cellDim_;}

    private:
      mutable Mesh mesh_;

      mutable int cellDim_;



    };
  }
}

#endif
