/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLFILTERBASE_H
#define SUNDANCE_CELLFILTERBASE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellSet.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
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
    class CellFilterBase : public TSFExtended::Handleable<CellFilterBase>,
                           public TSFExtended::Printable,
                           public TSFExtended::Describable,
                           public Noncopyable,
                           public TSFExtended::ObjectWithVerbosity<CellFilterBase>
                           
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

      /** Write to XML */
      virtual XMLObject toXML() const = 0 ;

      /** Defines how this filter is */
      virtual bool lessThan(const CellFilterBase* other) const = 0 ;

      /** Return the name of the type. Used in ordering. */
      virtual string typeName() const = 0 ;
      
      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(ostream& os) const {os << toXML();}
      //@}

      /** \name Describable interface */
      //@{
      /** Print to a stream */
      virtual string describe() const 
      {return typeName();}
      //@}
    
    protected:

      /** */
      virtual CellSet internalGetCells(const Mesh& mesh) const = 0 ;

    private:
      /** cache of previously computed cell sets */
      mutable Map<int, CellSet> cellSetCache_;

    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
