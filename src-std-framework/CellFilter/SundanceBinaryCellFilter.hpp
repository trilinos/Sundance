/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BINARYCELLFILTER_H
#define SUNDANCE_BINARYCELLFILTER_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellFilter.hpp"

namespace SundanceStdFwk
{
  using namespace Teuchos;
  using namespace TSFExtended;
  
  namespace Internal
  {
    /** 
     * BinaryCellFilter implements cell filters that do
     * binary logical operations on cell sets.
     */
    class BinaryCellFilter : public CellFilterBase
    {
    public:

      /** */
      enum CellFilterOpType {Union, Intersection, Difference};

      /** Empty ctor */
      BinaryCellFilter(const CellFilter& left, const CellFilter& right,
                       const CellFilterOpType& op);

      /** virtual dtor */
      virtual ~BinaryCellFilter(){;}

      /** Return the dimension of the cells that will be identified
       * by this filter when acting on the given mesh */
      virtual int dimension(const Mesh& mesh) const ;

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** Ordering */
      virtual bool lessThan(const CellFilterStub* other) const ;

      /** Return the name of the type. Used in ordering. */
      virtual string typeName() const {return "BinaryOpCellFilter";}

      /* */
      GET_RCP(CellFilterStub);
    protected:

      /** Get my cells for the given mesh */
      virtual CellSet internalGetCells(const Mesh& mesh) const ;

    private:
      /** */
      string opName() const ;

      /** The operation I perform */
      CellFilterOpType op_;

      /** My left operand  */
      CellFilter left_;

      /** My right operand */
      CellFilter right_;
    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
