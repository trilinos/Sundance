/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLREORDERERIMPLEMBASE_H
#define SUNDANCE_CELLREORDERERIMPLEMBASE_H


#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceNoncopyable.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include <typeinfo>

namespace SundanceStdMesh
{
  namespace Internal
  {
    class MeshBase;

    /**
     * Abstract interface for the low-level objects that 
     * implement cell reordering. 
     * 
     * <h4> Adding a new reordering algorithm </h4>
     *
     * To add a new reordering algorithm, you should create a new
     * subclass of CellReordererImplemBase. The only method you will
     * need to implement is
     * \code
     * virtual int advance(int currentLID) const 
     * \endcode
     * which should provide the maximal cell LID found after
     * the <tt>currentLID.</tt>
     * Depending on the algorithm , you may also want to override
     * the methods
     * \code
     * virtual int begin() const 
     * virtual int end() const 
     * \endcode
     * which return the index of the first cell to be processed,
     * and a past-the-end index. 
     */
    class CellReordererImplemBase 
      : public TSFExtended::ObjectWithVerbosity<CellReordererImplemBase>
    {
    public:
      /** Construct with a pointer to a mesh */
      CellReordererImplemBase(const MeshBase* mesh);
      
      /** virtual dtor */
      virtual ~CellReordererImplemBase(){;}

      /** return a descriptive string */
      virtual string typeName() const {return typeid(*this).name();}
    
      /** */
      virtual int advance(int currentLID) const = 0 ;
      
      /** */
      virtual int begin() const {return 0;}
      
      /** */
      virtual int end() const ;
    protected:
      /** */
      const MeshBase* mesh() const {return mesh_;}

    private:
      /** Unmanaged pointer to a mesh. The mesh will contain a smart
       * pointer to this reorderer, so to avoid closed reference
       * graphs we store a raw pointer here.*/
      const MeshBase* mesh_;
      
    };

  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
