/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MATLABWRITER_H
#define SUNDANCE_MATLABWRITER_H


#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"

namespace SundanceStdMesh
{
  /**
   * MatlabWriter writes a 1D mesh to a matlab file
   */
  class MatlabWriter : public FieldWriterBase
  {
  public:
    /** */
    MatlabWriter(const string& filename="") 
      : FieldWriterBase(filename) {;}
    
    /** virtual dtor */
    virtual ~MatlabWriter(){;}

    /** */
    virtual void write() const ;

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<FieldWriterBase> getRcp() {return rcp(this);}


  private:
#endif /* DOXYGEN_DEVELOPER_ONLY */
  };
}




#endif
