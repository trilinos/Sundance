/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_VERBOSEFIELDWRITER_H
#define SUNDANCE_VERBOSEFIELDWRITER_H


#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"

namespace SundanceStdMesh
{
  /**
   * VerboseFieldWriter dumps every imaginable bit of information
   * about a mesh. Used for debugging mesh data structures. 
   */
  class VerboseFieldWriter : public FieldWriterBase
  {
  public:
    /** */
    VerboseFieldWriter(const string& filename="") 
      : FieldWriterBase(filename) {;}
    
    /** virtual dtor */
    virtual ~VerboseFieldWriter(){;}

    /** */
    virtual void write() const ;

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<FieldWriterBase> getRcp() {return rcp(this);}
#endif /* DOXYGEN_DEVELOPER_ONLY */
  };
}




#endif
