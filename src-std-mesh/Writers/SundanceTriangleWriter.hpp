/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TRIANGLEWRITER_H
#define SUNDANCE_TRIANGLEWRITER_H


#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"

namespace SundanceStdMesh
{
  /**
   * TriangleWriter writes a mesh or fields to a file
   * in Shewchuk's Triangle format.
   */
  class TriangleWriter : public FieldWriterBase
  {
  public:
    /** */
    TriangleWriter(const string& filename="", int indexOffset=0) 
      : FieldWriterBase(filename), indexOffset_(indexOffset) {;}
    
    /** virtual dtor */
    virtual ~TriangleWriter(){;}

    /** */
    virtual void write() const ;

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<FieldWriterBase> getRcp() {return rcp(this);}
#endif /* DOXYGEN_DEVELOPER_ONLY */

  protected:
    /** */
      void writePoints(const string& filename) const ;

    /** */
    void writeCells(const string& filename) const ;
    
    /** */
    void writeEdges(const string& filename) const ;
    
    /** */
    void writeFaces(const string& filename) const ;
    
    /** */
    void writeHeader(const string& filename) const ;
    
    /** */
    void writeParallelInfo(const string& filename) const ;
    
    int indexOffset_;
  };
}




#endif
