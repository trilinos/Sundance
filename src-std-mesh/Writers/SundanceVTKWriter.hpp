/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_VTKWRITER_H
#define SUNDANCE_VTKWRITER_H


#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"

namespace SundanceStdMesh
{
  /**
   * VTKWriter writes a mesh or fields to a VTK file
   */
  class VTKWriter : public FieldWriterBase
  {
  public:
    /** */
    VTKWriter(const string& filename="") 
      : FieldWriterBase(filename) {;}
    
    /** virtual dtor */
    virtual ~VTKWriter(){;}

    /** */
    virtual void write() const ;

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<FieldWriterBase> getRcp() {return rcp(this);}


  private:
    /** */
    void lowLevelWrite(const string& filename, bool isPHeader) const ;

    /** */
    void writePoints(ostream& os, bool isPHeader) const ;

    /** */
    void writeCells(ostream& os) const ;

    /** */
    void writePointData(ostream& os, bool isPHeader) const ;

    /** */
    void writeCellData(ostream& os, bool isPHeader) const ;

    /** */
    void writeDataArray(ostream& os, const string& name,
                        const RefCountPtr<FieldBase>& expr, bool isPHeader, bool isPointData) const ;
#endif /* DOXYGEN_DEVELOPER_ONLY */
  };
}




#endif
