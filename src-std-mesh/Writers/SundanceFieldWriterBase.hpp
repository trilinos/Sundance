/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FIELDWRITERBASE_H
#define SUNDANCE_FIELDWRITERBASE_H


#ifndef DOXYGEN_DEVELOPER_ONLY


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceFieldBase.hpp"

namespace SundanceStdMesh
{
  namespace Internal
  {
    /**
     * FieldWriterBase is a base class for objects that write fields
     * and/or meshes to a stream. 
     */
    class FieldWriterBase : public TSFExtended::Handleable<FieldWriterBase>,
                            public TSFExtended::ObjectWithVerbosity<FieldWriterBase>
    {
    public:
      /** */
      FieldWriterBase(const string& filename);

      /** virtual dtor */
      virtual ~FieldWriterBase(){;}

      /** */
      void addMesh(const Mesh& mesh);

      /** add a comment */
      virtual void addCommentLine(const string& line) ;

      /** add a field, tagging it with the given string as a name */
      virtual void addField(const string& name, 
                            const RefCountPtr<FieldBase>& field) ;

      /** */
      virtual void write() const = 0 ;

      /**  */
      virtual void impersonateParallelProc(int nProc, int rank) ;

    protected:
      /** */
      int nProc() const ;

      /** */
      int myRank() const ;

      /** */
      const string& filename() const {return filename_;}

      /** */
      const Mesh& mesh() const {return mesh_;}

      /** Indicate whether the given writer subtype does anything special
       * for vector fields. Default is false, in which case
       * vectors are simply written as a list of scalars.
       */
      virtual bool supportsSpecializedVectors() const {return false;}

      const Array<string>& comments() const {return comments_;}
      Array<string>& comments() {return comments_;}

      const Array<RefCountPtr<FieldBase> >& pointScalarFields() const {return pointScalarFields_;}
      Array<RefCountPtr<FieldBase> >& pointScalarFields() {return pointScalarFields_;}

      const Array<RefCountPtr<FieldBase> >& cellScalarFields() const {return cellScalarFields_;}
      Array<RefCountPtr<FieldBase> >& cellScalarFields() {return cellScalarFields_;}

      const Array<string>& pointScalarNames() const {return pointScalarNames_;}
      Array<string>& pointScalarNames() {return pointScalarNames_;}

      const Array<string>& cellScalarNames() const {return cellScalarNames_;}
      Array<string>& cellScalarNames() {return cellScalarNames_;}

      const Array<RefCountPtr<FieldBase> >& pointVectorFields() const {return pointVectorFields_;}
      Array<RefCountPtr<FieldBase> >& pointVectorFields() {return pointVectorFields_;}

      const Array<RefCountPtr<FieldBase> >& cellVectorFields() const {return cellVectorFields_;}
      Array<RefCountPtr<FieldBase> >& cellVectorFields() {return cellVectorFields_;}

      const Array<string>& pointVectorNames() const {return pointVectorNames_;}
      Array<string>& pointVectorNames() {return pointVectorNames_;}

      const Array<string>& cellVectorNames() const {return cellVectorNames_;}
      Array<string>& cellVectorNames() {return cellVectorNames_;}

      virtual void writeCommentLine(const string& line) const {;}

    private:
      string filename_;

      Mesh mesh_;

      int nProc_;

      int myRank_;

      int meshID_;

      Array<string> comments_;

      Array<RefCountPtr<FieldBase> > pointScalarFields_;
      Array<RefCountPtr<FieldBase> > cellScalarFields_;
      Array<RefCountPtr<FieldBase> > pointVectorFields_;
      Array<RefCountPtr<FieldBase> > cellVectorFields_;
      Array<string> pointScalarNames_;
      Array<string> cellScalarNames_;
      Array<string> pointVectorNames_;
      Array<string> cellVectorNames_;
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif