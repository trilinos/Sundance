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

#ifndef SUNDANCE_FIELDWRITERBASE_H
#define SUNDANCE_FIELDWRITERBASE_H


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

      /** set the numerical value to be written at cells on which
       * a field is undefined. */
      void setUndefinedValue(const double& x) {undefinedValue_ = x;}

    protected:

      /** */
      double undefinedValue() const {return undefinedValue_;}
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

      double undefinedValue_;
    };
  }
}



#endif
