// @HEADER
// @HEADER

#include "PlayaMPIDataType.hpp"

namespace Playa
{
using Teuchos::RCP;
using Teuchos::rcp;

MPIDataType::MPIDataType(const std::string& name)
  : name_()
#ifdef HAVE_MPI
  , mpiType_(rcp(new MPI_Datatype()))
#endif
{}

#ifdef HAVE_MPI
MPIDataType::MPIDataType(const std::string& name, 
  const RCP<MPI_Datatype>& mpiType)
  : name_(name), mpiType_(mpiType){}


MPI_Datatype* MPIDataType::ptr() 
{
  TEUCHOS_TEST_FOR_EXCEPT(mpiType_.get()==0);
  return mpiType_.get();
}

const MPI_Datatype& MPIDataType::handle() const 
{
  TEUCHOS_TEST_FOR_EXCEPT(mpiType_.get()==0);
  return *(mpiType_.get());
}
#endif



std::stack<MPIDataType>& MPIDataType::typeRegistry()
{
  static std::stack<MPIDataType> rtn;

  return rtn;
}

void MPIDataType::clearTypeRegistry()
{
  while(!typeRegistry().empty())
  {
#ifdef HAVE_MPI
    MPIDataType t = typeRegistry().top();

    int ierr = MPI_Type_free(t.ptr());

    TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::runtime_error,
      "Error code=" << ierr << " detected in MPI_Type_free()");
#endif
    typeRegistry().pop();
  }
}

void MPIDataType::registerType(const MPIDataType& dataType)
{
  typeRegistry().push(dataType);
}

MPIDataType MPIDataType::intType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI int", rcp(new MPI_Datatype(MPI_INT)));
#else
  static MPIDataType rtn("MPI int");
#endif
  return rtn;
}



MPIDataType MPIDataType::floatType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI float", rcp(new MPI_Datatype(MPI_FLOAT)));
#else
  static MPIDataType rtn("MPI float");
#endif
  return rtn;
}



MPIDataType MPIDataType::doubleType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI double", rcp(new MPI_Datatype(MPI_DOUBLE)));
#else
  static MPIDataType rtn("MPI double");
#endif
  return rtn;
}



MPIDataType MPIDataType::doubleIntPairType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI double/int pair", rcp(new MPI_Datatype(MPI_DOUBLE_INT)));
#else
  static MPIDataType rtn("MPI double/int pair");
#endif
  return rtn;
}




MPIDataType MPIDataType::charType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI char", rcp(new MPI_Datatype(MPI_CHAR)));
#else
  static MPIDataType rtn("MPI char");
#endif
  return rtn;
}



	
} // namespace Playa

