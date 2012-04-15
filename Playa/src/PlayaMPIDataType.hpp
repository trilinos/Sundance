// @HEADER
// @HEADER

#ifndef PLAYA_MPI_DATA_TYPE_H
#define PLAYA_MPI_DATA_TYPE_H


#include "PlayaDefs.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include <stack>

#ifdef HAVE_MPI
#include "mpi.h"
#endif


namespace Playa
{
using Teuchos::RCP;
using Teuchos::rcp;

/** */
class MPIDataType
{
public:
  /** */
  MPIDataType(const std::string& name);

#ifdef HAVE_MPI
  /** */
  MPIDataType(const std::string& name, const RCP<MPI_Datatype>& mpiType);
#endif

  /** */
  const std::string& name() const {return name_;}

#ifdef HAVE_MPI
  /** */
  MPI_Datatype* ptr() ;

  /** */
  const MPI_Datatype& handle() const ;
#endif

  /** */
  static MPIDataType intType() ;
  
  /** */
  static MPIDataType floatType() ;
  
  /** */
  static MPIDataType doubleType() ;
  
  /** */
  static MPIDataType doubleIntPairType() ;
  
  /** */
  static MPIDataType charType() ;

  /** */
  static void registerType(const MPIDataType& dataType);
      
  static std::stack<MPIDataType>& typeRegistry();

  static void clearTypeRegistry();

private:
  std::string name_;
#ifdef HAVE_MPI
  RCP<MPI_Datatype> mpiType_;
#endif
};


	
} // namespace Playa

#endif
