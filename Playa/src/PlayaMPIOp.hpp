// @HEADER
// @HEADER

#ifndef PLAYA_MPI_OP_H
#define PLAYA_MPI_OP_H


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
class MPIOp
{
public:
  /** */
  MPIOp(const std::string& name);

#ifdef HAVE_MPI
  /** */
  MPIOp(const std::string& name, const RCP<MPI_Op>& mpiType);
#endif

  /** */
  const std::string& name() const {return name_;}

#ifdef HAVE_MPI
  MPI_Op* ptr() ;
#endif

#ifdef HAVE_MPI
  const MPI_Op& handle() const ;
#endif


  /** */
  static MPIOp sumOp() ;
  
  /** */
  static MPIOp minOp() ;
  
  /** */
  static MPIOp maxOp() ;
  
  /** */
  static MPIOp minlocOp() ;
  
  /** */
  static MPIOp maxlocOp() ;
  
  /** */
  static MPIOp productOp() ;

  /** */
  static void registerOp(const MPIOp& op);
      
  static std::stack<MPIOp>& opRegistry();

  static void clearOpRegistry();
private:
  std::string name_;
#ifdef HAVE_MPI
  RCP<MPI_Op> mpiOp_;
#endif
};


	
} // namespace Playa

#endif
