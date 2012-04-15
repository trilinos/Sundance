// @HEADER
// @HEADER

#include "PlayaMPIOp.hpp"

namespace Playa
{
using Teuchos::RCP;
using Teuchos::rcp;

MPIOp::MPIOp(const std::string& name)
  : name_()
#ifdef HAVE_MPI
  , mpiOp_(rcp(new MPI_Op()))
#endif
{}

#ifdef HAVE_MPI
MPIOp::MPIOp(const std::string& name, 
  const RCP<MPI_Op>& mpiOp)
  : name_(name), mpiOp_(mpiOp){}


MPI_Op* MPIOp::ptr() 
{
  TEUCHOS_TEST_FOR_EXCEPT(mpiOp_.get()==0);
  return mpiOp_.get();
}

const MPI_Op& MPIOp::handle() const 
{
  TEUCHOS_TEST_FOR_EXCEPT(mpiOp_.get()==0);
  return *(mpiOp_.get());
}
#endif



std::stack<MPIOp>& MPIOp::opRegistry()
{
  static std::stack<MPIOp> rtn;

  return rtn;
}

void MPIOp::clearOpRegistry()
{
  while(!opRegistry().empty())
  {
#ifdef HAVE_MPI
    MPIOp t = opRegistry().top();

    int ierr = MPI_Op_free(t.ptr());

    TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::runtime_error,
      "Error code=" << ierr << " detected in MPI_Type_free()");
#endif
    opRegistry().pop();
  }
}

void MPIOp::registerOp(const MPIOp& opType)
{
  opRegistry().push(opType);
}


MPIOp MPIOp::sumOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI sum", rcp(new MPI_Op(MPI_SUM)));
#else
  static MPIOp rtn("MPI sum");
#endif
  return rtn;
}


MPIOp MPIOp::minOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI min", rcp(new MPI_Op(MPI_MIN)));
#else
  static MPIOp rtn("MPI min");
#endif
  return rtn;
}

MPIOp MPIOp::maxOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI max", rcp(new MPI_Op(MPI_MAX)));
#else
  static MPIOp rtn("MPI max");
#endif
  return rtn;
}

MPIOp MPIOp::minlocOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI minloc", rcp(new MPI_Op(MPI_MINLOC)));
#else
  static MPIOp rtn("MPI minloc");
#endif
  return rtn;
}

MPIOp MPIOp::maxlocOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI maxloc", rcp(new MPI_Op(MPI_MAXLOC)));
#else
  static MPIOp rtn("MPI maxloc");
#endif
  return rtn;
}

MPIOp MPIOp::productOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI prod", rcp(new MPI_Op(MPI_PROD)));
#else
  static MPIOp rtn("MPI prod");
#endif
  return rtn;
}


	
} // namespace Playa

