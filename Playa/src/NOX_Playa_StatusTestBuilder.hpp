// $Id$ 
// $Source$ 

//@HEADER
//   
//@HEADER

#ifndef NOX_Playa_STATUSTESTBUILDER_H
#define NOX_Playa_STATUSTESTBUILDER_H

#include "NOX_Abstract_Group.H"	
#include "NOX_StatusTest_Generic.H"	
#include "NOX_Common.H"         
#include "Teuchos_RefCountPtr.hpp"         
#include "Teuchos_ParameterList.hpp"         



namespace NOX 
{
  namespace NOXPlaya 
  {
    class StatusTestBuilder
    {
    public:
      static Teuchos::RCP<NOX::StatusTest::Generic> 
      makeStatusTest(const Teuchos::ParameterList& params) ;

    };

  } // namespace Playa
} // namespace NOX


#endif
