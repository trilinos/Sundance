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

#ifndef SUNDANCE_EVALVECTOR_H
#define SUNDANCE_EVALVECTOR_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceNoncopyable.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    using namespace Teuchos;

    /**
     *
     */
    class EvalVector : public Noncopyable,
                       public TSFExtended::ObjectWithVerbosity<EvalVector>
    {
      friend class EvalManager;
      friend class TempStack;

    private:

      /** */
      EvalVector(TempStack* s);

      /** */
      EvalVector(TempStack* s, const RefCountPtr<Array<double> >& data,
                 const string& str);


    public:
      /** 
       * EvalVector has a nontrivial destructor. Upon destruction, 
       * the vector's underlying data object is not destroyed, but rather
       * is put back on the stack of temporary vectors. 
       */
      ~EvalVector();

      /** \name Mathematical operations */
      //@{

      /** */
      void add_SV(const double& alpha, 
                  const EvalVector* B) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this + alpha*B*C
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void add_SVV(const double& alpha,
                   const EvalVector* B,
                   const EvalVector* C) ;

      /** */
      void add_V(const EvalVector* A) ;

      /** */
      void add_S(const double& alpha);

      /**
       * Perform the operation 
       * \f[ 
       * this = this + A*B
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void add_VV(const EvalVector* A,
                  const EvalVector* B) ;


      /**  
       * Perform a scaled addition with another vector,
       * \f[ 
       * this = \alpha this + \beta C
       * \f]
       * The operation is done in-place, overwriting the old values of the
       * vector. 
       */
      void multiply_S_add_SV(const double& alpha, 
                             const double& beta,
                             const EvalVector* C) ;

      /** Scale and add a constant to this vector. 
       * The operation is done in-place, overwriting the old values of
       * the vector. Each element x[i] is updated as:
       * \f[
       * this = alpha * this + beta
       * \f]
       */
      void multiply_S_add_S(const double& alpha,
                            const double& beta) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*A + B*C*D
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_V_add_VVV(const EvalVector* A,
                              const EvalVector* B,
                              const EvalVector* C,
                              const EvalVector* D) ;


      /**
       * Perform the operation 
       * \f[ 
       * this = this*A + beta*C*D
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_V_add_SVV(const EvalVector* A,
                              const double& beta,
                              const EvalVector* C,
                              const EvalVector* D) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*A + beta*C
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_V_add_SV(const EvalVector* A,
                             const double& beta,
                             const EvalVector* C) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*A*B
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_VV(const EvalVector* A,
                       const EvalVector* B) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*alpha*B
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_SV(const double& alpha,
                       const EvalVector* B) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*A
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_V(const EvalVector* A) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*alpha
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_S(const double& alpha) ;

      /**
       *
       */
      void setTo_S_add_SVV(const double& alpha,
                           const double& beta,
                           const EvalVector* C,
                           const EvalVector* D);

      /**
       *
       */
      void setTo_S_add_VV(const double& alpha,
                          const EvalVector* B,
                          const EvalVector* C);

      /**
       *
       */
      void setTo_S_add_SV(const double& alpha,
                          const double& beta,
                          const EvalVector* C);

      /** 
       *
       */
      void setTo_S_add_V(const double& alpha,
                         const EvalVector* B);


      /**
       *
       */
      void setTo_V(const EvalVector* A);

      /**
       *
       */
      void setTo_VV(const EvalVector* A,
                    const EvalVector* B);

      /**
       *
       */
      void setTo_SV(const double& alpha,
                    const EvalVector* B);

      /**
       *
       */
      void setTo_SVV(const double& alpha,
                     const EvalVector* B,
                     const EvalVector* C);

      


      /**
       * Set every element to a constant value
       */
      void setToConstant(const double& alpha) ;

      /** 
       * Apply a unary function
       */
      void applyUnaryOperator(const UnaryFunctor* func, 
                              Array<RefCountPtr<EvalVector> >& opDerivs);
      
      
      /** */
      RefCountPtr<EvalVector> clone() const ;

      /** */
      void resize(int n);

      /** */
      int length() const {return data_->size();}
      
      /** */
      void print(ostream& os) const ;

      /** */
      const double * const start() const {return &((*data_)[0]);}

      /** */
      double * const start() {return &((*data_)[0]);}

      const string& str() const {return str_;}

      void setString(const string& str) {str_ = str;}

      inline static bool& shadowOps() {static bool rtn = false; return rtn;}

      bool isValid() const {return data_.get() != 0 && s_ != 0;}
      //@}

      

      inline static double& totalFlops() {static double rtn = 0; return rtn;}

    private:

      inline static void addFlops(const double& flops) {totalFlops() += flops;}

      mutable TempStack* s_;

      RefCountPtr<Array<double> > data_;

      string str_;

    };




    //     inline void EvalVector::setToConstant(const double& alpha) 
    //     {
    //       int n = data_->size();
    //       if (n > 0)
    //         {
    //           double* x = &((*data_)[0]);
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] = alpha;
    //             }
    //         }
    //       if (shadowOps()) str_ = Teuchos::toString(alpha);
    //     }


    //     inline void EvalVector::add_SV(const double& alpha,
    //                             const EvalVector* B)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Bx = B->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] += alpha*Bx[i];
    //             }

    //           addFlops(2*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + str_ + "+" 
    //             + Teuchos::toString(alpha) + "*" + B->str_ + ")";
    //         }
    //     }

    //     inline void EvalVector::add_S(const double& alpha)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] += alpha;
    //             }
    //           addFlops(n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + str_ + "+" 
    //             + Teuchos::toString(alpha) + ")";
    //         }
    //     }


    //     inline void EvalVector::add_V(const EvalVector* A)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Ax = A->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] += Ax[i];
    //             }
    //           addFlops(n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + str_ + " + " +  A->str_ + ")";
    //         }
    //     }

    //     inline void EvalVector::add_SVV(const double& alpha,
    //                              const EvalVector* B,
    //                              const EvalVector* C)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Bx = B->start();
    //           const double* const Cx = C->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] += alpha*Bx[i]*Cx[i];
    //             }
    //           addFlops(3*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + str_ + " + " + Teuchos::toString(alpha) + "*"
    //             + B->str() + "*" + C->str() + ")";
    //         }
    //     }

    //     inline void EvalVector::add_VV(const EvalVector* A,
    //                             const EvalVector* B)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Ax = A->start();
    //           const double* const Bx = B->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] += Ax[i]*Bx[i];
    //             }
    //           addFlops(2*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + str_ + " + " + A->str() + "*" + B->str() + ")";
    //         }
    //     }


    //     inline void EvalVector::multiply_S_add_SV(const double& alpha, 
    //                                        const double& beta,
    //                                        const EvalVector* C)
    //     {
    //       int n = data_->size();
  
    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Cx = C->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] *= alpha;
    //               x[i] += beta*Cx[i];
    //             }
    //           addFlops(3*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + Teuchos::toString(alpha) + "*" + str_ + "+"
    //             + Teuchos::toString(beta) + "*" + C->str_ + ")";
    //         }
    //     }


    //     inline void EvalVector::multiply_S_add_S(const double& alpha, 
    //                                       const double& beta)
    //     {
    //       int n = data_->size();
  
    //       if (n > 0)
    //         {
    //           double* const x = start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] *= alpha;
    //               x[i] += beta;
    //             }
    //           addFlops(2*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + Teuchos::toString(alpha) + "*" + str_
    //             + " + " + Teuchos::toString(beta) + ")";
    //         }
    //     }

    //     inline void EvalVector::multiply_V(const EvalVector* A)
    //     {
    //       int n = data_->size();
  
    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Ax = A->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] *= Ax[i];
    //             }
    //           addFlops(n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = str() + "*" + A->str();
    //         }
    //     }

    //     inline void EvalVector::multiply_V_add_VVV(const EvalVector* A,
    //                                         const EvalVector* B,
    //                                         const EvalVector* C,
    //                                         const EvalVector* D)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Ax = A->start();
    //           const double* const Bx = B->start();
    //           const double* const Cx = C->start();
    //           const double* const Dx = D->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] *= Ax[i];
    //               x[i] += Bx[i]*Cx[i]*Dx[i];
    //             }
    //           addFlops(4*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + str() + "*" + A->str() + " + " 
    //             + B->str() + "*" + C->str() + "*" + D->str() + ")";
    //         }
    //     }

    //     inline void EvalVector::multiply_V_add_SVV(const EvalVector* A,
    //                                         const double& beta,
    //                                         const EvalVector* C,
    //                                         const EvalVector* D)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Ax = A->start();
    //           const double* const Cx = C->start();
    //           const double* const Dx = D->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] *= Ax[i];
    //               x[i] += beta*Cx[i]*Dx[i];
    //             }
    //           addFlops(4*n);
    //         }

    //       if (shadowOps())
    //         {
    //           if (beta != 0.0)
    //             {
    //               str_ = "(" + str() + "*" + A->str();
    //             }
    //           else
    //             {
    //               str_ = str() + "*" + A->str();
    //             }
    //           if (beta != 0.0)
    //             {
    //               str_ += " + ";
    //               if (beta != 1.0)
    //                 {
    //                   str_ += Teuchos::toString(beta) + "*";
    //                 }
    //               str_ += C->str() + "*" + D->str();
    //             }
    //         }
    //     }

    //     inline void EvalVector::multiply_V_add_SV(const EvalVector* A,
    //                                        const double& beta,
    //                                        const EvalVector* C)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Ax = A->start();
    //           const double* const Cx = C->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] *= Ax[i];
    //               x[i] += beta*Cx[i];
    //             }
    //           addFlops(3*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + str() + "*" + A->str() + " + " + Teuchos::toString(beta)
    //             + "*" + C->str() + ")";
    //         }
    //     }

    //     inline void EvalVector::multiply_VV(const EvalVector* A,
    //                                  const EvalVector* B)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Ax = A->start();
    //           const double* const Bx = B->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] *= Ax[i]*Bx[i];
    //             }
    //           addFlops(2*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = str() + "*" + A->str() + "*" + B->str();
    //         }
    //     }


    //     inline void EvalVector::multiply_SV(const double& alpha,
    //                                  const EvalVector* B)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
      
    //           const double* const Bx = B->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] *= alpha*Bx[i];
    //             }
    //           addFlops(2*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = Teuchos::toString(alpha) + "*" + str_ + "*" + B->str();
    //         }
    //     }

    //     inline void EvalVector::multiply_S(const double& alpha)
    //     {
    //       int n = data_->size();

    //       if (n > 0)
    //         {
    //           double* const x = start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] *= alpha;
    //             }
    //           addFlops(n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = Teuchos::toString(alpha) + "*" + str_;
    //         }
    //     }


    //     inline void EvalVector::setTo_S_add_SVV(const double& alpha,
    //                                      const double& beta,
    //                                      const EvalVector* C,
    //                                      const EvalVector* D)
    //     {
    //       int n = C->data_->size();
    //       resize(n);

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Cx = C->start();
    //           const double* const Dx = D->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] = alpha + beta*Cx[i]*Dx[i];
    //             }
    //           addFlops(3*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + Teuchos::toString(alpha) + " + " 
    //             + Teuchos::toString(beta) + "*" 
    //             + C->str() + "*" + D->str() + ")";
    //         }
    //     }

    //     inline void EvalVector::setTo_S_add_VV(const double& alpha,
    //                                     const EvalVector* B,
    //                                     const EvalVector* C)
    //     {
    //       int n = B->data_->size();
    //       resize(n);

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Bx = B->start();
    //           const double* const Cx = C->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] = alpha + Bx[i]*Cx[i];
    //             }
    //           addFlops(2*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + Teuchos::toString(alpha) + " + " 
    //             + B->str() + "*" + C->str() + ")";
    //         }
    //     }

    //     inline void EvalVector::setTo_S_add_SV(const double& alpha,
    //                                     const double& beta,
    //                                     const EvalVector* C)
    //     {
    //       int n = C->data_->size();
    //       resize(n);

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Cx = C->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] = alpha + beta*Cx[i];
    //             }
    //           addFlops(2*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + Teuchos::toString(alpha) + " + " 
    //             + Teuchos::toString(beta) + "*" 
    //             + C->str() + ")";
    //         }
    //     }



    //     inline void EvalVector::setTo_S_add_V(const double& alpha,
    //                                    const EvalVector* B)
    //     {
    //       int n = B->data_->size();
    //       resize(n);

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Bx = B->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] = alpha + Bx[i];
    //             }
    //           addFlops(n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = "(" + Teuchos::toString(alpha) + " + " 
    //             + B->str() + ")";
    //         }
    //     }


    //     inline void EvalVector::setTo_V(const EvalVector* A)
    //     {
    //       int n = A->data_->size();
    //       resize(n);

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Ax = A->start();
      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] = Ax[i];
    //             }
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = A->str();
    //         }
    //     }



    //     inline void EvalVector::setTo_SVV(const double& alpha,
    //                                const EvalVector* B,
    //                                const EvalVector* C)
    //     {
    //       int n = B->data_->size();
    //       resize(n);

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Bx = B->start();
    //           const double* const Cx = C->start();

      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] = alpha*Bx[i]*Cx[i];
    //             }
    //           addFlops(2*n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = Teuchos::toString(alpha) + "*" 
    //             + B->str() + "*" + C->str();
    //         }
    //     }

    //     inline void EvalVector::setTo_VV(const EvalVector* A,
    //                               const EvalVector* B)
    //     {
    //       int n = A->data_->size();
    //       resize(n);

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Ax = A->start();
    //           const double* const Bx = B->start();

      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] = Ax[i]*Bx[i];
    //             }

    //           addFlops(n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = A->str() + "*" + B->str();
    //         }
    //     }

    //     inline void EvalVector::setTo_SV(const double& alpha,
    //                               const EvalVector* B)
    //     {
    //       int n = B->data_->size();
    //       resize(n);

    //       if (n > 0)
    //         {
    //           double* const x = start();
    //           const double* const Bx = B->start();

      
    //           for (int i=0; i<n; i++)
    //             {
    //               x[i] = alpha*Bx[i];
    //             }
    //           addFlops(n);
    //         }

    //       if (shadowOps())
    //         {
    //           str_ = Teuchos::toString(alpha) + "*" 
    //             + B->str();
    //         }
    //     }
  }


}

namespace std
{
  inline ostream& operator<<(ostream& os, 
                             const SundanceCore::Internal::EvalVector& vec)
  {
    vec.print(os);
    return os;
  }
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif
