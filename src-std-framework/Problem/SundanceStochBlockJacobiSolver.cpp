#include "SundanceStochBlockJacobiSolver.hpp"
#include "Sundance.hpp"

void
StochBlockJacobiSolver::solve(const Array<LinearOperator<double> >& KBlock,
  const Array<Vector<double> >& fBlock,
  Array<Vector<double> >& xBlock) const
{
  int L = KBlock.size();
  int P = pcBasis_.nterms();
  int Q = fBlock.size();

  /*
   * Solve the equations using block Gauss-Jacobi iteration
   */
  Array<Vector<double> > uPrev(P);
  Array<Vector<double> > uCur(P);

  for (int i=0; i<P; i++)
  {
    uPrev[i] = fBlock[0].copy();
    uCur[i] = fBlock[0].copy();
    uPrev[i].zero();
    uCur[i].zero();
  }

  if (verbosity_) Out::root() << "starting Jacobi loop" << endl;
  bool converged = false;
  for (int iter=0; iter<maxIters_; iter++)
  {
    if (verbosity_) Out::root() << "Jacobi iter=" << iter << endl;
    bool haveNonConvergedBlock = false;
    double maxErr = -1.0;
    int numNonzeroBlocks = 0;
    for (int i=0; i<P; i++)
    {
      if (verbosity_) Out::root() << "Iter " << iter << ": block row i=" << i << " of " << P << " ..." << ends;
      Vector<double> b = fBlock[0].copy();
      b.zero();
      int nVecAdds = 0;
      for (int j=0; j<Q; j++)
      {
        double c_ij0 = pcBasis_.expectation(i,j,0);
        if (fabs(c_ij0) > 0.0) 
        {
          b = b + c_ij0 * fBlock[j];
          nVecAdds++;
        }
        if (j>=L) continue; 
        Vector<double> tmp = fBlock[0].copy();
        tmp.zero();
        bool blockIsNeeded = false;
        for (int k=0; k<P; k++)
        {
          if (j==0 && k==i) continue;
          double c_ijk = pcBasis_.expectation(i,j,k);
          if (fabs(c_ijk) > 0.0)
          {
            tmp = tmp + c_ijk * uPrev[k];
            nVecAdds++;
            blockIsNeeded = true;
          }
        }
        numNonzeroBlocks += blockIsNeeded;
        b = (b - KBlock[j]*tmp);
        nVecAdds++;
      }
      b = b * (1.0/pcBasis_.expectation(i,i,0));
      if (verbosity_) Out::root() << "num vec adds = " << nVecAdds << endl;
      diagonalSolver_.solve(KBlock[0], b, uCur[i]);
      double err = (uCur[i]-uPrev[i]).norm2();
      if (err > convTol_) haveNonConvergedBlock=true;
      if (err > maxErr) maxErr = err;
    }

    /* update solution blocks */
    for (int i=0; i<P; i++) uPrev[i] = uCur[i].copy();
      
    /* done all block rows -- check convergence */
    if (!haveNonConvergedBlock)
    {
      if (verbosity_) Out::root() << "=======> max err=" << maxErr << endl;
      if (verbosity_) Out::root() << "=======> converged! woo-hoo!" << endl;
      if (verbosity_) Out::root() << "estimated storage cost: " 
                  << setprecision(3) << 100*((double) L)/((double) numNonzeroBlocks) 
                  << " percent of monolithic storage" << endl;
      converged = true;
      break;
    }
    else
    {
      if (verbosity_) Out::root() << "maxErr=" << maxErr << ", trying again" << endl;
    }
  }

  TEST_FOR_EXCEPT(!converged);
  xBlock = uCur;
}
