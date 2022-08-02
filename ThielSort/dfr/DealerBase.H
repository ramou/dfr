/* $Header:  $ */ 

#ifndef DealerBaseH
#define DealerBaseH

/********************************************************************
 *
 *                            DealerBase.H
 *
 * DealerBase provides access to virtual dealing functions.
 ********************************************************************/

// DealerBase                                                     DealerBase
template<typename ELEM> class DealerBase {

 public:
  // Functions to support dealing from a bucket vector DataPool
  virtual void ProcessBlk(const ELEM *e_vec, int n) = 0;

};  // DealerBase

#endif
