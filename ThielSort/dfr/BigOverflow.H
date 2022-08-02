/* $Header:  $ */ 

#ifndef BigOverflowH
#define BigOverflowH

#include "BktRemapper.H"

/********************************************************************
 *
 *                            BigOverflow.H
 *
 * Data structures and functions used to squeeze almost all of the
 * (very large) overflow areas back into the primary bkt buffer,
 * avoiding the need to allocate expensively large overflow buffers.
 *
 * Used by Overflow.H.
 ********************************************************************/

// BigOverflow                                                  BigOverflow
template <typename Elem_T, int N_Buckets, int Rt_Per_Mu>
class BigOverflow {
 public:

  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  //BigOverflow();
  //~BigOverflow();

  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  // Transfer functions to the imp code in BigOverflow.C.
  /**
   * Populate Bkt[] from the vectors of bkt counts and compute
   * TotalPN, TotalOfN, and ReqdPadding.
   **/
  void InitBktVec(const CountType *p_bkt_ns, const CountType *of_bkt_ns) {
    InitBktVec_imp(p_bkt_ns, of_bkt_ns);
  }
  /**
   * Determine which chunks will go to the OF area and which will go to
   * the primary area.  Create both area maps by assigning values to the
   * offset for the chunks to be stored there. Return the required size
   * of the OF area.
   **/
  CountType CreateAreaMaps() { return CreateAreaMaps_imp(); }

  /**
   * On entry, p_addr_v[] contains the original addrs of the
   * pbuckets. The bkt lens are already stored in Bkt[]. This function
   * will move each pbkt to its new location, store the remapped addrs
   * of each pbkt in p_addr_v and each ofbkt into of_addr_v. p_bfr_p
   * is the start of the primary area and of_bfr_p is the start of the
   * OF area. If a chunk is empty, nullptr is stored.
   **/
  void SetOfChunkAddrs(Elem_T** p_addr_v, Elem_T** of_addr_v,
		       Elem_T* p_bfr_p, Elem_T* of_bfr_p) {
    
    for (int b = 0; b < MaxNBuckets; b++) {
      Bkt[b].InitOldAddr(p_addr_v[b]); 
      p_addr_v[b] = Bkt[b].PBkt(p_bfr_p, of_bfr_p);
      of_addr_v[b] = Bkt[b].OfBkt(p_bfr_p, of_bfr_p);
    }
    // Addrs are set up. Now move the pbkts. First move all pbkts
    // which don't have to be moved (normally, only the first, and
    // empty ones) or for which the target space is in the OF area,
    // since is always available.
    MoveEasyBkts(p_addr_v);
    ShiftBktsLeft(p_addr_v);
    ShiftBktsRight(p_addr_v);
  }
private:

  // Allocate a sentinel at either end of this vector. Address Via Bkt
  BktRemapper<Elem_T> Bkt[MaxNBuckets];
  
  // Total reqd space, including padding.
  CountType TotalPN;
  CountType TotalOfN;
  // ReqdPadding is computed by InitBktVec to reflect the actual padding
  // required for each bucket (of both kinds).
  int ReqdPadding;  // Max total padding needed for all chunks.
  
#include "BigOverflow.C"
  
};  // BigOverflow

#endif

