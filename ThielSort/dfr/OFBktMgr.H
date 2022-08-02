/* $Header:  $ */ 

#ifndef OFBktMgrH
#define OFBktMgrH

/********************************************************************
 *
 *                            OFBktMgr.H
 *
 * Create and manage the list of bkts which overflowed this pass. Uses
 * the store_n and store_fence vectors to locate and count bkts
 * which are overflowing. Allows ExtractNext() to return the next OF bkt
 * (in some undefined order) and TestAndClear() to test and unmark OF
 * bkts. Uses simd ops to help identify bkts. 
 *
 * This class doesn't know about the record type. It just uses
 * unsigned CountType / simd overlay vectors to record and manipulate
 * info about bkts. We use unsigned so we can use a logical right
 * shift to get a 0 or 1 result when subtracting store_n[b] from
 * store_fence[b] in MarkOfBkts().
 *
 * The simd overlay template SimdOverlay<> is declared here as public
 * and can be used by the HybridFlusher.
 ********************************************************************/

// OFBktMgr                                                   OFBktMgr
class OFBktMgr {
public:
  using BaseT = /*unsigned*/ CountType;
  static constexpr int ShiftCnt = (sizeof(BaseT) * 8) -1;
  
  /**
   * BaseT vectors overlaid with simd vectors allowing manipulation by
   * simd instructions. For indices and counts describing buckets and
   * ladles. Template parameters:
   * 
   * Base_N - Number of BaseT elements in the vector. If 0, Base_N
   *   is set to 2 * NPerSimd (Allowing a 1 simd wide summary and a
   *   sentinel).
   *
   * Log_Su_Len - Log of simd size in bytes: 4 means 16 byte simds, 5
   *   means 32, and 6 means 64. Defaults to machine specific info
   *   from makefile.
   **/
  // Simd overlay definitions for int indices and counts in buckets.
  template<int Base_N, int Log_Su_Len = SIMD_SHIFT_SIZE>
  struct SimdOverlay {
    using SimdT = SimdSupport::VectorType<BaseT, Log_Su_Len>;

    static constexpr int LogSuLen = Log_Su_Len;
    static constexpr int NPerSimd = SimdT::NPerSimd;
    // size of base vector
    static constexpr int BSize = Base_N ? Base_N : 2 * NPerSimd;
    // And size of simd vector.
    static constexpr int SdSize = (Base_N + NPerSimd - 1) / NPerSimd;
    union {
      BaseT BV[BSize];  // Base Vector
      SimdT SdV[SdSize];  // Simd vector
    };
    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor

    SimdOverlay() {}
    SimdOverlay(BaseT val) {
      for (int i = 0; i < BSize; i++) BV[i] = val;
    }
  };
  using PerBktIntT = SimdOverlay<MaxNBuckets>;
  using OneSimdIntT = SimdOverlay<0>;
  using SimdT = PerBktIntT::SimdT;
  
  static constexpr int BSize = PerBktIntT::BSize;
  static constexpr int SdSize = PerBktIntT::SdSize;
  static constexpr int NPerSimd = PerBktIntT::NPerSimd;
  
  OFBktMgr() {  }
  /**
   * Uses store_n and store_fence to identify the bkts where
   * store_next[b] > store_fence[b] and returns the number of bkts
   * where this happens. Prepares a list of marked bkts for scanning
   * via ExtractNext(). The marked bkt IDs will be returned 1 at a
   * time. The order in which they are returned is undefined.
   **/  
  int MarkOfBkts(PerBktIntT *store_n, PerBktIntT *store_fence) {
    if (USE_SIMD_PREP) {
      return MarkOfBkts_simd_imp (store_n, store_fence);
    } else {
      return MarkOfBkts_imp(store_n, store_fence);
    }
  }
  // Version when looking for NAvail < 0.
  int MarkOfBkts(PerBktIntT *n_avail) {
    return MarkOfBkts_imp(n_avail);
  }
  /**
   * If bkt b is marked, unmark it and return true. Otherwise, return
   * false. 
   *
   * Use if b might be flushable to see if it is actually
   * overflowing. The unmark means that after checking the bkts known
   * to be potentially flushable, any marked bkts (findable with
   * ExtractNext()) are new overflows.
   **/  
  bool TestAndClear(int b) {
    bool t = TagVec.BV[b];
    if (t) {
      StripeSums.BV[StripeNdx(b)]--;
      TagVec.BV[b] = 0;
    }
    return t;
  }
  /**
   * Reset the list for scanning and return the first bkt id in the
   * list (which may not be the lowest numbered bkt).
   **/
  int ExtractNext() {
    return ResetList_imp();
  }
  /**
   * Return the bkt id after prev_b in the list (and unmark it) or -1
   * if there are no more.
   **/
  int ExtractNext(int prev_b) { return ExtractNext_imp(prev_b); }
private:

  PerBktIntT TagVec;
  // Simd ops only use first simd of StripeSums, but we want a
  // sentinel after its end.
  OneSimdIntT StripeSums{1};  // Set sentinel in ~[1]. ~[0] doesn't care.

#include "OFBktMgr.c"
} __attribute__ ((aligned (CACHE_LINE_LEN)));  // OFBktMgr

#endif
