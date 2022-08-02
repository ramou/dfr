/* $Header:  $ */ 

#ifndef IsConstBaseH
#define IsConstBaseH

/********************************************************************
 *
 *                            IsConstBase.H
 *
 * Base class for insertion sort variations. The functions declared
 * here are just utility and search functions that mostly do not
 * modify the vector to be sorted. The exceptions are
 * BaseInit() and SetupTopSentinel(), which have limited
 * scope and will modify the original source at the ends before the
 * main sort starts. See their comments.
 *
 * Thus IsConstBase can be the base class for both the inplace IS and
 * the copy IS.
 ********************************************************************/

// IsConstBase                                                IsConstBase
template<typename Elem_T>class IsConstBase {
protected:

  static const int Nf = 8;
  // Thresholds control which sort is used for "large" SeBlock.
  static const int IsThreshold = 16;
  static const int FrThreshold = 4096;
  
  Elem_T TmMask;
  Elem_T SmallestRec;  // Bottom sentinel. But dups may exist.
  
  Elem_T TopMask(Elem_T val) { return val & TmMask; }
  Elem_T TopMask(const Elem_T *p) { return TmMask & *p; }

  static const int CntCats = 5;
  long Limits[CntCats] =     {ISNS_SHORT_RUN, Nf, 20, 50, 250};
  long Counts[CntCats + 1] = {0,  0,  0,   0,  0, 0};
  long Sums[CntCats + 1] = {0,  0,  0,   0,  0, 0};

  // Total number of buckets (from TmMask). For comparing with N.
  long TotalNBkts = 0;
  long TotalN = 0;
  // If Acpb > AcpbThreshold, use LRBracketedSort() (LR means
  // LongRun). Else use DoBracketedSort() (optimized for mostly fully
  // sorted).
  //double Acpb = 0.0;  // Average count per bucket.
  //double AcpbThreshold = ACPB_THRESHOLD;
  
  // For debugging.
  // DestVec is the addr of the original dest vect.
  struct Between {
    Elem_T* DestVec;
    Elem_T* TargetAdr;
    long TargetNdx;
  
    void InitDestOrigin(Elem_T* dest_v, long target_ndx) {
      DestVec = dest_v;
      TargetNdx = target_ndx;
      TargetAdr = DestVec + target_ndx;
    }
    
    // Call ErrTag iff (Target - range_pm) < p < (Target + range_pm). 
    void IsNear(Elem_T* p, int range_pm) {
      if (((TargetAdr - range_pm) < p) && (p < TargetAdr + range_pm))
	ErrTag("Between hit");
    }
  };

  //#define IS_FINDING_INDEX
#ifdef IS_FINDING_INDEX
  Between DbgNdxFinder;

  // Will call ErrTag() when a test address is between
  // (DestVec+TargetNdx - range_pm) and (DestVec+TargetNdx + range_pm).
  void IsNear(Elem_T* p, int range_pm) { DbgNdxFinder.IsNear(p, range_pm); }
#else
  // Any stray calls will be optimized out.
  void IsNear(Elem_T* p, int range_pm) {}
#endif  
  
  void ResetCounts() {
    for (int i = 0; i <= CntCats; i++) Counts[i] = Sums[i] = 0;
  }
  
  void LogCnts(int n) {
#ifdef IS_COUNTING
    static_assert(CntCats == 5);
    if (n <= Limits[1]) { Counts[(n > Limits[0])    ]++; return; }
    if (n <= Limits[3]) { Counts[(n > Limits[2]) + 2]++; return; }
    Counts[(n > Limits[4]) + 4]++;
#endif
  }

  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  /**
   * Store tm_mask for later use. Sort first (or only) SEB to put
   * bottom sentinel in place (and in SmallestRec). Return addr of
   * first unsorted record (which may be src_e if the entire sort is 1
   * SEB).
   **/  
  Elem_T *BaseInit(Elem_T *src_p, Elem_T *src_e,
		   Elem_T tm_mask, Elem_T *dest_v) {
    static_assert(Nf < IsThreshold, "Nf too big");
    TotalN = src_e - src_p;
    TmMask = tm_mask;
    //FindTtlNBkts(tm_mask);
    ResetCounts();
    if (TotalN == 0) return src_e;  // ???
    if (TopMask(src_p) != TopMask(src_e[-1])) {  // Not all 1 SEB.
      src_e = FindNextSeb<1>(src_p, src_e);
      SortLargeSeBlock(src_p, src_e);
      SmallestRec = *src_p;
      return src_e;
    }
    // All 1 SEB.
    SortLargeSeBlock(src_p, src_e);
    SmallestRec = *src_p;
    return src_e;
  }

  ~IsConstBase() {
#ifdef IS_COUNTING
    if (TimingOption <= 1) return;
    std::cout << "IS counts: fwd<= " << Limits[0] << ": " << Counts[0];
    for (int i = 1; i < CntCats; i++)
      std::cout << ", <= " << Limits[i] << ": " << Counts[i];
    std::cout << ", > " << Limits[CntCats-1] << ": " << Counts[CntCats];
    std::cout << '\n';
    //std::cout << "TmMask = " << std::hex << TmMask << std::dec << '\n';
#endif    
  }


  //////////////////////////////////////////////////////////////////////
  //                   Other functions.
// Find functions                                       Find functions  

  /**
   *    Bidirectional find functions to locate an edge of an SeBlock
   *
   * Given an ascending or descending sequence of SeBlocks from src_p
   * to src_e (the search sequence), where TopMask(src_p) !=
   * TopMask(src_e), these functions return a pointer p to the value
   * in the search sequence that is nearest to src_p such that
   * TopMask(p) != TopMask(src_p). By definition, p != src_p, but it
   * may == src_e. Note that src_e is "past the end" of the search
   * sequence. 
   *
   * *** This means that given the seb (lb <= p < ub (all pointers)),
         searching upwards uses ub as src_e, but searching downwards
         needs lb-1 as src_e. Caveat emptor.
   *
   * We don't ever evaluate TopMask(src_e), since it is
   * (logically) != TopMask(src_p) by definition.
   *
   * The calling sequence is similar for both:
   *
   * Template parameter Inc_Val determines the search direction:
   * Inc_Val == 1: search upwards (src_e > src_p)
   * Inc_Val == -1: search downwards (src_e < src_p)
   *
   * Note that Inc_Val is logically redundant, since it can be
   * determined based on whether src_p is > or < src_e, but it allows
   * compile time selection of a conditional in 1 loop, and may allow
   * other optimizations as well.
   *
   * Elem_T *src_p - points into the Seblock in question.
   * Elem_T *src_e - marks the end of the sequence to search.
   *
   * NOTE: These are somewhat slow. Intended for use when a "near"
   * edge isn't found. They don't depend on sentinels.
   **/

  /**
   * Return a pointer p to the first element in the search sequence
   * such that TopMask(p) != TopMask(src_p). Pointer src_e should be
   * "close" for speed, so most calls here should be from FindNextSeb().
   **/
  template<int Inc_Val>
  Elem_T *FindSebBndry(Elem_T *src_p, Elem_T *src_e) {
    Elem_T *p = FindSebBndry_imp<Inc_Val>(src_p, src_e, TopMask(src_p));
#ifdef OCD_TESTING
    if (TopMask(p - Inc_Val) != TopMask(src_p)) ErrTag("FindSebBndry1");
    if (TopMask(p) == TopMask(src_p)) {
      if (p != src_e) ErrTag("FindSebBndry2");
      // else OK
    }
#endif
    return p;
  }
  /**
   * Return a pointer p to the nearest element to src_p in the search
   * sequence such that TopMask(p) != TopMask(src_p). We "improve"
   * src_e before calling FindSebBndry().
   **/
  template<int Inc_Val>
  Elem_T *FindNextSeb(Elem_T *src_p, Elem_T *src_e) {
    Elem_T sebTmv = TopMask(src_p);
    static const int SkipVal = 8;  // Starting inc for search.
    int skip = SkipVal * Inc_Val;
    Elem_T * p = src_p + skip;
    // Searching closer locations first should give best average
    // results in a semisorted vector, and avoids looking at remote
    // memory locations that won't be in cache and won't be needed
    // again soon.
    while ((Inc_Val > 0) ? (p < src_e) : (p > src_e)) {
      if (TopMask(p) != sebTmv) goto found;
      src_p = p;  // Still in initial seb
      skip = 2 * skip;
      p += skip;
    }
    p = src_e;  // Put p back in range.
  found:
    return FindSebBndry<Inc_Val>(src_p, p);
  }
  /**
   * Locate start of last SeBlock which starts before (src_e - Nf). We
   * need Nf+1 slots so that if p is sentinel addr, p[Nf] is still in
   * the sentinel. The top sentinel will be placed there. Sort the
   * part of the vector from the start of the sentinel area. We don't
   * actually put the sentinel in place. Doesn't assume bottom
   * sentinel.
   **/  
  Elem_T *SortTopSentinelArea(Elem_T *src_p, Elem_T *src_e) {
    // Need (Nf + 2) to have (Nf +1) slots before end of valid space.
    Elem_T *sentP = src_e - (Nf + 2);  // Initial sentinel start.
    // If last rec not in same seb as *sentP, we reduce src_e until it is.
    Elem_T *oldEnd = src_e;
    Elem_T sebID = TopMask(sentP);
    while (sebID != TopMask(--src_e)) ;  // *src_e is last rec in seb.
    src_e++;  // Now at end of seb
    // Sort the recs after new src_e (if more than 1)
    if (oldEnd > (src_e + 1)) NoSentIsSort(src_e, oldEnd);
    // Make sentP start of seb.
    sentP = FindNextSeb<-1>(sentP, src_p - 1);
    // Need + 1 unless we are at start of sequence to get addr of
    // space, not away side of edge.
    if (sentP != src_p) sentP++;
    if ((src_e - sentP) <= IsThreshold) {  // small (expected) case
      NoSentIsSort(sentP, src_e);
    } else {  // Bigger seb, and sentinels aren't in place.
      SortLargeSeBlock(sentP, src_e);
    }
    return sentP;
  }
  // Sort a Seblock that is too big for insertion sort, or has no
  // sentinels in place. Let std::sort() figure out if it is really
  // small.  ***** Put in IS sort when appropriate? NO. When called
  // from SortSeBlock() we know its too big for IS. ****
  void SortLargeSeBlock(Elem_T *src_p, Elem_T *src_e) {
    bool DBGtag = false;
    if ((src_e - src_p) < FrThreshold) {
      std::sort(src_p, src_e);
    } else {  // needs a FR recurse.
      SdTemplate<Elem_T> SdObj/*(src_v, src_n)*/;
      ///////////////////////////////////////////////////////////////////
      // Assumes radix digits are bytes in unsigned Elem_Ts (longs)!.
      int nEffectiveDigits = (TmMask == 0) ? (sizeof(TmMask) - 1) : 0;
      if (nEffectiveDigits == 0) {  // Not the whole Elem_T
	// Can do better if we see fld defn!
	unsigned long digitMask = 0xff;
	while ((digitMask & TmMask) == 0) {
	  nEffectiveDigits++;
	  digitMask <<= 8;
	}
      }
      RecursiveLevel++;
      SdObj.DoSort(src_p, src_e-src_p, nEffectiveDigits, TmMask);
      RecursiveLevel--;
      //std::sort(src_p, src_e);
      //DBGtag = true;
      //std::cout << "deleting SortLargeSeBlock::SdObj" << std::endl;
    }
    //if (DBGtag) std::cout << "Scope exited. delete OK" << std::endl;
  }
  /**
   * too_small < p[0], so sift it down. p[1] is empty (was too_small).
   * Bottom sentinel is in place.
   **/
  INLINE_ATT void SiftDown(Elem_T too_small, Elem_T *p) {
    p[1] = p[0];
    while (*--p > too_small) p[1] = *p;
    p[1] = too_small;
  }
  /**
   * too_small < p[0], so sift it down. p[1] is empty (was too_small).
   * No bottom sentinel is in place. bottom_p is lowest addr to look
   * at.
   **/
  INLINE_ATT void SiftDown(Elem_T too_small, Elem_T *p, Elem_T *bottom_p) {
    p[1] = p[0];
    for (p--; (bottom_p <= p) && (*p > too_small); p--) p[1] = *p;
    p[1] = too_small;
  }
  
  /**
   * Simple in place insertion sort with no sentinels.
   **/  
  void NoSentIsSort(Elem_T *src_p, Elem_T *src_e) {
    Elem_T *p = src_p + 1;
    while (p < src_e) {
      Elem_T t = *p;
      if (t < p[-1]) SiftDown(t, p - 1, src_p);
      p++;
    }
  }
private:
  void FindTtlNBkts(unsigned long mask) {
    while ((mask & 255) == 0) mask = mask >> 8;
    if ((mask & 15) == 0) mask >> 4;
    if ((mask & 3) == 0) mask >> 2;
    if ((mask & 1) == 0) mask >> 1;
    TotalNBkts = mask + 1;
    //Acpb = double(TotalN) / TotalNBkts;
  }
#include "IsConstBase.C"
  
};  // IsConstBase

#endif
