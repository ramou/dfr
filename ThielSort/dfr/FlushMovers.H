/* $Header:  $ */ 

#ifndef FlushMoversH
#define FlushMoversH

/********************************************************************
 *
 *                            FlushMovers.H
 *
 * We need to try several ways to move data from the ladles to the
 * buckets, and different methods may work best for different
 * archetectures. I want to keep the code which actually moves data
 * separate from the more complicated code to deal with the ladle
 * status so that the various move classes are small and hold only
 * logic directly required by moving n_e items from ldl_p to
 * block_end_p. This leads to a 3 level design:
 * 
 * - the flusher / LadleManager decides which ladles need to be
 *   flushed and calls functions in class FMoverT to do the flush. It
 *   also handles the addressing manipulation for split moves and
 *   overflows.
 *
 * - class FMoverT functions handle address translation where needed
 *   and call functions in the selected mover class.
 *
 * - the mover class is only concerned with the move.
 *
 * Class FMoverT is defined in this file and the various mover classes
 * are selected with the value in SELECT_FLUSH_MOVER. 
 *
 * This file #includes the small FlushMover classes that contain
 * functions to move data as needed duing the flushing operation. Each
 * one declares a class that moves data of a certain kind (e.g.,
 * interleaved in the ladle or not, and of a certain size) in a
 * particular way (e.g., using simd ops or not). For the moment, we
 * expect that data type to divide the simd type evenly.
 *
 * Much of the configuration is determined by the selection of the
 * FlushMover / ladle manager pair. The amount of interleaving, direct
 * or indirect ladling, and the alignment must all be consistent
 * between the 3 levels.
 *
 * Note that the alignment is now named RtPerMu
 * (RecordTypePerMoveUnit) and is a static constexpr member of the
 * individual FlushMovers. At first glance, the 2 ideas don't seem
 * synonymous. We may want to repect CL alignment even if we are using
 * an MU of half a CL. But if we are respecting CL alignment, if we
 * have 32 byte "MUS", we are only moving them in pairs, and the true
 * MU is really a CL. 
 *
 * We implement this by using NativeRtPerMu to specify the minimum
 * value compatable with the mover (normally 1 or RtPerCL) and
 * ExtendedRtPerMu to be the value supplied by the configuration
 * (Extended_RtPerMu). We choose the value for RtPerMu as follows:
 *
 * - ExtendedRtPerMu <= 1:  RtPerMu = NativeRtPerMu. 
 * - ExtendedRtPerMu = n*NativeRtPerMu: RtPerMu = ExtendedRtPerMu.
 * - Anything else generates a compile time error.
 *
 * Each FlushMover class is templated with the following parameters:
 *
 * SRec_T : The type of the sort record.
 *
 * //Il_N : The ladle interleave value.
 * **NOTE: Now global const IlN, from config constant IL_VALUE.  
 *
 * Class IlElemT = IlElem<SRec_T, IlN> is the ladle array class and
 * means most code can ignore interleaving.
 *
 * Each class contains:
 *
 * - typedef BktMuT : the bulk move unit, overlaid on the bucket space.
 * - static const int RtPerMu : N sort recs in a move unit.
 * 
 * - void FlushMus(IlElemT* ldl_p, BktMuT* bkt_end_p, int n_e);
 *   Flush n_e elements from ldl_p to bkt_end_p. n_e must be a
 *   multiple of RtPerMu. It may be 0. 
 *
 ********************************************************************/

#include "IlElem.H"

namespace FlushMovers {

  // Helper for FindRPMHelper() below
  constexpr int FindRPMHelper(int rec_len, int low_order_0_cnt) {
    return (rec_len & 1) ? low_order_0_cnt :
      FindRPMHelper(rec_len >> 1, low_order_0_cnt + 1);
  }
  /**
   * Given a record length and the simd log len, find log of RtPerMu
   * (the minimum number of records in an Mu such that an Mu will
   * respect the simd boundry). To force cache line alignment, use 6
   * instead of the simd log. But this may give a bigger RtPerMu, and
   * thus more leftovers.
   **/
  constexpr int FindLogRtPerMu(int rec_len, int simd_log) {
    return simd_log - FindRPMHelper(rec_len, 0);
  }

  // Look at native RtPerMu and config ExtendedRtPerMu and compute the
  // effective RtPerMu.
  constexpr int FindRtPerMu(int native_n) {
    return (ExtendedRtPerMu <= native_n) ? native_n :
      ((ExtendedRtPerMu % native_n) == 0) ? ExtendedRtPerMu : -1;
  }
  
  template<typename SRec_T, int Selector_Id, int Il_N = IL_VALUE>
  class MoverClass;

  /**
   * copy n_e < Rt_Per_Mu elements from src_p to dest_p using an
   * unwound loop of max size (Rt_Per_Mu - 1). Callable from movers.
   **/  
  template<typename SRec_T, int Rt_Per_Mu>
  void SlideTail(SRec_T* src_p, SRec_T* dest_p, int n_e) {
    //if (n_e >= RtPerMu) ErrTag(); // Nope
    // for written to allow unwinding (with internal exit). 
    for (int i = 0; i < Rt_Per_Mu; i++) {
      if (i >= n_e) return;
      dest_p[i] = src_p[i];
    }
  }
  
#include "Flush1By1.H"
#include "FlushSimd.H"
  //#include "FastSimd.H"
  
#if 0
  // exclude some variant flushers
#include "FlushParaCL.H"
#include "FlushNT8.H"
#include "FlushBinPair.H"
#endif
  
#include "FlushStatus.H"
#include "MuAligner.H"
  
  // FMoverT                                                       FMoverT
  template<typename SRec_T> class FMoverT {
  public:
    ////////////////////////////////////////////////////////////////
    // Declare mover template chosen from the LegoDefs.H lists
    typedef MoverClass<SRec_T, SELECT_FLUSH_MOVER> MoverT;
    typedef IlElem<SRec_T, IlN> IlElemT;
    typedef typename MoverT::MuT MuT;  // For simd, simd type else SRec_T
    typedef MuAligner<SRec_T, MuT> AlignerT;
    
    static const constexpr int RtPerMu = MoverT::RtPerMu;
    static const int NRecPerPage = SystemPageSize / sizeof(SRec_T);

    typedef FlushStatus<SRec_T, MuT, RtPerMu> FlushStatusT;

    // Make visible to Flushers.
    FlushStatusT FlushDesc[MaxNBuckets];
    
    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor

    //////////////////////////////////////////////////////////////////////
    //                   Other functions.

    static const char *WhoAmI() { return MoverT::WhoAmI(); }
    // readers/writers
    MoverT *GetMoverObj() { return &Mover; }
    void InitLdlAddr(int bkt_id, SRec_T *ldl_p, SRec_T *bkt_p)
    {
      FlushDesc[bkt_id].InitLdlAddr(ldl_p, bkt_p);
      BktPfPos[bkt_id] = 0; 
    }
    // This "init" is really used to store target overflow address.
    void InitBktAvail(int bkt_id, SRec_T *p) {
      FlushDesc[bkt_id].InitBktAvail(p);
    }
    // Remove leftovers. IlN isn't present.
    int FLen(int ttl_len) { return ttl_len & ~(IlN-1); }
    /**
     * Flush ladle ldl0_id and ldl1_id, two consecutive ladles in a
     * bin.  We expect that the lengths ldl?_len include leftovers but
     * don't reflect IlN, and probably won't be the same.
     **/
    void FlushLdlBinPair(int ldl0_id, int ldl1_id) {
      FlushStatusT* const ldlP0 = FlushDesc + ldl0_id;
      FlushStatusT* const ldlP1 = FlushDesc + ldl1_id;
      int l0 = ldlP0->Len();
      int l1 = ldlP1->Len();
      // Try passing MuT* instead of SRec_T* to FlushMuPairs
      static const int RtPerMuT = sizeof(MuT) / sizeof(SRec_T);
      int diff = l0 - l1;
      MuT* bkt0End = ldlP0->MuBktEnd(l0);
      MuT* bkt1End = ldlP1->MuBktEnd(l1);
      // l0p needed for *both* shortens so we have an even numbered address.
      IlElemT* l0p = ldlP0->IlLadle();
      if (diff < 0) {  // Must shorten ldl1 down to l0
#if 0
	ldlP1->ShortenLo(&Mover, l0, -diff, bkt1End);
	bkt1End += diff / RtPerMuT;  // Cause its negative 
#else
        bkt1End = Mover.template Shorten<false>(l0p + l0, bkt1End, -diff);
#endif
      } else if (l1 < l0) {  // Must shorten ldl0 down to l1
#if 0
        ldlP0->ShortHo(&Mover, l1, diff, bkt0End);
        bkt0End -= diff / RtPerMuT;
#else
	bkt0End = Mover.template Shorten<true>(l0p + l1, bkt0End, diff);
#endif
	l0 = l1;  
      }	 // else already have same length
      if (l0) {
        Mover.FlushMuPairs(ldlP0->IlLadle(), bkt0End, bkt1End, l0);
      }
      ldlP0->SlideTail();
      ldlP1->SlideTail();
    }
    /**
     * Mark ladle bkt_id as already flushed.
     **/
    void MarkFlushed(int bkt_id) { FlushDesc[bkt_id].StoreIlLen(0); }

    // Flusher functions: flush flush_n recs from (possibly
    // interleaved) ladle to bucket space using the designated flush
    // mover. 
    /**
     * Flush flush_n items starting at FlushDesc[bkt_id].Ladle() to
     * FlushDesc[bkt_id].Bkt(), inc FlushDesc[bkt_id].Bkt(), and slide
     * leftover_n items to start of ladle. Most flushes done through
     * here.
     **/
    void FlushLdl(int bkt_id, int flush_n, int leftover_n) {
      FlushDesc[bkt_id].FlushToMvr(&Mover, flush_n);
      if (RtPerMu > 1)  // Else can't be any leftovers otherwise.
	SlideTail(bkt_id, FlushDesc[bkt_id].IlLadle() + flush_n, leftover_n);
    }
    void FlushLdl(int bkt_id) {
      FlushDesc[bkt_id].FlushLdl(&Mover);
    }
    /**
     * Flush flush_n items starting at FlushDesc[bkt_id].Ladle() to
     * FlushDesc[bkt_id].Bkt() and inc
     * FlushDesc[bkt_id].Bkt(). Used to flush leading part of ladle
     * when no room for all.
     **/
    IlElemT* FlushLdlHead(int bkt_id, int flush_n) {
      if (__builtin_expect(flush_n > 0, true)) 
	FlushDesc[bkt_id].FlushToMvr(&Mover, flush_n);
      return FlushDesc[bkt_id].IlLadle() + flush_n;
    }

    /**
     * Move functions can be used directly by flusher for newly
     * overflowing buckets, and to move overflow blocks to final
     * overflow area.
     **/
    IlElemT*FlushLdlTail(IlElemT *ldl_p, SRec_T* block_end_p, int flush_n) {
      Mover.FlushMus (ldl_p, ldl_p + flush_n,
		      reinterpret_cast<MuT *>(block_end_p), flush_n);
      return ldl_p + flush_n;
    }

    // Copy n_e = leftovers (< 1 MU) from end of ladle to target addr,
    // which is not incremented. Target is ...

    // ... start of ladle[bkt_id]. 
#if 0
    void SlideTail(int bkt_id, IlElemT* src_p, int n_e) {
      if (RtPerMu == 1) return;  // Never any tail to flush
      // NOTE: Not current!!
      SRec_T *destP = FlushDesc[bkt_id].Ladle();
      for (int i = 0; i < n_e; i+=IlN) destP[i] = src_p[i];
    }
#else
    void SlideTail(int bkt_id, IlElemT* src_p, int n_e) {
      if (RtPerMu == 1) return;  // Never any tail to flush
      if (n_e == 0) return;
      //if (n_e >= RtPerMu) ErrTag(); Nope
      IlElemT *destP = FlushDesc[bkt_id].IlLadle();
      // for written to allow unwinding (with internal exit). Does it???
      for (int i = 0; i < RtPerMu; i++) {
	if (i >= n_e) return;
	destP[i] = src_p[i];
      }
    }
#endif    
    // ... to bucket space during Close(). Note that during Close(),
    // ladle holds just leftovers, so the "end of the ladle" is really
    // at the start of the ladle.
    void FlushTail(int bkt_id, int n_e) {
      if (RtPerMu == 1) return;  // Never any tail to flush
      const IlElemT *srcP = FlushDesc[bkt_id].IlLadle();
      SRec_T *destP = FlushDesc[bkt_id].Bkt();
      for (int i = 0; i < n_e; i++) {
	//if (srcP[i] == 0) ErrTag("0 found in FlushTail");
	*destP++ = srcP[i];
      }
    }
    ///////////////// Prefetch helpers
    /**
     * Prefetch space to hold n SRec_Ts being flushed to bucket[bkt_id].
     **/
    void PfBktSpace(int bkt_id, int n) {
      // Avoid clutter in the calling loop by testing "don't do it"
      // condition here.
      if ((BktNToPrime < 0) && (SystemPageSize < 0)) return;
      // Following 2 stmts should go away if unused.
      const SRec_T *p = FlushDesc[bkt_id].Bkt();
      const SRec_T *pPlusN = p + n;  
      if (BktNToPrime >= 0) {
	RevPrefetchVec<1, BktLocality>(p, pPlusN - 1);
      }
      if (SystemPageSize > 0) {
	// Nope. Neither at home or on titan
	int curPfPos = BktPfPos[bkt_id] - n;
	if (__builtin_expect(curPfPos < 0, false)) {  // Time for another pf.
	  __builtin_prefetch(pPlusN + NRecPerPage, 0, 0);
	  curPfPos += NRecPerPage;
	}
	BktPfPos[bkt_id] = curPfPos;
      }
    }
    void PfBktSpace(int bkt_id) { FlushDesc[bkt_id].PfBktSpace(); }
    /**
     * Prefetch space holding the n SRec_Ts being flushed from
     * ladle[bkt_id].
     *
     * If IlN > 1, we should avoid multiple fetches of the same parts
     * of a bin. In principle we want to prefetch the biggest n for
     * the ladles in the bin, only once for each bin. For now, just
     * prefetch the first ladle in each bin, and ignore the reset.
     **/
    void PfLadle(int bkt_id, int n) {
      // Avoid clutter in the calling loop by testing "don't do it"
      // condition here.
      if (LdlNToPrime < 0) return;
      // For interleaved case, don't do (mostly) redundant prefetch. If
      // using DirectBinFlusher, following test is wrong, since flusher
      // handles this already, and test is redundant at best, will
      // skip useful pf at worst.
      
      FlushDesc[bkt_id].PfLadle(n);
    }  // PfFTarget

    //
  private:
    MoverT Mover;
    // Used to space out prefetches to force page reads in bkts.
    int BktPfPos[MaxNBuckets];
    
  };  // FMoverT

}; // FlushMovers

#endif
#if 0
////////////// Copied from PfBktSpace() to remove clutter. Comments
////////////// worth saving for reference. Am trying diff approach.
if (SystemPageSize > 0) {
  // Do prefetch from next page in this bucket to make page
  // available for later flush ops. Doesn't apply to this
  // one. The prefetch is useless except as a non-blocking way
  // to force the page to be read if needed. This is effective
  // (my desttop, anyway) on first pass of first run of 3, where
  // activating the allocated memory seems to be costly. But it
  // costs on other flushes. Pays for itself for the first,
  // anomolous, run. Just cost in later runs. Would it help to
  // save the prefetch locations for each bucket and only issue
  // prefetch once per page instead of every flush?? Should
  // help. Worth a try. Later.
  // OPTIMIZE??
  __builtin_prefetch(pPlusN + (SystemPageSize / sizeof(SRec_T)),
		     0, 0);
 }

#endif
