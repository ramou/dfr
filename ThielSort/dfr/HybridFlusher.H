/* $Header:  $ */ 

#ifndef HybridFlusherH
#define HybridFlusherH

/********************************************************************
 *
 *                            HybridFlusher.H
 *
 * This "flusher" only flushes overflowing data to the overflow
 * area. The dealer deals most recs directly to the target
 * vector. Each bucket holds room for BucketLen[b] recs, including the
 * padding of BatchSize recs which protects against overflow within a
 * batch. I refer to position (BucketLen[b] - BatchSize) as a fence,
 * and the space between the fence and the end of the bkt space as the
 * fenced area.
 *
 * When the dealer drops data past the fence, we change the target
 * addr to the ladle, which works like DirectKissFlusher. Data in the
 * ladle is flushed to either the target bucket (if there is room in
 * the fenced area) or the true overflow area. Only the data flushed
 * to the overflow area needs to be redealt
 *
 * This flusher requires IL_VALUE == 1, since we don't want to
 * interleave the data in the target buffer. We use wide indices to
 * permit deal targets to be in separate vectors (the target vector or
 * the ladle). So we require WIDE_INDICES == true.
 *
 * See discussion in HybridFlusher.Notes for details.
 ********************************************************************/

#include "OFBktMgr.H"

// HybridFlusher                                    HybridFlusher

template<typename SRec_T>
class FlusherClass<SRec_T, DHYBRID_FLUSHER> {
  // Very little needed from LadleManager. Perhaps some prefetch stuff.

  typedef IlElemTemplate<SRec_T> IlElemT;
  typedef SRec_T BucketT;
  typedef FlushMovers::FMoverT<SRec_T> FlushMoverT;
  typedef typename FlushMoverT::MoverT MoverT;
  using AlignerT = typename FlushMoverT::AlignerT;
  typedef DataPoolTemplate<SRec_T> DataPoolType;

  static constexpr int RtPerMu = FlushMoverT::RtPerMu;
  // Estimated bkt sizes <= Tiny are reset to 0.
  static constexpr int Tiny = (Trigger >= RtPerMu) ? Trigger : RtPerMu;

  using PerBktIntT = OFBktMgr::PerBktIntT;

public:

  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  FlusherClass(CountType src_n, TimeWrapper *tw_p = nullptr) :
    SortN(src_n),
    OFArea(src_n, FMover.GetMoverObj()), TimeWpr(tw_p), BktFinder{} {
    
    if (tw_p) {
      CloseTimerID = tw_p->AllocateFnTimer(", (sub)Sort Close time = ");
    }
  }
  ~FlusherClass() {}

  // Called to specify whether the next pass(s) is counted (i.e., has
  // exact bkt sized specified).
  void SetCountedPass(bool v) { CountedPass = v; }
  
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.
  static const char *WhoAmI() { return "HybridFlusher"; }
  static const char *WhoIsLM() { return "none"; }

  /**
   * Init the ladles for this pass. Parameters:
   * dc_p : pointer to dealer core for this pass.
   * batch_size : Size of a batch, and thus the size of the needed overflow
   *              buffer in each bucket. Allows us to place the fence.
   *              If batch_size == 0, this is a counted deal, so there can
   *              be no overflow, and the dealt data will exactly fill the
   *              buckets.
   * safe_ofa : addr of source vector. Used as overflow of sufficient
   *            size (after preface area is filled).
   * flush_tmr : Allows separate timing of flush step.
   * dp_p :      Target DataPool. Not currently used here.
   *
   **/
  template<typename Dealer_Core>
  void InitFlusher(Dealer_Core *dc_p, int batch_size, SRec_T *safe_ofa,
		   int n_buckets, StopWatch *flush_tmr,
		   DataPoolType *dp_p) {
    NBuckets = n_buckets;
    FlushTimer = flush_tmr;
    BatchSize = batch_size;
    TrueOfBkts.Reset();
    FenceOfBkts.Reset();
    OFArea.OverflowInit(safe_ofa, batch_size);
  }  
// InitBkts                                               InitBkts  
  /**
   * Allocate buckets in bkt_buffer[] for the next pass according to
   * bkt_lens[]. parameters:
   *
   * dc_p : pointer to dealer core for this pass. 
   *
   * bkt_buffer : the buffer to use. It must be long enought to hold
   *    the set of buckets defined in bkt_lens[].
   *
   * bkt_lens : bkt_lens[b] contains the requested length in
   *   bkt_buffer of bucket[b]. If CountedPass the requested length is
   *   used. Otherwise, actual bkt sizes in bkt_buffer are adjusted:
   *
   * - change tiny bkts to have length 0.
   * - round each length down to x * RtPerMu
   * - round bkt_buffer up to an aligned addr (and subtract RtPerMu
   *   from the first non-tiny bkt).
   *
   **/
  template<typename Dealer_Core>
  void InitBkts(Dealer_Core *dc_p, BucketT *bkt_buffer,
		CountType *bkt_lens) {
    if (CountedPass) {
      // bkt_lens are exact, and we respect them. 
      BktBuffer = bkt_buffer;
      InitCountedBkts(bkt_lens);
    } else {  // Estimated bkt sizes.

      // Prepare ladle buffer and set up to init LdlStart[].
      // A single ladle can hold at least an entire batch and is
      // rounded up to n*CacheLineLen long.
      int oneLadleLen = (Trigger + BatchSize + RtPerMu + CacheLineLen - 1);
      oneLadleLen -= oneLadleLen % CacheLineLen;
      AlignLadleSpace.GetAdjustedAdr(&LadleBuffer, oneLadleLen * NBuckets);
      int ldlOffset = 0;
      
      // Prepare bkt buffer and set up to init bkt info.
      BktBuffer = AdjustBktLens(bkt_buffer, bkt_lens);
      CountType bktOffset = 0;

      // Init ladles and bkts.
      for (int b = 0; b < NBuckets; b++) {
	LdlStart[b] = LadleBuffer + ldlOffset;
	ldlOffset += oneLadleLen;
	bktOffset = InitEstimatedBktB(b, bkt_lens[b], bktOffset);
      }
    }
    dc_p->InitDc(&StoreNext[0], nullptr, nullptr);
    //BC::PfBin.PfRebase();
  }
  /**
   * Flush SRec_Ts from ladle to bucket or OF.
   **/  
  void Flush() {
    if (CountedPass) return;
    FlushTimer->Start();
    OFArea.CheckOfArea();  // don't always need it, but...
    for (int b = 0; b < NBuckets; b++)
      StoreN.BV[b] = StoreNext[b] - StoreStart[b];
    if (BktFinder.MarkOfBkts(&StoreN, &StoreFence)) {
      FlushTrueOF();
      FlushFenceOF();
      // Any remaining bkts with the OF condition are newly detected
      // fence OF bkts.
      for (int b = BktFinder.ExtractNext(); b >= 0;
         b = BktFinder.ExtractNext(b)) {
	NewOverflow(b);
      }
    }
    PfCurStorePoint();
    FlushTimer->Stop();
  }
// Close                                                               Close
  /**
   * Flush all SRec_Ts from ladles.
   **/
  template<typename DBA, typename Pool_T>void Close(Pool_T *pool_p) {
    TimeWpr->ReStartTimer(CloseTimerID);
    // The dealing is done, so the recs in the current merged overflow
    // buffer have been dealt. Make the space available for the next
    // pass.
    MergedOfMem.MarkFree();
    FlushAll();  // Flush everything, including leftovers
    // Create primary buffer info.
    CountType BktN[MaxNBuckets];
    SRec_T *BktAddr[MaxNBuckets];
    for (int b = 0; b < NBuckets; b++) {
      BktAddr[b] = BktBuffer + BktStart[b];
      BktN[b] = BktNext[b] - BktStart[b];
    }
#if 0
    ConfirmBktsAtClose<DBA>(BktAddr, BktN);
#endif
    if (CountedPass) {  // Data dealt to bkt vector. No flush, no OF.
      for (int b = 0; b < NBuckets; b++) {
	pool_p->DefineChunk(b, BktBuffer + BktStart[b],
			    BktEnd[b] - BktStart[b], nullptr, 0);
      }
    } else {  // Now deal with overflows.
    
      CountType nOverflowing = SumOfCounts();
      OFArea.template CreateChunks<DBA>
	(TrueOfBkts.GetList(), TrueOfBkts.GetN(), BktAddr,  BktN,
	 OflowN, nOverflowing, RtPerMu, pool_p);
    }
    /**
     * Now the Close() is done. That means that the buffer space for
     * the ladles and the preface overflow buffers can be unlocked.
     **/
    AlignLadleSpace.MarkFree();
    PrefOfMem.MarkFree();
    TimeWpr->StopTimer(CloseTimerID);
  }
  

private:
  template<typename DBA>
  void ConfirmBktsAtClose(SRec_T **BktAddr, CountType *BktN) {
    CountType totalBktN = 0;  // DEBUG
    for (int b = 0; b < NBuckets; b++) {
      if (BktN[b] < 0) ErrTag("negative bkt n??");  // DEBUG
      if (int nerrs = ConfirmChunk<DBA>(BktAddr[b], BktN[b], b)) {
        std::cout << "In Close(): " << nerrs << " bad deal bytes in bkt[" <<
          b << "] out of " << BktN[b] << '\n';
      }
      totalBktN += BktN[b];  // DEBUG
    }
    std::cout << "In Close(): " << totalBktN << " in bkts";
    if (CountedPass) {
      std::cout << " for counted pass. ";
    } else {
      CountType nOverflowing = SumOfCounts();
      totalBktN += nOverflowing;
      std::cout << " and " << nOverflowing << " overflowing. Total = " <<
        totalBktN;
    }
    if (totalBktN == SortN) std::cout << " OK\n";
    else std::cout << "Expected " << SortN << '\n';

  }
  /**
   * Confirm that all n items starting at chunk_v have bkt_id in their
   * DBA field. Return number of errors found.
   **/
  template<typename DBA>
  int ConfirmChunk(SRec_T *chunk_v, int n, int bkt_id) {
    int nErrs = 0;
    for (int i = 0; i < n; i++) {
      if (bkt_id != DBA::GetDealB(chunk_v + i)) nErrs++;
    }
    return nErrs;
  }
  /**
   * bkt_lens[] are estimates. Round bkt_buffer up to an MU boundry
   * and reduce an estimate enough to ensure we won't overflow the
   * buffer.  Return the possibly changed bkt_buffer.
   **/
  SRec_T *AdjustBktLens(SRec_T *bkt_buffer, CountType *bkt_lens) {
    SRec_T *bktBfr = AlignerT::AlignMuTUp(bkt_buffer);  // Align buffer.
    if (bktBfr == bkt_buffer) return bkt_buffer;  // Already was aligned.
    // Shorten an estimated bkt size enough to make up for the shorter bfr.
    int n = bktBfr - bkt_buffer;
    for (int b = 0; b < NBuckets; b++) {
      if (bkt_lens[b] >= n) {
	bkt_lens[b] -= n;
	return bktBfr;
      }
    }
    ErrTag("No estimates bigger than alignment adjustment????", true);
    return bktBfr;  // Keep compiler happy. ErrTag() will cause abort.
  }
  /**
   * At end of flush, prefetch current store point for each ladle.  
   **/
  void PfCurStorePoint() {
    static const int UnwindN = 8;
    SRec_T **p = StoreNext;
    for (int b = 0; b < NBuckets; b += UnwindN) {
      for (int i = 0; i < UnwindN; i++) {
	__builtin_prefetch(p[i], 1, 3);
      }
      p += UnwindN;
    }
  }
// FlushAll                                                       FlushAll
  /**
   * End of the pass. Flush everything except leftovers by adjusting
   * StoreFence to trigger flushing if > (RtPerMu - 1) recs are in the
   * ladle. Then flush any leftovers.
   **/
  void FlushAll() {
    // Reset fences for ladles to trigger flush unless only leftovers.
    SetMinFence(TrueOfBkts);
    SetMinFence(FenceOfBkts);
    //CountDealt();  // DEBUG
    Flush();
    // Flush leftovers, if any.
    if (RtPerMu > 1) {
      // Check fence OF bkts for leftovers
      for (int bndx = FenceOfBkts.GetN() - 1; bndx >= 0; bndx--) {
	int b = FenceOfBkts.GetItem(bndx);
	if (int len = StoreNext[b] - StoreStart[b]) {  // Something to flush
	  if (BktNext[b] == BktEnd[b]) {  // bkt is exactly full: make trueOF
	    FenceOfBkts.Delete(bndx);
	    TrueOfBkts.Append(b);
	    // Actual flush done in next loop
	  } else {  
	    SRec_T *p = BktBuffer + BktNext[b];
	    for (int i = 0; i < len; i++) p[i] = StoreStart[b][i];
	    BktNext[b] += len;  // Include leftovers
	    StoreNext[b] = StoreStart[b];
	  }
	}
      }
    //CountDealt();  // DEBUG
      // Now check true OF bkts
      for (int bndx = TrueOfBkts.GetN() - 1; bndx >= 0; bndx--) {
        int b = TrueOfBkts.GetItem(bndx);
        if (int len = StoreNext[b] - StoreStart[b]) {
          SRec_T *p = OFArea.GetOfBlockP(RtPerMu);
	  for (int i = 0; i < len; i++) p[i] = StoreStart[b][i];
	  OflowN[b] += len;  // count leftovers
	  StoreNext[b] = StoreStart[b];
        }
      }
    }
    //CountDealt();  // DEBUG
    // Set BktNext[b] for bkts that haven't overflowed at all.
    for (int b = 0; b < NBuckets; b++) {
      if (StoreNext[b] != StoreStart[b]) {  
	// since all ladles are emptied in above loops, we're dealing
	// directly to this bkt, so haven't crossed the fence yet. so
	// BktNext hasn't been set.
	BktNext[b] = StoreNext[b] - BktBuffer;
      }
    }
    //CountDealt();  // DEBUG
  }
  // All data has been flushed from ladles. Count what's in the bkts
  // and the overflow and print it out. Should be SortN.
  void CountDealt() {  // DEBUG only
    CountType n = 0;
    CountType nOF = 0;
    for (int b = 0; b < NBuckets; b++) {
      if (StoreStart[b] == LdlStart[b]) {
	// We are storing to the ladle. Do we have fence OF or true OF?
	if (OflowN[b])  nOF += OflowN[b];  // True OF found.
	n += (BktNext[b] - BktStart[b]) + (StoreNext[b] - StoreStart[b]);
      } else {  // Still storing direct to bkt.
	n += (StoreNext[b] - StoreStart[b]);
      }
    }
    CountType total = n + nOF;
    std::cout << "FlushAll():Counted " << n << " recs dealt to bkts, " <<
      nOF << " to overflow.\n";
    std::cout << "Total = " << total;
    if (total == SortN) std::cout << ": OK\n";
    else std::cout << ": ***Error! Expected =" << SortN << '\n';
  }
  // InitEstimatedBktB                                     InitEstimatedBktB
  /**
   * Define bkt[b] based on estimated length len. Adjust len as
   * needed.  bkt_offset is the aligned start of this bkt. Return the
   * end of this bkt (also aligned).
   **/
  CountType InitEstimatedBktB(int b, int len, CountType bkt_offset) {
    if (len < Tiny) len = 0;
    else if (RtPerMu > 1) {
      len -= len % RtPerMu;  // Make it a multiple of RtPerMu
    }
    StoreStart[b] = StoreNext[b] = BktBuffer + bkt_offset;
    BktStart[b] = BktNext[b] = bkt_offset;
    BktEnd[b] = bkt_offset = bkt_offset + len;
    if (BatchSize < len) StoreFence.BV[b] = len - BatchSize;
    else StartFenced(b, len);
    OflowN[b] = 0;
    return bkt_offset;
  }
  /**
   * Bkt[b] is so small no recs can be dealt directly to it, so it
   * starts out fenced.
   **/
  void StartFenced(int b, int len) {
    StoreStart[b] = StoreNext[b] = LdlStart[b];
    StoreFence.BV[b] = Trigger;
    if (len == 0) {  // True OF
      TrueOfBkts.Append(b);
    } else {
      FenceOfBkts.Append(b);
    }
  }
  // InitCountedBkts                                          InitCountedBkts
  /**
   * With counted bkt sizes we ignore alignment (since dealing is 1 by
   * 1) and use exactly what we are given. No OF, no problems. Not
   * much context needed: StoreNext, BktStart, and BktEnd
   **/
  void InitCountedBkts(CountType *bkt_lens) {
    SRec_T *bktP = BktBuffer;
    CountType bktOffset = 0;
    for (int b = 0; b < NBuckets; b++) {
      BktStart[b] = bktOffset;
      StoreNext[b] = bktP + bktOffset;
      BktEnd[b] = bktOffset = bktOffset + bkt_lens[b];
      // StoreFence.BV[b] not used.
      OflowN[b] = 0;  // just for safety. Unused.
    }
  }
  // CopyRecs                                                  CopyRecs
  /**
   * Flush (n_e - leftovers) recs from ldl_p to dest_p. If there are
   * any leftovers, slide them to ldl_p.  
   * *** Reset StoreNext[b] = ldl_p + leftovers.
   *
   * ldl_p : addr of this ladle
   * dest_p : addr to flush to
   * n_e : number or recs in ladle, including any leftovers.
   * return n actually flushed to dest_p (not counting leftovers).
   **/
  int CopyRecs(int b, SRec_T *ldl_p, SRec_T *dest_p, int n_e) {
    // MoverT::CopyRecs ignores leftovers.
    MoverT::CopyRecs(ldl_p, dest_p, n_e);
    if (RtPerMu > 1) {  // Leftovers possible
      if (int leftOvers = n_e % RtPerMu) {  // leftovers exist
	// Slide leftovers to start of ladle
	FlushMovers::SlideTail<SRec_T, RtPerMu>
	  (ldl_p + (n_e - leftOvers), ldl_p, leftOvers);
	StoreNext[b] = ldl_p + leftOvers;
	return n_e - leftOvers;
      }
    }  // Else no leftovers.
    StoreNext[b] = ldl_p;
    return n_e;
  }  // CopyRecs
  // FlushTrueOF                                  FlushTrueOF
  /**
   * Scan the list of ladles that have been tagged to be flushed to
   * the true overflow area. For each ladle b in the list where
   * StoreNext[b] is past StoreFence.BV[b], move recs from the ladle
   * into OFArea.  
   **/
  void FlushTrueOF() {
    for (int bNdx = TrueOfBkts.GetN() - 1; bNdx >= 0; bNdx--) {
      int b = TrueOfBkts.GetItem(bNdx);
      if (BktFinder.TestAndClear(b)) {  // must flush this one.
	int len = StoreN.BV[b]; 
	int leftovers = (RtPerMu > 1) ? len % RtPerMu : 0;
	SRec_T *destP = OFArea.GetOfBlockP(len - leftovers);
	OflowN[b] += CopyRecs(b, LdlStart[b], destP, len);
      }
    }
  }  // FlushTrueOF
  // FlushFenceOF                                  FlushFenceOF
  /**
   * Scan the list of ladles that have been tagged to be flushed to
   * the fenced area. For each ladle b in the list where
   * StoreNext[b] > StoreFence.BV[b], move recs from the ladle to
   * the fenced area. If fenced area overflows, split the flush and
   * move the bktID from the fenced OF list to the true OF list.
   **/
  void FlushFenceOF() {
    static_assert(IlN == 1);  // For this function, at least.
    for (int bNdx = FenceOfBkts.GetN() - 1; bNdx >= 0; bNdx--) {
      int b = FenceOfBkts.GetItem(bNdx);
      if (BktFinder.TestAndClear(b)) LadleToFencedArea(bNdx, b);
    }
  }  // FlushFenceOF
  // See HybridFlusher.Notes for more detail.

  //                       Count vectors

  /**
   * Describe the bkts. BktBuffer and each individual bkt is MU
   * aligned unless this is a counted pass. BktBuffer, BktStart, and
   * BktEnd are const during a pass. BktNext[b] starts out ==
   * BktStart[b] but is updated if needed at the end of each batch to
   * point to the next free slot. If the bkt is full, BktNext[b] ==
   * BktEnd[b].
   **/
  SRec_T *BktBuffer = nullptr;  
  CountType BktStart[MaxNBuckets];
  CountType BktNext[MaxNBuckets];
  CountType BktEnd[MaxNBuckets];

  /**
   * Pointers to the ladles. LadleBuffer and each individual ladle is
   * MU aligned, and is const during a pass.
   **/
  SRec_T *LadleBuffer = nullptr;
  SRec_T *LdlStart[MaxNBuckets];

  /**
   * Describe the store containers. Store[b] is tied to either
   * ladle[b] or bkt[b]. Values are MU aligned iff the base vector
   * is. These values change as values are dealt.
   **/
  // Addr to deal the next record to. Points into a bkt or ladle.
  SRec_T *StoreNext[MaxNBuckets];
  // Addr of start of the bkt or ladle where StoreNext[b] is pointing.
  // Allows us to compute the n within the active store area.
  SRec_T *StoreStart[MaxNBuckets];
  // StoreN.BV[b] = StoreNext[b] - StoreStart[b]. Set at the end of
  // each pass. Only reflect previous pass' deal.
  PerBktIntT StoreN;
  // When StoreN.BV[b] > StoreFence.BV[b], we have (fence or true) OF.
  PerBktIntT StoreFence;
  
  
  OFBktMgr BktFinder;
  
  // OflowN keeps running total of recs dealt to OF area for each bkt.
  // Updated at end of each batch.
  int OflowN[MaxNBuckets];

  // Lists of buckets by OF status
  BucketList<MaxNBuckets, false> TrueOfBkts;  // True OF
  BucketList<MaxNBuckets, false> FenceOfBkts;  // Flush ladles to fenced area

  Overflow<SRec_T, MaxNBuckets> OFArea;
  StopWatch *FlushTimer;

  // Code to actually move stuff lives here.
  FlushMoverT FMover;

  TimeWrapper *TimeWpr = nullptr;
  int NBuckets;
  int BatchSize;
  int CloseTimerID = -1;
  CountType SortN = 0;
  bool CountedPass = false;  // Set as needed from SetCountedPass().
  
#include "HybridFlusher.C"
  // DirectFlushCommon.C not compatible!
  //#include "DirectFlushCommon.C"
  
};  // HybridFlusher

#endif
