/******************************************************************
 
                        HybridFlusher.C

*******************************************************************/
/**
 * Ladle[b] must be flushed to the fenced area in bkt[b]. We may need
 * to overflow into the true OF area.
 **/
void LadleToFencedArea(int b_ndx, int b) {
  SRec_T* ldlAddr = LdlStart[b];
  int len = StoreN.BV[b];
  CountType destNdx = BktNext[b];
  CountType bktEnd = BktEnd[b];
  if ((destNdx + len) <= bktEnd) {  // Flush into fenced area
    BktNext[b] += CopyRecs(b, ldlAddr, BktBuffer + destNdx, len);
  } else {// we have new bkt with true OF
    FenceOfBkts.Delete(b_ndx);
    TrueOfBkts.Append(b);
    // Must flush to some combination of bucket and OF area.
    int leftovers = len % RtPerMu;
    int flushLen = len -leftovers;  // want leftovers in only for some code
    int nToBkt = bktEnd - destNdx;
    if (nToBkt == 0) {  // bkt exactly full, so flush all to OF area.
      SRec_T *destP = OFArea.GetOfBlockP(flushLen);
      OflowN[b] += CopyRecs(b, ldlAddr, destP, len);
    } else if (flushLen == nToBkt) {  // Exactly fills bkt
      BktNext[b] += CopyRecs(b, ldlAddr, BktBuffer + destNdx, len);
    } else {  // Must flush nToBkt to bkt, rest (except leftovers) to OF
      // First fill up rest of fenced area
      MoverT::CopyRecs(ldlAddr, BktBuffer + destNdx, nToBkt);
      BktNext[b] = bktEnd;
      // Move rest to OF area.
      flushLen = flushLen - nToBkt;  // Len for OF area.
      OflowN[b] += flushLen;
      MoverT:: CopyRecs(ldlAddr + nToBkt, OFArea.GetOfBlockP(flushLen),
			flushLen);
      if ((RtPerMu > 1) && leftovers) {
	// Slide any leftovers to start of ladle
	FlushMovers::SlideTail<SRec_T, RtPerMu>(StoreNext[b] - leftovers,
			       ldlAddr, leftovers);
	ldlAddr += leftovers;
      }  //else no leftovers.
      StoreNext[b] = ldlAddr;
    }
  }

}
// NewOverflow                                                  NewOverflow
/**
 * bkt[b] is a new fence overflow. It cannot overflow all the way to a
 * true OF, since that would imply that BatchSize is more than the
 * size of a fenced area. Marking a bkt as overflowing before the
 * first item is dealt doesn't go through here.
 **/
void NewOverflow(int b) {
  FenceOfBkts.Append(b);
  // Init BktNext and move leftovers to front of ladle.
  SRec_T* ldlAddr = StoreStart[b] = LdlStart[b];
  if (RtPerMu > 1) {  // Move any leftovers to start of ladle.
    SRec_T *alignedSn = AlignerT::AlignMuTDown(StoreNext[b]);
    BktNext[b] = alignedSn - BktBuffer;
    int leftovers = StoreNext[b] - alignedSn;
    for (int i = leftovers-1; i >= 0; i--) ldlAddr[i] = alignedSn[i];
    ldlAddr += leftovers;
  } else {
    BktNext[b] = StoreNext[b] - BktBuffer;    
  }
  StoreNext[b] = ldlAddr;
  StoreFence.BV[b] = Trigger;
}
/**
 * Adjust the fence for each ladle in bkt_list to minimum so that
 * after the next flush only the leftovers will remain in the ladle.
 **/
void SetMinFence(const BucketList<MaxNBuckets, false> &bkt_list) {
  for (int bndx = bkt_list.GetN() - 1; bndx >= 0; bndx--) {
    int b = bkt_list.GetItem(bndx);
    StoreFence.BV[b] = RtPerMu - 1;  // Flush if >= RtPerMu recs
  }
}
/**
 * Count the number of true overflows by summing OflowN[b] for OF bkts.
 **/
CountType SumOfCounts() {
  int ttlOverflows = 0;  // Find total overflow count
  for (int bNdx = TrueOfBkts.GetN() - 1; bNdx >= 0; bNdx--) {
    int b = TrueOfBkts.GetItem(bNdx);
    ttlOverflows += OflowN[b];
  }
  return ttlOverflows;
}
#if 0
// Set up to use a variant of this?? Probably!!!!!!!!!!!!!

// MainFlush                                                    MainFlush
/**
 * Called at most once per batch to flush any ladles that need
 * flushing. Designed to support prefetching the ladle and the bkt
 * space ( not using NT instrs), with independent pf leads.
 **/
NO_INLINE void MainFlush() {
  // Now flush targeted ladles. Loop done in 5 pieces, some of which
  // may be empty:
  // 1:{ pf bkt space }, (BktNToPrime - LdlNToPrime) 
  // 2{ pf bkt space, pf ldl space }, LdlNToPrime
  // 3{ pf bkt space, pf ldl space, flush a ladle }, (NBuckets - BktNToPrime)
  // 4:{ pf ldl space, flush a ladle }, (BktNToPrime - LdlNToPrime)
  // 5:{ flush a ladle } LdlNToPrime

  static_assert((BktNToPrime >= LdlNToPrime) || (BktNToPrime <= 0),
		"BktNToPrime < LdlNToPrime in MainFlush()??");

  // Convert *NToPrime to 0 if < 0 so loop ends don't get confused by
  // negatives. The pf is ignored in the pf function when NToPrime < 0.
  // Also, if BktNToPrime <= 0, make BktLeader0 = LdlLeader0.
  constexpr int LdlLeader0 = (LdlNToPrime < 0) ? 0 : LdlNToPrime;
  constexpr int BktLeader0 =
    (BktNToPrime < LdlLeader0) ? LdlLeader0 : BktNToPrime;

  // Now that the static_assert is done, make sure loops end properly
  // on short flush.
  int NFlushing = Flushables.GetN();
  int BktLeader = (BktLeader0 <= NFlushing) ? BktLeader0 : NFlushing;
  int LdlLeader = (LdlLeader0 <= NFlushing) ? LdlLeader0 : NFlushing;

  int Loop1N = BktLeader - LdlLeader;
  int Loop2N = LdlLeader;
  int Loop3N = NFlushing - BktLeader;
  int Loop4N = BktLeader - LdlLeader;
  int Loop5N = LdlLeader;

  // Loop indices for 3 loops through the  data.
  int bktPfNdx = 0;
  int ldlPfNdx = 0;
  int flushNdx = 0;

  for (int i = 0; i < Loop1N; i++) {
    int b = Flushables.GetItem(i);
    FMover.PfBktSpace(b, CurN.IV[b]);
  }
  bktPfNdx += Loop1N;
  for (int i = 0; i < Loop2N; i++) {
    int b = Flushables.GetItem(bktPfNdx + i);
    FMover.PfBktSpace(b, CurN.IV[b]);
    b = Flushables.GetItem(i);
    FMover.PfLadle(b, CurN.IV[b]);
  }
  bktPfNdx += Loop2N; ldlPfNdx += Loop2N;
  for (int i = 0; i < Loop3N; i++) {
    int b = Flushables.GetItem(bktPfNdx + i);
    FMover.PfBktSpace(b, CurN.IV[b]);
    b = Flushables.GetItem(ldlPfNdx + i);
    FMover.PfLadle(b, CurN.IV[b]);
    b = Flushables.GetItem(i);
    FMover.FlushLdl(b, CurN.IV[b], LeftOvers.IV[b]);
  }
  ldlPfNdx += Loop3N; flushNdx += Loop3N;
  for (int i = 0; i < Loop4N; i++) {
    int b = Flushables.GetItem(ldlPfNdx + i);
    FMover.PfLadle(b, CurN.IV[b]);
    b = Flushables.GetItem(flushNdx + i);
    FMover.FlushLdl(b, CurN.IV[b], LeftOvers.IV[b]);
  }
  flushNdx += Loop4N;
  for (int i = 0; i < Loop5N; i++) {
    int b = Flushables.GetItem(flushNdx + i);
    FMover.FlushLdl(b, CurN.IV[b], LeftOvers.IV[b]);
  }
}  // MainFlush
#endif
