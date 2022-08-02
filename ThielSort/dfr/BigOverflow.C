/******************************************************************
 
                        BigOverflow.C

This file is included into BigOverflow.H.

*******************************************************************/
// Compute the amount to add to n to make it 0 mod RtPerMu. Add that
// value to ReqdPadding and return the padded n.
CountType PadN(CountType n) {
  int t = Rt_Per_Mu - (n % Rt_Per_Mu);
  ReqdPadding += t;
  return n + t;
}

void InitBktVec_imp(const CountType *p_bkt_ns, const CountType *of_bkt_ns) {
  TotalPN = TotalOfN = ReqdPadding = 0;
  for (int b = 0; b < N_Buckets; b++) {
    Bkt[b].InitBktClass();
    CountType pN = p_bkt_ns[b];
    if (pN > 0) {
      TotalPN += Bkt[b].PbN = PadN(pN);
    }
    CountType ofN = of_bkt_ns[b];
    if (ofN > 0) TotalOfN += Bkt[b].OfN = PadN(ofN);
    else Bkt[b].OfN = 0;
  }    
}
CountType CreateAreaMaps_imp() {
  // Current free offset in the 2 areas
  CountType POffset = 0;
  CountType OfOffset = 0;
  int b = 0;  // Running index for next 2 loops.
  if (TotalPN > ReqdPadding) {  // Move some pbkts to OF area
    /**
     * Scan through the pbkts and tag them to be moved until we have
     * tagged enough. Tag all ofbkts and the rest of the pbkts to
     * reside in the primary area.
     **/
    // First bucket stays in primary area (so we don't have to move
    // the pbkt), unless we need to use it to have enough.
    if ((TotalPN - Bkt[0].PbN) >= ReqdPadding) {  // We don't need to move [0].
      Bkt[b].SetOffsets(POffset, POffset, 0);
      b++;
    }  // Else tag [0] in loop.
    // First loop: mark all pbkts for OF area til we've marked enough.
    for (CountType padding = ReqdPadding; padding > 0; b++) {
      Bkt[b].SetOffsets(OfOffset, POffset, 1);
      padding -= Bkt[b].PbN;
    }
    // We have marked at least ReqdPadding recs to be moved to OF
    // area. The rest of the chunks will be in the primary area.
  } else { 
    /**
     * Not enough recs in pbkts. So we'll just move enough
     * overflows. By definition, if there aren't very many recs in the
     * primary space, there have to be lots in the overflow area. I
     * won't try to optimize the selection here, since it is hard to
     * imagine it really happening. So this code should work, but may
     * in principle allow a much bigger OF area than is needed.
     **/
    // First loop: mark all ofbkts for OF area til we've marked enough.
    for (CountType padding = ReqdPadding; padding > 0; b++) {
      Bkt[b].SetOffsets(POffset, OfOffset, 2);
      padding -= Bkt[b].OfN;
    }
  }
  // Second loop: mark all chunks for primary area.
  for (; b < N_Buckets; b++) {
    Bkt[b].SetOffsets(POffset, POffset, 0);
  }
  return OfOffset;
}
// Moves to OF area always safe.
void MoveEasyBkts(Elem_T** new_addr_v) {
  for (int b = 0; b < N_Buckets; b++) {
    Bkt[b].EasyMove(new_addr_v[b]);
  }
}
// Move pbkts that must slide left (towards start of space).
void ShiftBktsLeft(Elem_T** new_addr_v) {
  for (int b = 0; b < N_Buckets; b++) {
    Bkt[b].LeftMove(new_addr_v[b]);
  }
}
// Move pbkts that must slide right (towards end of space).
void ShiftBktsRight(Elem_T** new_addr_v) {
  for (int b = N_Buckets - 1; b >= 0; b--) {
    Bkt[b].RightMove(new_addr_v[b]);
  }
}


#if 0
/**
 * Check to see if the target space in the primary area is unused. If
 * so, move the bucket to the remapped location and return
 * 1. Otherwise do nothing and return 0.  Since the target space and
 * the old space are both in the primary area, we can compare them.
 **/
int MoveIfPossible(Elem_T** new_addr_v, int first, int last) {
  nMoved = 0;
  for (int b = first; b <= last; b++) {
    destP = new_addr_v[b];
    destEnd = destP + Bkt[b].PbN;
    if (destP > Bkt[b].OldPbAdr) {  // Sliding right
      int bb = last;
      while (destP < Bkt[bb].OldPbAdr) bb--;
      // Bkt[bb].OldPbAdr <= destP < Bkt[bb+1].OldPbAdr. Since
      // Bkt[last] is the last Bkt that can cause overlap, we don't
      // care about ~[bb+1] if bb is last.
      if (destP >= Bkt[last].OldPbAdr) { 
      int bb = b+1;
    int test = Bkt[b+1].OverlapStatus(destP, destEnd);
    if (test < 0) {  // Moving right.
      
      for (int bb = b + 1; bb <= last; bb++) {
	test = Bkt[bb].OverlapStatus(destP, destEnd);
    nMoved += Bkt[b].(new_addr_v[b]);
  }
  return nMoved;
}
#endif
