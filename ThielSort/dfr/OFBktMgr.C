// Imp without using simd.

int NBkts;  // Number of bkts left to find in current stripe

// Extract stripe stripe index from ndx.
static int StripeNdx(int ndx) { return ndx & (NPerSimd - 1); }
  
// Return index of next non-empty stripe >= s or NPerSimd if no more. 0 <= s
// < NPerSimd on entry. s <= return value <= NPerSimd on exit.
//////////////////////// Search needs more speed.
int NextStripe(int s) {
  for (; StripeSums.BV[s] == 0; s++) ;  // Halt at sentinel
  if (s < NPerSimd) {
    NBkts = StripeSums.BV[s];
    StripeSums.BV[s] = 0;  // Don't refind this stripe.
  }  // else don't kill sentinel!
  return s; 
}
// NBkts > 0, so we have a live bkt in current stripe. Find it, unmark
// it, and return it.
int NxtInStripe(int ndx) {
  while (TagVec.BV[ndx] == 0) ndx += NPerSimd;
  NBkts--;
  TagVec.BV[ndx] = 0;
  return ndx;
}
// Return bkt id of first OF bkt or -1 if none. May not be smallest
// id!
int ResetList_imp() {
  int ndx = NextStripe(0);
  if (ndx >= NPerSimd) return -1;  // Empty!
  return NxtInStripe(ndx);
}

int ExtractNext_imp(int ndx) {
  ndx += NPerSimd;
  if (NBkts <= 0) {  // No more in this stripe.
    if ((ndx = NextStripe(StripeNdx(ndx))) >= NPerSimd) return -1;
  }
  return NxtInStripe(ndx);
}
using CtDiffT = unsigned CountType;
CountType MakeTag(CountType fence_val, CountType n) {
  CtDiffT diff = fence_val - n;
  return diff >> ShiftCnt;
}
CountType MakeTag(CtDiffT avail) { return avail >> ShiftCnt; }

int MarkOfBkts_imp(PerBktIntT *store_n, PerBktIntT *store_fence) {
  // Sign bit of (fence - next) == 1 iff we have overflow.
  // Sum of sign bits is count of OF bkts.
  for (int b = 0; b < NPerSimd; b++)
    StripeSums.BV[b] = TagVec.BV[b] =
      MakeTag(store_fence->BV[b],store_n->BV[b]);
  for (int b = NPerSimd; b < BSize; b++) {
    StripeSums.BV[StripeNdx(b)] += TagVec.BV[b] =
      MakeTag(store_fence->BV[b],store_n->BV[b]);
  }
  int nmarked = StripeSums.BV[0];
  for (int s = 1; s < NPerSimd; s++) nmarked += StripeSums.BV[s];
  return nmarked;
}
int MarkOfBkts_simd_imp(PerBktIntT *store_n, PerBktIntT *store_fence) {

  SimdT n1 = store_n->SdV[0];
  SimdT f1 = store_fence->SdV[0];
  SimdT LOBits = 1;
  SimdT t1 = f1 - n1;
  n1 = store_n->SdV[1];
  f1 = store_fence->SdV[1];
  SimdT sum1 = (t1 >> ShiftCnt) & LOBits;
  TagVec.SdV[0] = sum1;
  
  for (int i = 1; i < SdSize - 1; i++) {
    t1 = f1 - n1;
    n1 = store_n->SdV[i+1];
    f1 = store_fence->SdV[i + 1];
    t1 = (t1 >> ShiftCnt) & LOBits;
    sum1 = sum1 + t1;
    TagVec.SdV[i] = t1;
  }
  // Last n1, f1 still to be used
  t1 = f1 - n1;
  t1 = (t1 >> ShiftCnt) & LOBits;
  StripeSums.SdV[0] = sum1 + t1;
  TagVec.SdV[SdSize-1] = t1;
  int nmarked = StripeSums.BV[0];
  for (int s = 1; s < NPerSimd; s++) nmarked += StripeSums.BV[s];
  return nmarked;
}
  
int MarkOfBkts_imp(PerBktIntT *n_avail) {
  // Sign bit of n_avail == 1 iff we have overflow.
  // Sum of sign bits is count of OF bkts.
  for (int b = 0; b < NPerSimd; b++)
    StripeSums.BV[b] = TagVec.BV[b] = MakeTag(n_avail->BV[b]);
  for (int b = NPerSimd; b < BSize; b++) {
    StripeSums.BV[StripeNdx(b)] += TagVec.BV[b] =
      MakeTag(n_avail->BV[b]);
  }
  int nmarked = StripeSums.BV[0];
  for (int s = 1; s < NPerSimd; s++) nmarked += StripeSums.BV[s];
  return nmarked;
}
